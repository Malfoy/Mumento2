//! External sorting for minimizer tuples
//!
//! This module implements RAM-bounded external sorting for large datasets,
//! following the C++ sshash implementation precisely:
//!
//! 1. Each thread has a RAM-bounded buffer
//! 2. When buffer fills → parallel sort + flush to temp binary file
//! 3. After all tuples processed → k-way merge of temp files
//! 4. Merge uses buffered file I/O (not mmap) to avoid resident page bloat
//!
//! ## Memory Management
//!
//! Buffer size per thread = `(ram_limit_gib * GiB) / (2 * sizeof(tuple) * num_threads)`
//!
//! The factor of 2 accounts for temporary memory during parallel sort.
//!
//! After merge, the sorted file is accessed via sequential buffered readers
//! (`FileTuples`) rather than mmap, keeping RSS proportional to the I/O buffer
//! size (~4 MB) instead of the file size (~8-10 GB).

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use memmap2::Mmap;
use rayon::prelude::*;
use tracing::{debug, info};

use super::minimizer_tuples::MinimizerTuple;

/// Size of `MinimizerTupleExternal` in bytes (packed, no padding)
pub const TUPLE_SIZE_BYTES: usize = 18;

/// Bytes per GiB
pub const GIB: usize = 1024 * 1024 * 1024;

/// Default buffer size for sequential tuple readers (4 MB).
const READER_BUF_SIZE: usize = 4 * 1024 * 1024;

/// Packed minimizer tuple for disk I/O (matches C++ `#pragma pack(push, 2)`)
///
/// Layout: minimizer (8) + pos_in_seq (8) + pos_in_kmer (1) + num_kmers_in_super_kmer (1) = 18 bytes
#[repr(C, packed(2))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MinimizerTupleExternal {
    /// The minimizer hash value
    pub minimizer: u64,
    /// Position in the sequence
    pub pos_in_seq: u64,
    /// Position of the minimizer within the k-mer
    pub pos_in_kmer: u8,
    /// Number of k-mers in the super-k-mer
    pub num_kmers_in_super_kmer: u8,
}

impl MinimizerTupleExternal {
    /// Convert from internal MinimizerTuple
    pub fn from_internal(t: &MinimizerTuple) -> Self {
        Self {
            minimizer: t.minimizer,
            pos_in_seq: t.pos_in_seq,
            pos_in_kmer: t.pos_in_kmer,
            num_kmers_in_super_kmer: t.num_kmers_in_super_kmer,
        }
    }

    /// Convert to internal MinimizerTuple
    pub fn to_internal(&self) -> MinimizerTuple {
        MinimizerTuple {
            minimizer: self.minimizer,
            pos_in_seq: self.pos_in_seq,
            pos_in_kmer: self.pos_in_kmer,
            num_kmers_in_super_kmer: self.num_kmers_in_super_kmer,
        }
    }

    /// Read from bytes (unsafe, assumes correct alignment)
    ///
    /// # Safety
    /// Caller must ensure `bytes` points to a valid `MinimizerTupleExternal`
    #[inline]
    pub unsafe fn from_bytes(bytes: *const u8) -> Self {
        // SAFETY: read_unaligned handles packed/unaligned access
        unsafe { std::ptr::read_unaligned(bytes as *const Self) }
    }

    /// Write to bytes
    pub fn to_bytes(&self) -> [u8; TUPLE_SIZE_BYTES] {
        let mut buf = [0u8; TUPLE_SIZE_BYTES];
        unsafe {
            std::ptr::copy_nonoverlapping(
                self as *const Self as *const u8,
                buf.as_mut_ptr(),
                TUPLE_SIZE_BYTES,
            );
        }
        buf
    }

    /// Read one tuple from a reader. Returns None on EOF.
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Option<Self>> {
        let mut buf = [0u8; TUPLE_SIZE_BYTES];
        match reader.read_exact(&mut buf) {
            Ok(()) => Ok(Some(unsafe { Self::from_bytes(buf.as_ptr()) })),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
            Err(e) => Err(e),
        }
    }
}

impl PartialOrd for MinimizerTupleExternal {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MinimizerTupleExternal {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Copy fields to avoid taking references to packed struct fields
        let self_min = self.minimizer;
        let other_min = other.minimizer;
        let self_pos = self.pos_in_seq;
        let other_pos = other.pos_in_seq;

        match self_min.cmp(&other_min) {
            std::cmp::Ordering::Equal => self_pos.cmp(&other_pos),
            ord => ord,
        }
    }
}

/// External sorter for minimizer tuples
///
/// Manages RAM-bounded sorting with temp file spillover and k-way merge.
pub struct ExternalSorter {
    /// Temp directory for intermediate files
    tmp_dir: PathBuf,
    /// Run identifier (timestamp-based for uniqueness)
    run_id: u64,
    /// Atomic counter for temp file IDs
    num_files: AtomicU64,
    /// RAM limit in GiB
    ram_limit_gib: usize,
    /// Number of threads
    num_threads: usize,
    /// Verbose logging
    verbose: bool,
    /// Independent ordered minimizer ranges; zero keeps the original merge.
    prefix_bits: u32,
    /// If true, merged file cleanup is owned by FileTuples (don't delete in Drop)
    merged_file_handed_off: std::sync::atomic::AtomicBool,
}

impl ExternalSorter {
    /// Create a new external sorter
    pub fn new(
        tmp_dir: impl AsRef<Path>,
        ram_limit_gib: usize,
        num_threads: usize,
        verbose: bool,
    ) -> std::io::Result<Self> {
        let tmp_dir = tmp_dir.as_ref().to_path_buf();

        // Create temp directory if it doesn't exist
        fs::create_dir_all(&tmp_dir)?;

        // Generate unique run ID from timestamp
        let run_id = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos() as u64;

        Ok(Self {
            tmp_dir,
            run_id,
            num_files: AtomicU64::new(0),
            ram_limit_gib,
            num_threads,
            verbose,
            prefix_bits: std::env::var("REINDEER_PREFIX_SORT_BITS")
                .ok()
                .and_then(|v| v.parse::<u32>().ok())
                .unwrap_or(6)
                .min(8),
            merged_file_handed_off: std::sync::atomic::AtomicBool::new(false),
        })
    }

    /// Calculate buffer size per thread in number of tuples
    ///
    /// Formula: `(ram_limit * GiB) / (2 * sizeof(tuple) * num_threads)`
    /// The factor of 2 accounts for temporary memory during parallel sort.
    pub fn buffer_size_per_thread(&self) -> usize {
        let total_bytes = self.ram_limit_gib * GIB;
        let bytes_per_thread = total_bytes / (2 * self.num_threads.max(1));
        bytes_per_thread / TUPLE_SIZE_BYTES
    }

    /// Get path for a temp file by ID
    fn temp_file_path(&self, id: u64) -> PathBuf {
        self.tmp_dir.join(format!(
            "sshash.tmp.run_{}.minimizers.{}.bin",
            self.run_id, id
        ))
    }

    /// Get path for the final merged file
    fn merged_file_path(&self) -> PathBuf {
        self.tmp_dir
            .join(format!("sshash.tmp.run_{}.minimizers.bin", self.run_id))
    }

    /// Sort a buffer and flush to a temp file
    ///
    /// Returns the file ID. Thread-safe via atomic counter.
    pub fn sort_and_flush(&self, buffer: &mut Vec<MinimizerTupleExternal>) -> std::io::Result<u64> {
        if self.prefix_bits != 0 {
            return self.distribute_sort_and_flush(buffer);
        }
        // Parallel sort
        buffer.par_sort_unstable();
        buffer.dedup_by_key(|t| (t.minimizer, t.pos_in_seq));

        // Get unique file ID
        let file_id = self.num_files.fetch_add(1, Ordering::SeqCst);
        let path = self.temp_file_path(file_id);

        if self.verbose {
            debug!("Flushing {} tuples to {:?}", buffer.len(), path);
        }

        // Write to binary file
        let file = File::create(&path)?;
        let mut writer = BufWriter::with_capacity(1024 * 1024, file);

        for tuple in buffer.iter() {
            writer.write_all(&tuple.to_bytes())?;
        }

        writer.flush()?;
        buffer.clear();

        Ok(file_id)
    }

    /// Number of temp files created
    pub fn num_files(&self) -> u64 {
        self.num_files.load(Ordering::SeqCst)
    }

    /// Merge all temp files into final sorted output
    ///
    /// Returns statistics: (num_minimizers, num_positions, num_super_kmers)
    pub fn merge(&self) -> std::io::Result<MergeResult> {
        if self.prefix_bits != 0 {
            return self.merge_prefixes();
        }
        let num_files = self.num_files();

        if num_files == 0 {
            return Ok(MergeResult::default());
        }

        if num_files == 1 {
            // Just rename the single file
            let src = self.temp_file_path(0);
            let dst = self.merged_file_path();
            fs::rename(&src, &dst)?;
            return self.scan_merged_file();
        }

        // Multiple files: k-way merge using buffered I/O
        info!("Merging {} temp files...", num_files);

        let mut merger =
            FileMergingIterator::new((0..num_files).map(|id| self.temp_file_path(id)).collect())?;

        let merged_path = self.merged_file_path();
        let file = File::create(&merged_path)?;
        let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, file);

        let mut result = MergeResult::default();
        let mut prev_minimizer = u64::MAX;
        let mut prev_pos_in_seq = u64::MAX;

        while merger.has_next() {
            let tuple = merger.current();

            // Track statistics; the scheme is forward, so no (minimizer,
            // pos_in_seq) pair may occur twice among the sorted tuples —
            // everything downstream sizes the index by the tuple count, so a
            // violation must stop the build rather than corrupt it.
            if tuple.minimizer != prev_minimizer {
                prev_minimizer = tuple.minimizer;
                result.num_minimizers += 1;
                result.num_positions += 1;
            } else if tuple.pos_in_seq != prev_pos_in_seq {
                result.num_positions += 1;
            } else {
                merger.next()?;
                continue;
            }
            prev_pos_in_seq = tuple.pos_in_seq;
            result.num_super_kmers += 1;

            writer.write_all(&tuple.to_bytes())?;
            merger.next()?;

            if self.verbose && result.num_super_kmers % 100_000_000 == 0 {
                info!("Merged {} tuples...", result.num_super_kmers);
            }
        }

        writer.flush()?;
        drop(merger);

        // Remove temp files
        for id in 0..num_files {
            let _ = fs::remove_file(self.temp_file_path(id));
        }

        info!(
            "Merge complete: {} minimizers, {} positions, {} super-kmers",
            result.num_minimizers, result.num_positions, result.num_super_kmers
        );

        Ok(result)
    }

    /// Scan merged file to compute statistics (for single-file case)
    fn scan_merged_file(&self) -> std::io::Result<MergeResult> {
        let path = self.merged_file_path();
        let file = File::open(&path)?;
        let mut reader = BufReader::with_capacity(READER_BUF_SIZE, file);

        let mut result = MergeResult::default();
        let mut prev_minimizer = u64::MAX;
        let mut prev_pos_in_seq = u64::MAX;

        while let Some(tuple) = MinimizerTupleExternal::read_from(&mut reader)? {
            if tuple.minimizer != prev_minimizer {
                prev_minimizer = tuple.minimizer;
                result.num_minimizers += 1;
                result.num_positions += 1;
            } else if tuple.pos_in_seq != prev_pos_in_seq {
                result.num_positions += 1;
            } else {
                // See merge(): the forward scheme forbids duplicate pairs.
                let (bad_min, bad_pos) = (tuple.minimizer, tuple.pos_in_seq);
                return Err(std::io::Error::other(format!(
                    "the minimizer scheme is not forward: (minimizer {bad_min:#x}, \
                     pos_in_seq {bad_pos}) occurs in more than one super-kmer"
                )));
            }
            prev_pos_in_seq = tuple.pos_in_seq;
            result.num_super_kmers += 1;
        }

        Ok(result)
    }

    /// Read merged tuples into memory as internal MinimizerTuples
    ///
    /// Call this after `merge()` to get the final sorted tuples.
    /// For large datasets, prefer [`Self::open_merged_file`] to avoid full materialization.
    pub fn read_merged_tuples(&self) -> std::io::Result<Vec<MinimizerTuple>> {
        let path = self.merged_file_path();
        let file = File::open(&path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        let num_tuples = mmap.len() / TUPLE_SIZE_BYTES;
        let mut tuples = Vec::with_capacity(num_tuples);

        for i in 0..num_tuples {
            let offset = i * TUPLE_SIZE_BYTES;
            let ext = unsafe { MinimizerTupleExternal::from_bytes(mmap.as_ptr().add(offset)) };
            tuples.push(ext.to_internal());
        }

        Ok(tuples)
    }

    /// Open the merged file for sequential buffered access.
    ///
    /// Returns a [`FileTuples`] handle that provides sequential access via
    /// buffered readers. Each scan pass opens a fresh `BufReader`, keeping
    /// RSS proportional to the buffer size (~4 MB) instead of the file size.
    ///
    /// **Important**: The returned `FileTuples` takes ownership of cleanup
    /// responsibility for the merged file.
    pub fn open_merged_file(&self) -> std::io::Result<FileTuples> {
        let path = self.merged_file_path();
        let file_len = fs::metadata(&path)?.len() as usize;
        let num_tuples = file_len / TUPLE_SIZE_BYTES;

        // Mark that cleanup responsibility has been transferred to FileTuples
        self.merged_file_handed_off.store(true, Ordering::SeqCst);

        info!(
            "Opened merged file: {} tuples ({:.2} GB on disk)",
            num_tuples,
            file_len as f64 / GIB as f64
        );

        Ok(FileTuples {
            path: path.clone(),
            num_tuples,
        })
    }

    /// Remove the merged file (cleanup)
    pub fn remove_merged_file(&self) -> std::io::Result<()> {
        let path = self.merged_file_path();
        if path.exists() {
            fs::remove_file(path)?;
        }
        Ok(())
    }

    fn prefix_path(&self, prefix: usize, run: u64) -> PathBuf {
        self.tmp_dir.join(format!(
            "sshash.tmp.run_{}.prefix_{prefix}.{run}.bin",
            self.run_id
        ))
    }

    fn prefix_merged_path(&self, prefix: usize) -> PathBuf {
        self.tmp_dir.join(format!(
            "sshash.tmp.run_{}.prefix_{prefix}.merged.bin",
            self.run_id
        ))
    }

    fn distribute_sort_and_flush(
        &self,
        buffer: &mut Vec<MinimizerTupleExternal>,
    ) -> std::io::Result<u64> {
        use std::mem::MaybeUninit;
        let prefixes = 1usize << self.prefix_bits;
        let prefix_of =
            |tuple: &MinimizerTupleExternal| (tuple.minimizer >> (64 - self.prefix_bits)) as usize;
        let mut counts = vec![0; prefixes];
        for tuple in buffer.iter() {
            counts[prefix_of(tuple)] += 1;
        }
        let mut offsets = Vec::with_capacity(prefixes + 1);
        offsets.push(0usize);
        for &count in &counts {
            offsets.push(offsets.last().unwrap() + count);
        }
        let mut cursors = offsets[..prefixes].to_vec();
        let mut grouped: Vec<MaybeUninit<MinimizerTupleExternal>> =
            Vec::with_capacity(buffer.len());
        // SAFETY: uninitialized MaybeUninit values are valid; counting creates
        // disjoint, exact spans, each completely filled by the scatter below.
        unsafe {
            grouped.set_len(buffer.len());
        }
        for &tuple in buffer.iter() {
            let p = prefix_of(&tuple);
            grouped[cursors[p]].write(tuple);
            cursors[p] += 1;
        }
        assert!(
            cursors
                .iter()
                .zip(&offsets[1..])
                .all(|(actual, end)| actual == end)
        );
        let (length, capacity) = (grouped.len(), grouped.capacity());
        let ptr = grouped.as_mut_ptr().cast::<MinimizerTupleExternal>();
        std::mem::forget(grouped);
        // SAFETY: all slots initialized; layout/alignment are identical.
        let mut grouped = unsafe { Vec::from_raw_parts(ptr, length, capacity) };
        let run = self.num_files.fetch_add(1, Ordering::SeqCst);
        let mut remaining = grouped.as_mut_slice();
        let jobs: Vec<_> = counts
            .iter()
            .enumerate()
            .map(|(prefix, &count)| {
                let (part, rest) = std::mem::take(&mut remaining).split_at_mut(count);
                remaining = rest;
                (prefix, part)
            })
            .collect();
        jobs.into_par_iter()
            .try_for_each(|(prefix, part)| -> std::io::Result<()> {
                if part.is_empty() {
                    return Ok(());
                }
                part.par_sort_unstable();
                let file = File::create(self.prefix_path(prefix, run))?;
                let mut writer = BufWriter::with_capacity(256 * 1024, file);
                let mut previous = None;
                for tuple in part {
                    let pair = (tuple.minimizer, tuple.pos_in_seq);
                    if previous != Some(pair) {
                        writer.write_all(&tuple.to_bytes())?;
                        previous = Some(pair);
                    }
                }
                writer.flush()
            })?;
        buffer.clear();
        Ok(run)
    }

    fn merge_prefix(&self, prefix: usize, reader_bytes: usize) -> std::io::Result<MergeResult> {
        let paths: Vec<_> = (0..self.num_files())
            .map(|run| self.prefix_path(prefix, run))
            .filter(|path| path.exists())
            .collect();
        if paths.is_empty() {
            return Ok(MergeResult::default());
        }
        let mut merger = FileMergingIterator::with_buffer(paths.clone(), reader_bytes)?;
        let mut writer =
            BufWriter::with_capacity(1024 * 1024, File::create(self.prefix_merged_path(prefix))?);
        let mut result = MergeResult::default();
        let mut previous = None;
        while merger.has_next() {
            let tuple = merger.current();
            let pair = (tuple.minimizer, tuple.pos_in_seq);
            if previous != Some(pair) {
                if previous.is_none_or(|(minimizer, _)| minimizer != tuple.minimizer) {
                    result.num_minimizers += 1;
                }
                result.num_positions += 1;
                result.num_super_kmers += 1;
                writer.write_all(&tuple.to_bytes())?;
                previous = Some(pair);
            }
            merger.next()?;
        }
        writer.flush()?;
        drop(writer);
        drop(merger);
        // Reclaim input runs before concatenation, keeping temporary space bounded.
        for path in paths {
            fs::remove_file(path)?;
        }
        Ok(result)
    }

    fn merge_prefixes(&self) -> std::io::Result<MergeResult> {
        let prefixes = 1usize << self.prefix_bits;
        let runs = self.num_files() as usize;
        // Leave descriptors for GGCAT, logs and unrelated handles in this process.
        let descriptor_limit = fs::read_to_string("/proc/self/limits")
            .ok()
            .and_then(|text| {
                text.lines()
                    .find(|line| line.starts_with("Max open files"))
                    .and_then(|line| line.split_whitespace().nth(3))
                    .and_then(|n| n.parse::<usize>().ok())
            })
            .unwrap_or(1024);
        let concurrent = rayon::current_num_threads().min(prefixes).min(
            descriptor_limit
                .saturating_sub(256)
                .checked_div(runs + 1)
                .unwrap_or(0)
                .max(1),
        );
        let reader_bytes = 128 * 1024;
        let started = std::time::Instant::now();
        let mut total = MergeResult::default();
        for first in (0..prefixes).step_by(concurrent) {
            let results: Vec<_> = (first..(first + concurrent).min(prefixes))
                .into_par_iter()
                .map(|prefix| self.merge_prefix(prefix, reader_bytes))
                .collect();
            for result in results {
                let result = result?;
                total.num_minimizers += result.num_minimizers;
                total.num_positions += result.num_positions;
                total.num_super_kmers += result.num_super_kmers;
            }
        }
        let mut output = File::create(self.merged_file_path())?;
        for prefix in 0..prefixes {
            let path = self.prefix_merged_path(prefix);
            if path.exists() {
                std::io::copy(&mut File::open(&path)?, &mut output)?;
                fs::remove_file(path)?;
            }
        }
        output.flush()?;
        eprintln!(
            "SSHash prefix merge: {:.3}s; prefixes={prefixes}; concurrent={concurrent}; runs={runs}",
            started.elapsed().as_secs_f64()
        );
        Ok(total)
    }
}

impl Drop for ExternalSorter {
    fn drop(&mut self) {
        // Clean up any remaining temp files
        for id in 0..self.num_files() {
            let _ = fs::remove_file(self.temp_file_path(id));
            if self.prefix_bits != 0 {
                for prefix in 0..1usize << self.prefix_bits {
                    let _ = fs::remove_file(self.prefix_path(prefix, id));
                }
            }
        }
        if self.prefix_bits != 0 {
            for prefix in 0..1usize << self.prefix_bits {
                let _ = fs::remove_file(self.prefix_merged_path(prefix));
            }
        }
        // Only clean up merged file if it wasn't handed off to FileTuples
        if !self.merged_file_handed_off.load(Ordering::SeqCst) {
            let _ = fs::remove_file(self.merged_file_path());
        }
    }
}

// ---------------------------------------------------------------------------
// Sequential buffered tuple access (no mmap)
// ---------------------------------------------------------------------------

/// Handle to the merged external sort file for sequential buffered access.
///
/// Unlike the previous `MmapTuples`, this does NOT mmap the file. Instead,
/// each access pass opens a fresh `BufReader`, reading tuples sequentially
/// with a ~4 MB buffer. This keeps RSS proportional to the buffer size
/// instead of the file size (~8-10 GB for large datasets).
pub struct FileTuples {
    /// Path to the merged file
    path: PathBuf,
    /// Total number of tuples
    num_tuples: usize,
}

impl FileTuples {
    /// Total number of tuples in the merged file.
    #[inline]
    pub fn num_tuples(&self) -> usize {
        self.num_tuples
    }

    /// Path to the underlying merged file (for deriving temp file names).
    #[inline]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Create an iterator that yields one [`BucketScan`] per unique minimizer,
    /// scanning the file sequentially via buffered I/O.
    ///
    /// Each call opens a fresh `BufReader` — no state is retained between
    /// passes, and no mmap pages linger in memory.
    pub fn bucket_iter(&self) -> std::io::Result<FileBucketIter> {
        let file = File::open(&self.path)?;
        let reader = BufReader::with_capacity(READER_BUF_SIZE, file);
        Ok(FileBucketIter {
            reader,
            pos: 0,
            num_tuples: self.num_tuples,
            lookahead: None,
        })
    }

    /// Create a sequential tuple reader for forward-only access.
    ///
    /// Returns tuples one at a time via `read_next()`. Used by the
    /// sparse/skew index builder for its sequential pass over all tuples.
    pub fn sequential_reader(&self) -> std::io::Result<SequentialTupleReader> {
        let file = File::open(&self.path)?;
        let reader = BufReader::with_capacity(READER_BUF_SIZE, file);
        Ok(SequentialTupleReader { reader })
    }
}

impl Drop for FileTuples {
    fn drop(&mut self) {
        // Clean up the merged file when done
        let _ = fs::remove_file(&self.path);
    }
}

/// Sequential reader for tuples from the merged file.
///
/// Reads tuples one at a time via `read_next()`. All access is
/// strictly forward — no seeking or random access.
pub struct SequentialTupleReader {
    reader: BufReader<File>,
}

impl SequentialTupleReader {
    /// Read the next tuple, or None on EOF.
    #[inline]
    pub fn read_next(&mut self) -> std::io::Result<Option<MinimizerTupleExternal>> {
        MinimizerTupleExternal::read_from(&mut self.reader)
    }
}

/// Sequential iterator over buckets in the merged file via buffered I/O.
///
/// Each call to `next()` returns a [`BucketScan`] describing one bucket
/// (a contiguous group of tuples sharing the same minimizer value).
/// Owns a `BufReader` — no mmap pages are kept resident.
pub struct FileBucketIter {
    reader: BufReader<File>,
    /// Current tuple position index
    pos: usize,
    /// Total tuples in file
    num_tuples: usize,
    /// Lookahead tuple (the first tuple of the NEXT bucket, read while
    /// scanning the current bucket's boundary)
    lookahead: Option<MinimizerTupleExternal>,
}

impl Iterator for FileBucketIter {
    type Item = BucketScan;

    fn next(&mut self) -> Option<BucketScan> {
        if self.pos >= self.num_tuples {
            return None;
        }

        // Get the first tuple of this bucket — either from lookahead or by reading
        let first = if let Some(la) = self.lookahead.take() {
            la
        } else {
            match MinimizerTupleExternal::read_from(&mut self.reader) {
                Ok(Some(t)) => t,
                _ => return None,
            }
        };

        let start = self.pos;
        let minimizer = first.minimizer;
        let mut cached_size = 1usize;
        let mut prev_pos_in_seq = first.pos_in_seq;
        let mut num_kmers = first.num_kmers_in_super_kmer as u64;
        self.pos += 1;

        // Read tuples until we hit a different minimizer. The scheme is
        // forward, so every tuple is one super-kmer; a duplicate pos_in_seq
        // must abort — cached_size sizes the EF and offsets arrays, so a
        // violation here would silently corrupt the index.
        while self.pos < self.num_tuples {
            match MinimizerTupleExternal::read_from(&mut self.reader) {
                Ok(Some(t)) => {
                    if t.minimizer != minimizer {
                        // This tuple belongs to the next bucket — save as lookahead
                        self.lookahead = Some(t);
                        break;
                    }
                    let t_pos = t.pos_in_seq;
                    assert_ne!(
                        t_pos, prev_pos_in_seq,
                        "the minimizer scheme is not forward: (minimizer {minimizer:#x}, \
                         pos_in_seq {prev_pos_in_seq}) occurs in more than one super-kmer"
                    );
                    cached_size += 1;
                    prev_pos_in_seq = t_pos;
                    num_kmers += t.num_kmers_in_super_kmer as u64;
                    self.pos += 1;
                }
                _ => break,
            }
        }

        Some(BucketScan {
            minimizer,
            cached_size,
            start_tuple_idx: start as u64,
            num_tuples: self.pos - start,
            num_kmers,
        })
    }
}

/// Metadata about a single bucket, gathered during a sequential scan.
#[derive(Debug, Clone, Copy)]
pub struct BucketScan {
    /// The minimizer value for this bucket.
    pub minimizer: u64,
    /// Number of unique super-kmers (distinct `pos_in_seq` values).
    pub cached_size: usize,
    /// Index of the first tuple in the file.
    pub start_tuple_idx: u64,
    /// Number of tuples in this bucket.
    pub num_tuples: usize,
    /// Total k-mers across all tuples in this bucket.
    pub num_kmers: u64,
}

/// Result of merge operation
#[derive(Debug, Default, Clone, Copy)]
pub struct MergeResult {
    /// Number of distinct minimizers
    pub num_minimizers: u64,
    /// Total number of positions
    pub num_positions: u64,
    /// Number of super-k-mers
    pub num_super_kmers: u64,
}

/// K-way merge iterator using buffered file I/O.
///
/// Uses buffered readers instead of mmap to avoid keeping all temp file pages
/// resident in memory. Each file gets a ~1 MB read buffer.
///
/// Memory footprint: ~1 MB per file (buffer) + 18 bytes per file (current tuple).
struct FileMergingIterator {
    /// Buffered readers for each input file
    readers: Vec<BufReader<File>>,
    /// Current tuple from each file (None if exhausted)
    current_tuples: Vec<Option<MinimizerTupleExternal>>,
    /// Minimum (minimizer, sequence position, file ID). File ID preserves
    /// deterministic ordering of equal keys across spill files.
    heap: BinaryHeap<Reverse<(u64, u64, usize)>>,
}

impl FileMergingIterator {
    fn new(paths: Vec<PathBuf>) -> std::io::Result<Self> {
        Self::with_buffer(paths, 1024 * 1024)
    }

    fn with_buffer(paths: Vec<PathBuf>, buffer_bytes: usize) -> std::io::Result<Self> {
        let num_files = paths.len();
        let mut readers = Vec::with_capacity(num_files);
        let mut current_tuples = Vec::with_capacity(num_files);
        let mut heads = Vec::with_capacity(num_files);

        for (id, path) in paths.iter().enumerate() {
            let file = File::open(path)?;
            let mut reader = BufReader::with_capacity(buffer_bytes, file);
            let tuple = MinimizerTupleExternal::read_from(&mut reader)?;
            if let Some(tuple) = tuple {
                heads.push(Reverse((tuple.minimizer, tuple.pos_in_seq, id)));
            }
            current_tuples.push(tuple);
            readers.push(reader);
        }

        Ok(Self {
            readers,
            current_tuples,
            heap: BinaryHeap::from(heads),
        })
    }

    fn has_next(&self) -> bool {
        !self.heap.is_empty()
    }

    fn current(&self) -> MinimizerTupleExternal {
        let id = self.heap.peek().unwrap().0.2;
        self.current_tuples[id].unwrap()
    }

    fn next(&mut self) -> std::io::Result<()> {
        let Some(&Reverse((_, _, id))) = self.heap.peek() else {
            return Ok(());
        };
        let tuple = MinimizerTupleExternal::read_from(&mut self.readers[id])?;
        self.current_tuples[id] = tuple;
        if let Some(tuple) = tuple {
            // Replace the root in place: O(log files), without allocation.
            *self.heap.peek_mut().unwrap() = Reverse((tuple.minimizer, tuple.pos_in_seq, id));
        } else {
            self.heap.pop();
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    #[ignore = "isolated before/after external sorting benchmark"]
    fn paper_external_sort_trial() {
        let n = std::env::var("REINDEER_TRIAL_KEYS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(10_000_000usize);
        let threads = std::env::var("REINDEER_TRIAL_THREADS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(64usize);
        let directory = PathBuf::from(
            std::env::var("REINDEER_TRIAL_DIRECTORY").expect("isolated trial directory"),
        );
        assert!(
            !directory.exists(),
            "Refuse to overwrite a prior sorting trial"
        );
        let sorter = ExternalSorter::new(&directory, 2, threads, false).unwrap();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let run_size = n.div_ceil(threads).max(1);
        let tuple = |position: usize| MinimizerTupleExternal {
            minimizer: (position as u64 % 1_000_003).wrapping_mul(0x9e3779b97f4a7c15),
            pos_in_seq: position as u64,
            pos_in_kmer: (position % 19) as u8,
            num_kmers_in_super_kmer: (position % 7 + 1) as u8,
        };
        let started = std::time::Instant::now();
        pool.install(|| {
            (0..n)
                .step_by(run_size)
                .collect::<Vec<_>>()
                .into_par_iter()
                .for_each(|begin| {
                    let mut buffer = (begin..(begin + run_size).min(n)).map(tuple).collect();
                    sorter.sort_and_flush(&mut buffer).unwrap();
                })
        });
        let stats = pool.install(|| sorter.merge().unwrap());
        let seconds = started.elapsed().as_secs_f64();
        assert_eq!(stats.num_positions, n as u64);
        let file = sorter.open_merged_file().unwrap();
        let mut reader = file.sequential_reader().unwrap();
        let mut used = vec![false; n];
        let mut previous = None;
        while let Some(actual) = reader.read_next().unwrap() {
            let position = actual.pos_in_seq as usize;
            assert!(position < n && !used[position]);
            used[position] = true;
            assert_eq!(actual, tuple(position));
            let pair = (actual.minimizer, actual.pos_in_seq);
            assert!(previous.is_none_or(|prev| prev < pair));
            previous = Some(pair);
        }
        assert!(used.iter().all(|&occupied| occupied));
        eprintln!(
            "PAPER_TRIAL {{\"kind\":\"sort\",\"keys\":{n},\"threads\":{threads},\"construction_seconds\":{seconds},\"prefix_bits\":{},\"verified\":true}}",
            sorter.prefix_bits
        );
    }

    #[test]
    fn prefix_sort_matches_global_sort_with_duplicates_and_empty_ranges() {
        let mut input: Vec<_> = (0..12000u64)
            .map(|i| MinimizerTupleExternal {
                minimizer: (i % 1703).wrapping_mul(0x9e3779b97f4a7c15),
                pos_in_seq: i % 103,
                pos_in_kmer: (i % 17) as u8,
                num_kmers_in_super_kmer: 1,
            })
            .collect();
        input.extend_from_within(0..512);
        input.extend_from_within(3000..3512);
        let mut expected = input.clone();
        expected.sort_unstable();
        expected.dedup_by_key(|tuple| (tuple.minimizer, tuple.pos_in_seq));
        for bits in [0, 1, 6, 8] {
            let directory = TempDir::new().unwrap();
            let mut sorter = ExternalSorter::new(directory.path(), 1, 4, false).unwrap();
            sorter.prefix_bits = bits;
            // Empty spill plus many overlapping runs exercise per-prefix heap merging.
            sorter.sort_and_flush(&mut Vec::new()).unwrap();
            for chunk in input.chunks(311) {
                sorter.sort_and_flush(&mut chunk.to_vec()).unwrap();
            }
            let stats = sorter.merge().unwrap();
            assert_eq!(stats.num_positions, expected.len() as u64);
            assert_eq!(stats.num_super_kmers, expected.len() as u64);
            let actual = sorter.read_merged_tuples().unwrap();
            let actual_pairs: Vec<_> = actual
                .iter()
                .map(|t| {
                    (
                        t.minimizer,
                        t.pos_in_seq,
                        t.pos_in_kmer,
                        t.num_kmers_in_super_kmer,
                    )
                })
                .collect();
            assert_eq!(
                actual_pairs,
                expected
                    .iter()
                    .map(|t| (
                        t.minimizer,
                        t.pos_in_seq,
                        t.pos_in_kmer,
                        t.num_kmers_in_super_kmer
                    ))
                    .collect::<Vec<_>>()
            );
            let distinct = expected
                .iter()
                .map(|t| t.minimizer)
                .collect::<std::collections::HashSet<_>>()
                .len();
            assert_eq!(stats.num_minimizers, distinct as u64);
            drop(sorter);
            assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 0);
        }
    }

    #[test]
    fn prefix_sort_handles_maximum_keys_and_all_empty_spills() {
        for input in [
            Vec::new(),
            vec![MinimizerTupleExternal {
                minimizer: u64::MAX,
                pos_in_seq: u64::MAX,
                pos_in_kmer: 0,
                num_kmers_in_super_kmer: 1,
            }],
        ] {
            let directory = TempDir::new().unwrap();
            let mut sorter = ExternalSorter::new(directory.path(), 1, 2, false).unwrap();
            sorter.prefix_bits = 6;
            sorter.sort_and_flush(&mut input.clone()).unwrap();
            let stats = sorter.merge().unwrap();
            assert_eq!(stats.num_positions, input.len() as u64);
            assert_eq!(stats.num_minimizers, input.len() as u64);
            assert_eq!(sorter.read_merged_tuples().unwrap().len(), input.len());
        }
    }

    #[test]
    fn heap_merge_preserves_order_and_ties_across_many_runs() {
        let directory = TempDir::new().unwrap();
        let mut paths = Vec::new();
        let mut expected = Vec::new();
        for file_id in 0..129 {
            let path = directory.path().join(format!("run-{file_id}"));
            let mut output = File::create(&path).unwrap();
            let mut tuples = Vec::new();
            // Include empty runs, shared keys, and unequal run lengths.
            for i in 0..file_id % 17 {
                tuples.push(MinimizerTupleExternal {
                    minimizer: ((i * 19 + file_id) % 11) as u64,
                    pos_in_seq: (i % 5) as u64,
                    pos_in_kmer: file_id as u8,
                    num_kmers_in_super_kmer: 1,
                });
            }
            tuples.sort();
            for tuple in tuples {
                output.write_all(&tuple.to_bytes()).unwrap();
                expected.push((tuple, file_id));
            }
            paths.push(path);
        }
        // Independent global stable sort is the reference, not another merge.
        expected.sort_by_key(|(t, id)| (t.minimizer, t.pos_in_seq, *id));
        let mut merge = FileMergingIterator::new(paths).unwrap();
        let mut actual = Vec::new();
        while merge.has_next() {
            actual.push(merge.current());
            merge.next().unwrap();
        }
        assert_eq!(
            actual,
            expected.into_iter().map(|(t, _)| t).collect::<Vec<_>>()
        );
        merge.next().unwrap();
        assert!(!merge.has_next());
        assert!(!FileMergingIterator::new(Vec::new()).unwrap().has_next());
    }

    #[test]
    fn test_tuple_external_size() {
        // Verify packed size
        assert_eq!(
            std::mem::size_of::<MinimizerTupleExternal>(),
            TUPLE_SIZE_BYTES
        );
    }

    #[test]
    fn test_tuple_roundtrip() {
        let tuple = MinimizerTupleExternal {
            minimizer: 12345,
            pos_in_seq: 67890,
            pos_in_kmer: 5,
            num_kmers_in_super_kmer: 3,
        };

        let bytes = tuple.to_bytes();
        let recovered = unsafe { MinimizerTupleExternal::from_bytes(bytes.as_ptr()) };

        assert_eq!(tuple, recovered);
    }

    #[test]
    fn test_external_sorter_basic() {
        let tmp_dir = TempDir::new().unwrap();
        let sorter = ExternalSorter::new(tmp_dir.path(), 1, 2, false).unwrap();

        // Buffer size should be reasonable
        let buf_size = sorter.buffer_size_per_thread();
        assert!(buf_size > 0);

        // Test sort_and_flush
        let mut buffer: Vec<MinimizerTupleExternal> = vec![
            MinimizerTupleExternal {
                minimizer: 100,
                pos_in_seq: 10,
                pos_in_kmer: 1,
                num_kmers_in_super_kmer: 2,
            },
            MinimizerTupleExternal {
                minimizer: 50,
                pos_in_seq: 20,
                pos_in_kmer: 3,
                num_kmers_in_super_kmer: 1,
            },
        ];

        sorter.sort_and_flush(&mut buffer).unwrap();
        assert!(buffer.is_empty());
        assert_eq!(sorter.num_files(), 1);
    }

    #[test]
    fn test_external_sorter_merge() {
        let tmp_dir = TempDir::new().unwrap();
        let sorter = ExternalSorter::new(tmp_dir.path(), 1, 2, false).unwrap();

        // Create multiple temp files
        let mut buffer1: Vec<MinimizerTupleExternal> = vec![
            MinimizerTupleExternal {
                minimizer: 10,
                pos_in_seq: 1,
                pos_in_kmer: 0,
                num_kmers_in_super_kmer: 1,
            },
            MinimizerTupleExternal {
                minimizer: 30,
                pos_in_seq: 3,
                pos_in_kmer: 0,
                num_kmers_in_super_kmer: 1,
            },
        ];
        sorter.sort_and_flush(&mut buffer1).unwrap();

        let mut buffer2: Vec<MinimizerTupleExternal> = vec![
            MinimizerTupleExternal {
                minimizer: 20,
                pos_in_seq: 2,
                pos_in_kmer: 0,
                num_kmers_in_super_kmer: 1,
            },
            MinimizerTupleExternal {
                minimizer: 40,
                pos_in_seq: 4,
                pos_in_kmer: 0,
                num_kmers_in_super_kmer: 1,
            },
        ];
        sorter.sort_and_flush(&mut buffer2).unwrap();

        assert_eq!(sorter.num_files(), 2);

        // Merge
        let result = sorter.merge().unwrap();
        assert_eq!(result.num_super_kmers, 4);
        assert_eq!(result.num_minimizers, 4);

        // Read merged tuples
        let tuples = sorter.read_merged_tuples().unwrap();
        assert_eq!(tuples.len(), 4);

        // Verify sorted order
        assert_eq!(tuples[0].minimizer, 10);
        assert_eq!(tuples[1].minimizer, 20);
        assert_eq!(tuples[2].minimizer, 30);
        assert_eq!(tuples[3].minimizer, 40);
    }

    #[test]
    fn test_tuple_ordering() {
        let t1 = MinimizerTupleExternal {
            minimizer: 100,
            pos_in_seq: 50,
            pos_in_kmer: 0,
            num_kmers_in_super_kmer: 1,
        };
        let t2 = MinimizerTupleExternal {
            minimizer: 100,
            pos_in_seq: 60,
            pos_in_kmer: 0,
            num_kmers_in_super_kmer: 1,
        };
        let t3 = MinimizerTupleExternal {
            minimizer: 200,
            pos_in_seq: 10,
            pos_in_kmer: 0,
            num_kmers_in_super_kmer: 1,
        };

        assert!(t1 < t2); // Same minimizer, different pos
        assert!(t1 < t3); // Different minimizer
        assert!(t2 < t3); // Different minimizer
    }

    #[test]
    fn test_file_tuples_bucket_iter() {
        let tmp_dir = TempDir::new().unwrap();
        let sorter = ExternalSorter::new(tmp_dir.path(), 1, 2, false).unwrap();

        // Create a sorted file with 2 buckets: minimizer 10 (2 tuples), minimizer 20 (1 tuple)
        let mut buffer: Vec<MinimizerTupleExternal> = vec![
            MinimizerTupleExternal {
                minimizer: 10,
                pos_in_seq: 1,
                pos_in_kmer: 0,
                num_kmers_in_super_kmer: 3,
            },
            MinimizerTupleExternal {
                minimizer: 10,
                pos_in_seq: 2,
                pos_in_kmer: 0,
                num_kmers_in_super_kmer: 2,
            },
            MinimizerTupleExternal {
                minimizer: 20,
                pos_in_seq: 5,
                pos_in_kmer: 1,
                num_kmers_in_super_kmer: 4,
            },
        ];
        sorter.sort_and_flush(&mut buffer).unwrap();
        sorter.merge().unwrap();
        let ft = sorter.open_merged_file().unwrap();

        let buckets: Vec<BucketScan> = ft.bucket_iter().unwrap().collect();
        assert_eq!(buckets.len(), 2);
        assert_eq!(buckets[0].minimizer, 10);
        assert_eq!(buckets[0].cached_size, 2);
        assert_eq!(buckets[0].num_tuples, 2);
        assert_eq!(buckets[0].num_kmers, 5); // 3 + 2
        assert_eq!(buckets[1].minimizer, 20);
        assert_eq!(buckets[1].cached_size, 1);
        assert_eq!(buckets[1].num_tuples, 1);
        assert_eq!(buckets[1].num_kmers, 4);
    }

    #[test]
    fn test_sequential_reader() {
        let tmp_dir = TempDir::new().unwrap();
        let sorter = ExternalSorter::new(tmp_dir.path(), 1, 2, false).unwrap();

        let mut buffer: Vec<MinimizerTupleExternal> = vec![
            MinimizerTupleExternal {
                minimizer: 5,
                pos_in_seq: 10,
                pos_in_kmer: 0,
                num_kmers_in_super_kmer: 1,
            },
            MinimizerTupleExternal {
                minimizer: 15,
                pos_in_seq: 20,
                pos_in_kmer: 0,
                num_kmers_in_super_kmer: 2,
            },
        ];
        sorter.sort_and_flush(&mut buffer).unwrap();
        sorter.merge().unwrap();
        let ft = sorter.open_merged_file().unwrap();

        let mut reader = ft.sequential_reader().unwrap();
        let t1 = reader.read_next().unwrap().unwrap();
        let t1_min = t1.minimizer; // copy to avoid packed struct ref
        assert_eq!(t1_min, 5);
        let t2 = reader.read_next().unwrap().unwrap();
        let t2_min = t2.minimizer;
        assert_eq!(t2_min, 15);
        assert!(reader.read_next().unwrap().is_none());
    }
}
