//! One query pass per genome: SIMD lookup and exact alignment evidence are
//! computed together. Reduction never reads genome bases or compressed caches.
use super::*;
mod compact;
mod postings;
mod restriction;
mod shared;
pub use shared::{StateInfo, locate_shared_mums, read_state_info};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, BufWriter, Read, Write},
    ops::Range,
    sync::atomic::{AtomicUsize, Ordering},
    thread,
    time::Instant,
};

#[derive(serde::Serialize, serde::Deserialize)]
struct ContigMeta {
    name: String,
    start: u64,
    len: usize,
}

#[derive(serde::Serialize, serde::Deserialize)]
struct GenomeMeta {
    input: PathBuf,
    contigs: Vec<ContigMeta>,
}

struct PackedContig {
    forward: PackedSeqVec,
    reverse: PackedSeqVec,
    runs: Vec<Range<usize>>,
}

// A reduction without an early-exit branch lets LLVM vectorize the common
// A/C/G/T-only case. Ambiguous contigs retain the exact run-boundary scan.
fn all_dna(sequence: &[u8]) -> bool {
    sequence
        .iter()
        .fold(0u8, |bad, &base| bad | u8::from(!is_dna(base)))
        == 0
}

impl PackedContig {
    fn new(sequence: &[u8], k: usize) -> Self {
        let forward = PackedSeqVec::from_ascii(sequence);
        let reverse = forward.as_slice().to_revcomp();
        let mut runs = Vec::new();
        if all_dna(sequence) {
            if sequence.len() >= k {
                runs.push(0..sequence.len());
            }
            return Self {
                forward,
                reverse,
                runs,
            };
        }
        let mut pos = 0;
        while pos < sequence.len() {
            while pos < sequence.len() && !is_dna(sequence[pos]) {
                pos += 1;
            }
            let start = pos;
            while pos < sequence.len() && is_dna(sequence[pos]) {
                pos += 1;
            }
            if pos - start >= k {
                runs.push(start..pos);
            }
        }
        Self {
            forward,
            reverse,
            runs,
        }
    }

    fn len(&self) -> usize {
        self.forward.len()
    }

    fn bounds(&self, start: usize, len: usize) -> io::Result<Range<usize>> {
        let next = self.runs.partition_point(|run| run.start <= start);
        next.checked_sub(1)
            .and_then(|i| self.runs.get(i))
            .filter(|r| start + len <= r.end)
            .cloned()
            .ok_or_else(|| invalid("hit outside an unambiguous DNA run"))
    }
}

// A proto-MUM has many independently usable positions. Each observed exact
// range retains its proto coordinates and its physical genome occurrence.
#[derive(Clone, Copy, Debug)]
struct Fragment {
    begin: usize,
    end: usize,
    hit: Hit,
}

struct Reference {
    contigs: Vec<PackedContig>,
    meta: GenomeMeta,
    fragments: Vec<Vec<Fragment>>,
}

struct Anchor {
    hit: Hit,
    len: usize,
}

#[derive(Default, Clone, Copy)]
struct ScanStats {
    genomes: usize,
    contexts: usize,
    alignments: usize,
    cache_hits: usize,
    hot_cache_hits: usize,
    context_witness_hits: usize,
    reference_seed_hits: usize,
    absent_windows_skipped: usize,
    unique_alignment_windows_skipped: usize,
    extended_bases: u64,
}

impl ScanStats {
    fn add(&mut self, other: Self) {
        self.genomes += other.genomes;
        self.contexts += other.contexts;
        self.alignments += other.alignments;
        self.cache_hits += other.cache_hits;
        self.hot_cache_hits += other.hot_cache_hits;
        self.context_witness_hits += other.context_witness_hits;
        self.reference_seed_hits += other.reference_seed_hits;
        self.absent_windows_skipped += other.absent_windows_skipped;
        self.unique_alignment_windows_skipped += other.unique_alignment_windows_skipped;
        self.extended_bases += other.extended_bases;
    }
}

struct SeedPosting {
    pattern: u32,
    offset: u32,
    next: u32,
    reverse: bool,
}

struct FragmentIndex {
    options: LocateOptions,
    patterns: Vec<PackedContig>,
    // Compatibility index only outside SSHash's supported k=3..64 range.
    // Supported k always uses compact or posting_dictionary, including low k.
    heads: FxHashMap<u32, u32>,
    postings: Vec<SeedPosting>,
    posting_hashes: Vec<u32>,
    position_offsets: Vec<usize>,
    active: Vec<u64>,
    unique_runs: Vec<u32>,
    probe: Vec<u64>,
    compact: Option<compact::CompactProtoIndex>,
    posting_dictionary: Option<postings::PostingDictionary>,
}

impl FragmentIndex {
    fn new(source: ProtoIndex) -> io::Result<Self> {
        let mut this = Self {
            options: source.options.clone(),
            patterns: Vec::new(),
            heads: FxHashMap::default(),
            postings: Vec::new(),
            posting_hashes: Vec::new(),
            position_offsets: vec![0],
            active: Vec::new(),
            unique_runs: Vec::new(),
            probe: Vec::new(),
            compact: None,
            posting_dictionary: None,
        };
        this.extend(source)?;
        Ok(this)
    }

    fn extend(&mut self, source: ProtoIndex) -> io::Result<()> {
        let k = self.options.k;
        let first_pattern = self.patterns.len();
        let mut hashes = Vec::new();
        for (id, pattern) in source.patterns.iter().enumerate() {
            let id =
                u32::try_from(first_pattern + id).map_err(|_| invalid("too many proto-MUMs"))?;
            let dna: Vec<_> = (pattern.start..pattern.start + pattern.len)
                .map(|i| source.sequences.as_slice().get_ascii(i))
                .collect();
            let packed = PackedContig::new(&dna, k);
            for reverse in [false, true] {
                if reverse && self.options.forward_only {
                    continue;
                }
                let seq = if reverse {
                    &packed.reverse
                } else {
                    &packed.forward
                };
                hashes.clear();
                if !(3..=64).contains(&k) {
                    NtHasher::<false>::new(k)
                        .hash_kmers_simd(seq.as_slice(), 1)
                        .collect_into(&mut hashes);
                }
                for pos in 0..pattern.len - k + 1 {
                    u32::try_from(self.postings.len())
                        .map_err(|_| invalid("too many proto positions"))?;
                    let offset = if reverse { pattern.len - pos - k } else { pos };
                    let offset = u32::try_from(offset)
                        .map_err(|_| invalid("proto-MUM exceeds position range"))?;
                    if !(3..=64).contains(&k) {
                        self.posting_hashes.push(hashes[pos]);
                    }
                    self.postings.push(SeedPosting {
                        pattern: id,
                        offset,
                        next: NO_POSTING,
                        reverse,
                    });
                }
            }
            self.position_offsets
                .push(self.position_offsets.last().unwrap() + packed.len() - k + 1);
            self.patterns.push(packed);
        }
        // Extension is also used by reusable union states. Restore every raw
        // candidate there; pruning is local to a stateless all-genome search.
        let count = *self.position_offsets.last().unwrap();
        self.active = vec![u64::MAX; count.div_ceil(64)];
        if count % 64 != 0 {
            *self.active.last_mut().unwrap() = (1u64 << (count % 64)) - 1;
        }
        if (3..=64).contains(&k) {
            self.posting_dictionary = postings::PostingDictionary::build(self)?;
        }
        self.rebuild();
        Ok(())
    }

    fn active_positions(&self) -> usize {
        self.active.iter().map(|w| w.count_ones() as usize).sum()
    }

    fn mark_present(&self, present: &mut [u64], fragment: Fragment) {
        let start = self.position_offsets[fragment.hit.pattern as usize] + fragment.begin;
        let end = self.position_offsets[fragment.hit.pattern as usize] + fragment.end
            - self.options.k
            + 1;
        mark_position_range(present, start, end);
    }

    fn retain_present(&mut self, present: &[u64]) -> bool {
        let mut changed = false;
        for (active, &seen) in self.active.iter_mut().zip(present) {
            let new = *active & seen;
            changed |= new != *active;
            *active = new;
        }
        if changed {
            self.rebuild();
        }
        changed
    }

    fn rebuild(&mut self) {
        if self.compact.is_some() {
            return;
        }
        if (3..=64).contains(&self.options.k) && self.postings.is_empty() {
            self.heads.clear();
            self.unique_runs.clear();
            self.probe.clear();
            return;
        }
        if let Some(dictionary) = &mut self.posting_dictionary {
            dictionary.heads.fill(NO_POSTING);
            for (id, posting) in self.postings.iter_mut().enumerate() {
                let bit = self.position_offsets[posting.pattern as usize] + posting.offset as usize;
                posting.next = NO_POSTING;
                if self.active[bit / 64] & (1u64 << (bit % 64)) != 0 {
                    let head = &mut dictionary.heads[dictionary.slots[id] as usize];
                    posting.next = *head;
                    *head = id as u32;
                }
            }
            self.heads.clear();
            self.probe.clear();
            self.unique_runs = vec![0; self.postings.len()];
            for &head in &dictionary.heads {
                if head != NO_POSTING && self.postings[head as usize].next == NO_POSTING {
                    self.unique_runs[head as usize] = 1;
                }
            }
            self.extend_unique_runs();
            return;
        }
        self.heads = FxHashMap::with_capacity_and_hasher(
            self.active_positions()
                .saturating_mul(if self.options.forward_only { 1 } else { 2 }),
            Default::default(),
        );
        for (id, (posting, &hash)) in self
            .postings
            .iter_mut()
            .zip(&self.posting_hashes)
            .enumerate()
        {
            let bit = self.position_offsets[posting.pattern as usize] + posting.offset as usize;
            posting.next = NO_POSTING;
            if self.active[bit / 64] & (1u64 << (bit % 64)) != 0 {
                posting.next = self.heads.insert(hash, id as u32).unwrap_or(NO_POSTING);
            }
        }
        // A unique hash has exactly one possible proto position. After an
        // exact extension, its following unique positions cannot discover any
        // new alignment. Stop at repeats, collisions, proto/strand boundaries,
        // or the end of the verified range. Rebuild after additions: a new
        // proto may invalidate an old hash's uniqueness.
        self.unique_runs = vec![0; self.postings.len()];
        for &head in self.heads.values() {
            if self.postings[head as usize].next == NO_POSTING {
                self.unique_runs[head as usize] = 1;
            }
        }
        self.extend_unique_runs();
        let bits = (self.heads.len().max(64) * 16).next_power_of_two();
        self.probe = vec![0; bits / 64];
        for &hash in self.heads.keys() {
            let bit = hash as usize & (bits - 1);
            self.probe[bit / 64] |= 1 << (bit % 64);
        }
    }

    fn extend_unique_runs(&mut self) {
        for i in (0..self.postings.len().saturating_sub(1)).rev() {
            if self.unique_runs[i] != 0
                && self.postings[i].pattern == self.postings[i + 1].pattern
                && self.postings[i].reverse == self.postings[i + 1].reverse
            {
                self.unique_runs[i] += self.unique_runs[i + 1];
            }
        }
    }

    // SIMD hashes visit every seed position. A cached exact range consumes all
    // later seed hits on that alignment without extending or emitting it again.
    // For a fixed diagonal, proto positions are monotone during this contig
    // scan, so only its last range needs to be retained.
    fn scan(
        &self,
        query: &PackedContig,
        contig: u32,
        origin: u64,
        mut accept: impl FnMut(Fragment) -> io::Result<()>,
    ) -> io::Result<()> {
        if let Some(compact) = &self.compact {
            return compact.scan(self, query, contig, origin, accept);
        }
        if let Some(dictionary) = &self.posting_dictionary {
            return dictionary.scan(self, query, contig, origin, accept);
        }
        if self.heads.is_empty() {
            return Ok(());
        }
        let k = self.options.k;
        let mut hashes = Vec::new();
        let mut covered: FxHashMap<(u32, bool, i64), (usize, usize)> = FxHashMap::default();
        for run in &query.runs {
            for start in (run.start..run.end).step_by(self.options.chunk_bases) {
                let core_end = start.saturating_add(self.options.chunk_bases).min(run.end);
                let end = core_end.saturating_add(k - 1).min(run.end);
                if end - start < k {
                    continue;
                }
                let seq = query.forward.slice(start..end);
                // The dense seed hash probe already rejects absent candidates.
                // A minimizer prepass duplicates work on shared-genome chunks;
                // directly hashing every k-mer preserves the same exact checks.
                hashes.clear();
                NtHasher::<false>::new(k)
                    .hash_kmers_simd(seq, 1)
                    .collect_into(&mut hashes);
                let mut local = 0;
                while local < hashes.len().min(core_end - start) {
                    let hash = hashes[local];
                    let physical_q = start + local;
                    local += 1;
                    let bit = hash as usize & (self.probe.len() * 64 - 1);
                    if self.probe[bit / 64] & (1 << (bit % 64)) == 0 {
                        continue;
                    }
                    let Some(&head) = self.heads.get(&hash) else {
                        continue;
                    };
                    let unique_run = if self.postings[head as usize].next == NO_POSTING {
                        self.unique_runs[head as usize] as usize
                    } else {
                        0
                    };
                    let mut next = head;
                    while next != NO_POSTING {
                        let posting = &self.postings[next as usize];
                        next = posting.next;
                        let r = posting.offset as usize;
                        let q = if posting.reverse {
                            query.len() - physical_q - k
                        } else {
                            physical_q
                        };
                        let key = (posting.pattern, posting.reverse, q as i64 - r as i64);
                        if let Some(&(a, b)) = covered.get(&key) {
                            if a <= r && r + k <= b {
                                if unique_run > 1 {
                                    let remaining = if posting.reverse {
                                        r - a + 1
                                    } else {
                                        b - r - k + 1
                                    };
                                    local += unique_run.min(remaining) - 1;
                                }
                                continue;
                            }
                        }
                        let proto = &self.patterns[posting.pattern as usize];
                        let (qforward, qbackward, qb) = if posting.reverse {
                            (
                                &query.reverse,
                                &query.forward,
                                query.len() - run.end..query.len() - run.start,
                            )
                        } else {
                            (&query.forward, &query.reverse, run.clone())
                        };
                        // Hashes only select candidates; equality and all range
                        // boundaries are checked against packed DNA exactly.
                        if proto
                            .forward
                            .slice(r..r + k)
                            .cmp_lcp(&qforward.slice(q..q + k))
                            .1
                            != k
                        {
                            continue;
                        }
                        let right = proto
                            .forward
                            .slice(r + k..proto.len())
                            .cmp_lcp(&qforward.slice(q + k..qb.end))
                            .1;
                        let left = proto
                            .reverse
                            .slice(proto.len() - r..proto.len())
                            .cmp_lcp(&qbackward.slice(query.len() - q..query.len() - qb.start))
                            .1;
                        let (begin, end) = (r - left, r + k + right);
                        covered.insert(key, (begin, end));
                        let physical_start = if posting.reverse {
                            physical_q - right
                        } else {
                            physical_q - left
                        };
                        accept(Fragment {
                            begin,
                            end,
                            hit: Hit {
                                pattern: posting.pattern,
                                contig,
                                start: origin + physical_start as u64,
                                reverse: posting.reverse,
                            },
                        })?;
                        if unique_run > 1 {
                            let remaining = if posting.reverse {
                                r - begin + 1
                            } else {
                                end - r - k + 1
                            };
                            local += unique_run.min(remaining) - 1;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn is_empty(&self) -> bool {
        self.active_positions() == 0
    }

    fn seed_count(&self) -> usize {
        if self.compact.is_some() {
            *self.position_offsets.last().unwrap()
        } else if let Some(dictionary) = &self.posting_dictionary {
            dictionary.heads.len() / 2
        } else {
            self.heads.len()
        }
    }
}

fn mark_position_range(words: &mut [u64], start: usize, end: usize) {
    if start >= end {
        return;
    }
    let first = start / 64;
    let last = (end - 1) / 64;
    let left = u64::MAX << (start % 64);
    let right = u64::MAX >> (63 - ((end - 1) % 64));
    if first == last {
        words[first] |= left & right;
    } else {
        words[first] |= left;
        words[first + 1..last].fill(u64::MAX);
        words[last] |= right;
    }
}

fn header_name(parser: &FastxParser<CONFIG>) -> String {
    String::from_utf8_lossy(parser.get_header())
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

fn read_reference(index: &FragmentIndex, input: &Path) -> io::Result<Reference> {
    let mut reference = Reference {
        contigs: Vec::new(),
        meta: GenomeMeta {
            input: input.into(),
            contigs: Vec::new(),
        },
        fragments: (0..index.patterns.len()).map(|_| Vec::new()).collect(),
    };
    let data = FileInput::new(input)?;
    if data.first_byte() != 0 {
        let mut parser = FastxParser::<CONFIG>::from_input(data)?;
        let mut origin = 0;
        while parser.next().is_some() {
            let id =
                u32::try_from(reference.contigs.len()).map_err(|_| invalid("too many contigs"))?;
            let sequence = parser.get_dna_string();
            let packed = PackedContig::new(sequence, index.options.k);
            reference.meta.contigs.push(ContigMeta {
                name: header_name(&parser),
                start: origin,
                len: sequence.len(),
            });
            index.scan(&packed, id, origin, |fragment| {
                reference.fragments[fragment.hit.pattern as usize].push(fragment);
                Ok(())
            })?;
            origin += packed.len() as u64;
            reference.contigs.push(packed);
        }
    }
    for fragments in &mut reference.fragments {
        fragments
            .sort_unstable_by_key(|f| (f.begin, f.end, f.hit.contig, f.hit.start, f.hit.reverse));
    }
    Ok(reference)
}

// Maximal exact intervals in reference coordinates. These are position validity
// ranges, not a single accept/reject bit for a complete proto-MUM.
#[derive(Clone, Copy, Debug)]
struct Record {
    begin: u64,
    end: u64,
    hit: Hit,
}
const RECORD_BYTES: usize = 32;

struct ContextSink {
    shards: Vec<Vec<Record>>,
    reference_bases: u64,
    k: u64,
}

fn context_path(scratch: &Path, genome: usize, shard: usize) -> PathBuf {
    scratch.join(format!("contexts-{genome}-{shard}.zst"))
}

fn sort_encode_records(records: &mut Vec<Record>, bytes: &mut Vec<u8>) {
    records.sort_unstable_by_key(|r| (r.begin, r.end, r.hit.contig, r.hit.start, r.hit.reverse));
    records.dedup_by_key(|r| (r.begin, r.end, r.hit.contig, r.hit.start, r.hit.reverse));
    bytes.clear();
    bytes.reserve(records.len() * RECORD_BYTES);
    for record in records.iter() {
        bytes.extend_from_slice(&record.begin.to_le_bytes());
        bytes.extend_from_slice(&record.end.to_le_bytes());
        bytes.extend_from_slice(&record.hit.start.to_le_bytes());
        bytes.extend_from_slice(&record.hit.contig.to_le_bytes());
        bytes.extend_from_slice(&u32::from(record.hit.reverse).to_le_bytes());
    }
}

// Streaming batch searches use one sequential scratch file per worker/shard,
// rather than one compressed file per genome/shard. Scratch evidence is short
// lived and buffered; reusable state keeps its existing compressed format.
struct ContextWriter {
    shards: Vec<ContextBuffer>,
    bytes: Vec<u8>,
}

struct ContextBuffer {
    path: PathBuf,
    bytes: Vec<u8>,
}

impl ContextBuffer {
    fn flush(&mut self) -> io::Result<()> {
        if !self.bytes.is_empty() {
            std::fs::OpenOptions::new()
                .append(true)
                .open(&self.path)?
                .write_all(&self.bytes)?;
            self.bytes.clear();
        }
        Ok(())
    }
}

fn worker_context_path(scratch: &Path, worker: usize, shard: usize) -> PathBuf {
    scratch.join(format!("worker-{worker}-contexts-{shard}.bin"))
}

impl ContextWriter {
    fn new(scratch: &Path, worker: usize, shards: usize) -> io::Result<Self> {
        Ok(Self {
            shards: (0..shards)
                .map(|shard| {
                    let path = worker_context_path(scratch, worker, shard);
                    File::create(&path)?;
                    Ok(ContextBuffer {
                        path,
                        bytes: Vec::with_capacity(64 * 1024),
                    })
                })
                .collect::<io::Result<_>>()?,
            bytes: Vec::new(),
        })
    }

    fn finish(self) -> io::Result<()> {
        for mut shard in self.shards {
            shard.flush()?;
        }
        Ok(())
    }
}

impl ContextSink {
    fn new(shards: usize, reference: &Reference, k: usize) -> Self {
        Self {
            shards: (0..shards).map(|_| Vec::new()).collect(),
            reference_bases: reference
                .meta
                .contigs
                .iter()
                .map(|c| c.len as u64)
                .sum::<u64>()
                .max(1),
            k: k as u64,
        }
    }
    fn shard(&self, start: u64) -> usize {
        ((start as u128 * self.shards.len() as u128) / self.reference_bases as u128) as usize
    }
    fn push(&mut self, record: Record) {
        let first = self.shard(record.begin);
        let last = self.shard(record.end - self.k);
        for shard in first..=last {
            self.shards[shard].push(record);
        }
    }
    fn finish(mut self, scratch: &Path, genome: usize) -> io::Result<()> {
        let mut compressor = zstd::bulk::Compressor::new(1)?;
        let mut bytes = Vec::new();
        for (shard, records) in self.shards.iter_mut().enumerate() {
            sort_encode_records(records, &mut bytes);
            std::fs::write(
                context_path(scratch, genome, shard),
                compressor.compress(&bytes)?,
            )?;
        }
        Ok(())
    }

    fn finish_to(mut self, writer: &mut ContextWriter, genome: usize) -> io::Result<()> {
        for (shard, records) in self.shards.iter_mut().enumerate() {
            if records.is_empty() {
                continue;
            }
            sort_encode_records(records, &mut writer.bytes);
            let mut header = [0u8; 16];
            header[..8].copy_from_slice(&(genome as u64).to_le_bytes());
            header[8..].copy_from_slice(&(records.len() as u64).to_le_bytes());
            let target = &mut writer.shards[shard];
            target.bytes.extend_from_slice(&header);
            target.bytes.extend_from_slice(&writer.bytes);
            // Bound open descriptors by active workers rather than workers*shards.
            if target.bytes.len() >= 64 * 1024 {
                target.flush()?;
            }
        }
        Ok(())
    }
}

fn seed_at(fragment: Fragment, pos: usize, k: usize) -> Hit {
    Hit {
        start: fragment.hit.start
            + if fragment.hit.reverse {
                (fragment.end - pos - k) as u64
            } else {
                (pos - fragment.begin) as u64
            },
        ..fragment.hit
    }
}

fn consume_fragment(
    reference: &Reference,
    fragment: Fragment,
    query: &PackedContig,
    origin: u64,
    k: usize,
    cache: &mut AlignmentCache,
    stats: &mut ScanStats,
    sink: &mut ContextSink,
) -> io::Result<()> {
    for &source in &reference.fragments[fragment.hit.pattern as usize] {
        let pos = source.begin.max(fragment.begin);
        if pos + k > source.end.min(fragment.end) {
            continue;
        }
        let anchor = Anchor {
            hit: seed_at(source, pos, k),
            len: k,
        };
        let hit = seed_at(fragment, pos, k);
        let before = stats.alignments;
        let (left, right) =
            alignment_context(reference, &anchor, query, hit, origin, cache, stats)?;
        if before == stats.alignments {
            continue;
        }
        let begin = anchor.hit.start - if anchor.hit.reverse { right } else { left };
        let end = anchor.hit.start + k as u64 + if anchor.hit.reverse { left } else { right };
        sink.push(Record {
            begin,
            end,
            hit: Hit {
                start: hit.start - if hit.reverse { right } else { left },
                reverse: hit.reverse != anchor.hit.reverse,
                ..hit
            },
        });
        stats.contexts += 1;
    }
    Ok(())
}

fn reference_contexts(
    reference: &Reference,
    shards: usize,
    k: usize,
    writer: &mut ContextWriter,
) -> io::Result<ScanStats> {
    let mut sink = ContextSink::new(shards, reference, k);
    let mut stats = ScanStats {
        genomes: 1,
        ..Default::default()
    };
    // Group fragment occurrences by contig to reuse the same exact alignment
    // cache across all proto-MUMs, including repeat and reverse occurrences.
    for (id, query) in reference.contigs.iter().enumerate() {
        let mut cache = AlignmentCache::default();
        for fragments in &reference.fragments {
            for &fragment in fragments {
                if fragment.hit.contig as usize == id {
                    consume_fragment(
                        reference,
                        fragment,
                        query,
                        reference.meta.contigs[id].start,
                        k,
                        &mut cache,
                        &mut stats,
                        &mut sink,
                    )?;
                }
            }
        }
    }
    sink.finish_to(writer, 0)?;
    Ok(stats)
}

#[cfg(test)]
fn parallel_reference_contexts(
    reference: &Reference,
    shards: usize,
    k: usize,
    workers: usize,
    writer: &mut ContextWriter,
) -> io::Result<ScanStats> {
    let mut grouped: Vec<Vec<Fragment>> =
        (0..reference.contigs.len()).map(|_| Vec::new()).collect();
    for &fragment in reference.fragments.iter().flatten() {
        grouped[fragment.hit.contig as usize].push(fragment);
    }
    let next: Vec<_> = grouped.iter().map(|_| AtomicUsize::new(0)).collect();
    let workers = workers.min(
        grouped
            .iter()
            .map(Vec::len)
            .sum::<usize>()
            .div_ceil(256)
            .max(1),
    );
    let parts = thread::scope(|scope| -> io::Result<Vec<(ContextSink, ScanStats)>> {
        let mut handles = Vec::new();
        for _ in 0..workers {
            let grouped = &grouped;
            let next = &next;
            handles.push(scope.spawn(move || -> io::Result<_> {
                let mut sink = ContextSink::new(shards, reference, k);
                let mut stats = ScanStats::default();
                for (id, fragments) in grouped.iter().enumerate() {
                    let mut cache = AlignmentCache::default();
                    loop {
                        let start = next[id].fetch_add(256, Ordering::Relaxed);
                        if start >= fragments.len() {
                            break;
                        }
                        for &fragment in &fragments[start..(start + 256).min(fragments.len())] {
                            consume_fragment(
                                reference,
                                fragment,
                                &reference.contigs[id],
                                reference.meta.contigs[id].start,
                                k,
                                &mut cache,
                                &mut stats,
                                &mut sink,
                            )?;
                        }
                    }
                }
                Ok((sink, stats))
            }));
        }
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .map_err(|_| io::Error::other("reference worker panicked"))?
            })
            .collect()
    })?;
    let mut stats = ScanStats {
        genomes: 1,
        ..Default::default()
    };
    for (sink, part) in parts {
        sink.finish_to(writer, 0)?;
        stats.add(part);
    }
    Ok(stats)
}

fn scan_query(
    index: &FragmentIndex,
    reference: &Reference,
    restriction: Option<&restriction::RestrictionIndex>,
    input: &Path,
    id: usize,
    shards: usize,
    writer: &mut ContextWriter,
    mut present: Option<&mut [u64]>,
) -> io::Result<(GenomeMeta, ScanStats)> {
    let mut sink = ContextSink::new(shards, reference, index.options.k);
    let mut stats = ScanStats {
        genomes: 1,
        ..Default::default()
    };
    let mut meta = GenomeMeta {
        input: input.into(),
        contigs: Vec::new(),
    };
    let data = FileInput::new(input)?;
    if data.first_byte() != 0 {
        let mut parser = FastxParser::<CONFIG>::from_input(data)?;
        let mut origin = 0;
        while parser.next().is_some() {
            let contig =
                u32::try_from(meta.contigs.len()).map_err(|_| invalid("too many contigs"))?;
            let sequence = parser.get_dna_string();
            meta.contigs.push(ContigMeta {
                name: header_name(&parser),
                start: origin,
                len: sequence.len(),
            });
            if index.is_empty() {
                origin += sequence.len() as u64;
                continue;
            }
            let packed = PackedContig::new(sequence, index.options.k);
            let mut cache = AlignmentCache::default();
            if let Some(restriction) = restriction {
                restriction.scan(
                    index,
                    reference,
                    &packed,
                    contig,
                    origin,
                    &mut stats,
                    &mut sink,
                    present.as_deref_mut(),
                )?;
            } else {
                index.scan(&packed, contig, origin, |fragment| {
                    if let Some(bits) = present.as_deref_mut() {
                        index.mark_present(bits, fragment);
                    }
                    consume_fragment(
                        reference,
                        fragment,
                        &packed,
                        origin,
                        index.options.k,
                        &mut cache,
                        &mut stats,
                        &mut sink,
                    )
                })?;
            }
            origin += packed.len() as u64;
        }
    }
    sink.finish_to(writer, id)?;
    eprintln!(
        "[genome_scan]\t{}\t{} validity ranges\t{} alignments\t{} cached",
        id + 1,
        stats.contexts,
        stats.alignments,
        stats.cache_hits
    );
    Ok((meta, stats))
}

#[derive(Clone, Copy)]
struct Context {
    genome: usize,
    record: Record,
}
struct VerifiedMatch {
    len: u64,
    hits: Vec<Hit>,
}

fn interval_hit(record: Record, begin: u64, end: u64) -> Hit {
    Hit {
        start: record.hit.start
            + if record.hit.reverse {
                record.end - end
            } else {
                begin - record.begin
            },
        ..record.hit
    }
}

// Equal oriented mappings describe the same physical occurrence throughout
// their overlap. Use widened arithmetic for valid u64 coordinates near the limit.
fn same_mapping(a: Record, b: Record) -> bool {
    a.hit.contig == b.hit.contig
        && a.hit.reverse == b.hit.reverse
        && if a.hit.reverse {
            a.hit.start as u128 + a.end as u128 == b.hit.start as u128 + b.end as u128
        } else {
            a.hit.start as u128 + b.begin as u128 == b.hit.start as u128 + a.begin as u128
        }
}

const NO_CONTEXT: usize = usize::MAX;

// At one physical position there can be at most two oriented mappings (the
// palindrome case). Therefore the three furthest distinct mappings suffice:
// if all three reach a candidate's end, at least two physical copies exist.
// Discarded shorter mappings can never change that decision. Later records can
// promote their mapping again if they extend beyond one of the retained ends.
fn retain_occurrence(top: &mut [usize; 3], contexts: &[Context], candidate: usize) {
    let record = contexts[candidate].record;
    if top[0] != NO_CONTEXT && contexts[top[0]].record.end <= record.begin {
        *top = [candidate, NO_CONTEXT, NO_CONTEXT];
        return;
    }
    let existing = top
        .iter()
        .position(|&id| id != NO_CONTEXT && same_mapping(contexts[id].record, record));
    let mut slot = if let Some(slot) = existing {
        if contexts[top[slot]].record.end >= record.end {
            return;
        }
        slot
    } else {
        if top[2] != NO_CONTEXT && contexts[top[2]].record.end >= record.end {
            return;
        }
        2
    };
    while slot > 0
        && (top[slot - 1] == NO_CONTEXT || contexts[top[slot - 1]].record.end < record.end)
    {
        top[slot] = top[slot - 1];
        slot -= 1;
    }
    top[slot] = candidate;
}

// As starts advance, each genome contributes its furthest surviving end.
// Uniqueness uses bounded occurrence evidence, without rescanning old contexts.
fn reduce_ranges(
    contexts: &mut [Context],
    k: u64,
    genomes: usize,
    shard: usize,
    shards: usize,
    bases: u64,
) -> Vec<VerifiedMatch> {
    contexts.sort_unstable_by_key(|c| c.record.begin);
    let mut best = vec![[NO_CONTEXT; 3]; genomes];
    let mut present = 0;
    let mut previous_end = 0;
    let mut cursor = 0;
    let mut output = Vec::new();
    while cursor < contexts.len() {
        let begin = contexts[cursor].record.begin;
        while cursor < contexts.len() && contexts[cursor].record.begin == begin {
            let top = &mut best[contexts[cursor].genome];
            if top[0] == NO_CONTEXT {
                present += 1;
            }
            retain_occurrence(top, contexts, cursor);
            cursor += 1;
        }
        if present != genomes {
            continue;
        }
        let end = best
            .iter()
            .map(|top| contexts[top[0]].record.end)
            .min()
            .unwrap();
        if end <= previous_end {
            continue;
        }
        previous_end = end;
        if end < begin + k || ((begin as u128 * shards as u128) / bases as u128) as usize != shard {
            continue;
        }
        let mut hits = Vec::with_capacity(genomes);
        let mut ambiguous = false;
        for top in &best {
            let mut chosen = interval_hit(contexts[top[0]].record, begin, end);
            for &id in &top[1..] {
                if id == NO_CONTEXT || contexts[id].record.end < end {
                    break;
                }
                let candidate = interval_hit(contexts[id].record, begin, end);
                if candidate.contig != chosen.contig || candidate.start != chosen.start {
                    ambiguous = true;
                    break;
                }
                // A reverse-complement palindrome is one physical occurrence.
                if !candidate.reverse {
                    chosen = candidate;
                }
            }
            if ambiguous {
                break;
            }
            hits.push(chosen);
        }
        if !ambiguous {
            output.push(VerifiedMatch {
                len: end - begin,
                hits,
            });
        }
    }
    output
}

fn read_worker_contexts(
    input: &mut impl Read,
    contexts: &mut Vec<Context>,
    bytes: &mut Vec<u8>,
    genomes: usize,
) -> io::Result<()> {
    loop {
        let mut header = [0u8; 16];
        if input.read(&mut header[..1])? == 0 {
            return Ok(());
        }
        input.read_exact(&mut header[1..])?;
        let genome = usize::try_from(u64::from_le_bytes(header[..8].try_into().unwrap()))
            .map_err(|_| invalid("invalid context genome"))?;
        let count = usize::try_from(u64::from_le_bytes(header[8..].try_into().unwrap()))
            .map_err(|_| invalid("invalid context batch size"))?;
        if genome >= genomes {
            return Err(invalid("context genome outside input collection"));
        }
        let size = count
            .checked_mul(RECORD_BYTES)
            .ok_or_else(|| invalid("context batch too large"))?;
        bytes.clear();
        bytes
            .try_reserve(size)
            .map_err(|_| invalid("context batch too large"))?;
        bytes.resize(size, 0);
        input.read_exact(bytes)?;
        for record in bytes.chunks_exact(RECORD_BYTES) {
            contexts.push(Context {
                genome,
                record: Record {
                    begin: u64::from_le_bytes(record[..8].try_into().unwrap()),
                    end: u64::from_le_bytes(record[8..16].try_into().unwrap()),
                    hit: Hit {
                        pattern: 0,
                        start: u64::from_le_bytes(record[16..24].try_into().unwrap()),
                        contig: u32::from_le_bytes(record[24..28].try_into().unwrap()),
                        reverse: u32::from_le_bytes(record[28..32].try_into().unwrap()) != 0,
                    },
                },
            });
        }
    }
}

fn reduce_shard(
    reference: &Reference,
    scratch: &Path,
    shard: usize,
    shards: usize,
    workers: usize,
    genomes: usize,
    k: usize,
) -> io::Result<Vec<VerifiedMatch>> {
    let mut contexts = Vec::new();
    let mut bytes = Vec::new();
    for worker in 0..workers {
        let mut input = BufReader::with_capacity(
            64 * 1024,
            File::open(worker_context_path(scratch, worker, shard))?,
        );
        read_worker_contexts(&mut input, &mut contexts, &mut bytes, genomes)?;
    }
    let bases = reference
        .meta
        .contigs
        .iter()
        .map(|c| c.len as u64)
        .sum::<u64>()
        .max(1);
    let output = reduce_ranges(&mut contexts, k as u64, genomes, shard, shards, bases);
    Ok(output)
}
// For a fixed reference contig, orientation and diagonal, a maximal exact
// alignment interval can serve every seed inside it. No probabilistic LCEs.
type AlignmentCache = FxHashMap<(u32, bool, i64), AlignmentIntervals>;
#[cfg(test)]
type MonotonicAlignmentCache = FxHashMap<(u32, bool, i64), LastAlignmentInterval>;

trait AlignmentHistory: Default {
    fn preceding(&self, position: usize) -> Option<(usize, usize)>;
    fn insert(&mut self, begin: usize, end: usize);
}

#[cfg(test)]
#[derive(Default)]
struct LastAlignmentInterval {
    begin: usize,
    end: usize,
}
#[cfg(test)]
impl AlignmentHistory for LastAlignmentInterval {
    fn preceding(&self, position: usize) -> Option<(usize, usize)> {
        (self.begin <= position).then_some((self.begin, self.end))
    }
    fn insert(&mut self, begin: usize, end: usize) {
        self.begin = begin;
        self.end = end;
    }
}

#[derive(Default)]
enum AlignmentIntervals {
    #[default]
    Empty,
    Single(usize, usize),
    Multiple(BTreeMap<usize, usize>),
}

impl AlignmentHistory for AlignmentIntervals {
    fn preceding(&self, position: usize) -> Option<(usize, usize)> {
        match self {
            Self::Empty => None,
            Self::Single(begin, end) => (*begin <= position).then_some((*begin, *end)),
            Self::Multiple(ranges) => ranges.range(..=position).next_back().map(|(&a, &b)| (a, b)),
        }
    }
    fn insert(&mut self, begin: usize, end: usize) {
        match self {
            Self::Empty => *self = Self::Single(begin, end),
            Self::Single(old_begin, old_end) if *old_begin == begin => *old_end = end,
            Self::Single(old_begin, old_end) => {
                *self = Self::Multiple(BTreeMap::from([(*old_begin, *old_end), (begin, end)]));
            }
            Self::Multiple(ranges) => {
                ranges.insert(begin, end);
            }
        }
    }
}

#[inline(always)]
fn alignment_context(
    reference: &Reference,
    anchor: &Anchor,
    query: &PackedContig,
    hit: Hit,
    origin: u64,
    cache: &mut AlignmentCache,
    stats: &mut ScanStats,
) -> io::Result<(u64, u64)> {
    alignment_context_impl(reference, anchor, query, hit, origin, cache, stats)
}

// A fixed-length compact scan visits query seed starts monotonically. For a
// fixed reference contig, strand and diagonal, distinct exact maximal intervals
// are disjoint and visited in order. If the latest interval cannot cover the
// next seed, no older interval can. Arbitrarily ordered fragment scans retain
// the full history through alignment_context instead.
#[cfg(test)]
#[inline(always)]
fn alignment_context_monotonic(
    reference: &Reference,
    anchor: &Anchor,
    query: &PackedContig,
    hit: Hit,
    origin: u64,
    cache: &mut MonotonicAlignmentCache,
    stats: &mut ScanStats,
) -> io::Result<(u64, u64)> {
    alignment_context_impl(reference, anchor, query, hit, origin, cache, stats)
}

fn alignment_context_impl<History: AlignmentHistory>(
    reference: &Reference,
    anchor: &Anchor,
    query: &PackedContig,
    hit: Hit,
    origin: u64,
    cache: &mut FxHashMap<(u32, bool, i64), History>,
    stats: &mut ScanStats,
) -> io::Result<(u64, u64)> {
    let r = (anchor.hit.start - reference.meta.contigs[anchor.hit.contig as usize].start) as usize;
    let physical_q = (hit.start - origin) as usize;
    let reverse = anchor.hit.reverse != hit.reverse;
    let q = if reverse {
        query.len() - physical_q - anchor.len
    } else {
        physical_q
    };
    let key = (anchor.hit.contig, reverse, q as i64 - r as i64);
    let intervals = cache.entry(key).or_default();
    if let Some((begin, end)) = intervals.preceding(r) {
        if end >= r + anchor.len {
            stats.cache_hits += 1;
            let (left, right) = ((r - begin) as u64, (end - r - anchor.len) as u64);
            return Ok(if anchor.hit.reverse {
                (right, left)
            } else {
                (left, right)
            });
        }
    }
    let (left, right) = alignment_context_uncached(reference, anchor, query, hit, origin, stats)?;
    let (physical_left, physical_right) = if anchor.hit.reverse {
        (right, left)
    } else {
        (left, right)
    };
    intervals.insert(
        r - physical_left as usize,
        r + anchor.len + physical_right as usize,
    );
    Ok((left, right))
}

// Exact bounded extension, shared by cached fragment scans and the compact
// scanner. Its direct-mapped recent cache is bounded independently of input
// size; collisions cause recomputation and scratch dedup removes duplicate
// contexts. No unbounded per-diagonal map is needed on this path.
#[inline(always)]
fn alignment_context_uncached(
    reference: &Reference,
    anchor: &Anchor,
    query: &PackedContig,
    hit: Hit,
    origin: u64,
    stats: &mut ScanStats,
) -> io::Result<(u64, u64)> {
    let rcontig = anchor.hit.contig as usize;
    let rseq = &reference.contigs[rcontig];
    let r = (anchor.hit.start - reference.meta.contigs[rcontig].start) as usize;
    let physical_q = (hit.start - origin) as usize;
    let reverse = anchor.hit.reverse != hit.reverse;
    let q = if reverse {
        query.len() - physical_q - anchor.len
    } else {
        physical_q
    };
    let rb = rseq.bounds(r, anchor.len)?;
    let qb = query.bounds(physical_q, anchor.len)?;
    let qb = if reverse {
        query.len() - qb.end..query.len() - qb.start
    } else {
        qb
    };
    let (qforward, qbackward) = if reverse {
        (&query.reverse, &query.forward)
    } else {
        (&query.forward, &query.reverse)
    };
    // packed-seq compares 29 bases per machine word and stops exactly at the
    // first mismatch. Explicit run bounds prevent extending across ambiguity.
    let right = rseq
        .forward
        .slice(r + anchor.len..rb.end)
        .cmp_lcp(&qforward.slice(q + anchor.len..qb.end))
        .1;
    let left = rseq
        .reverse
        .slice(rseq.len() - r..rseq.len() - rb.start)
        .cmp_lcp(&qbackward.slice(query.len() - q..query.len() - qb.start))
        .1;
    stats.alignments += 1;
    stats.extended_bases += (left + right) as u64;
    Ok(if anchor.hit.reverse {
        (right as u64, left as u64)
    } else {
        (left as u64, right as u64)
    })
}

fn write_output_shard(
    mut matches: Vec<VerifiedMatch>,
    genomes: &[GenomeMeta],
    scratch: &Path,
    shard: usize,
) -> io::Result<(usize, usize)> {
    for item in &mut matches {
        if item.hits[0].reverse {
            for hit in &mut item.hits {
                hit.reverse = !hit.reverse;
            }
        }
    }
    matches.sort_unstable_by(|a, b| {
        (a.hits[0].start, a.len)
            .cmp(&(b.hits[0].start, b.len))
            .then_with(|| {
                a.hits
                    .iter()
                    .map(|h| (h.contig, h.start, h.reverse))
                    .cmp(b.hits.iter().map(|h| (h.contig, h.start, h.reverse)))
            })
    });
    let total = matches.len();
    matches.dedup_by_key(|m| (m.hits[0].start, m.len));
    let count = matches.len();
    let file = File::create(scratch.join(format!("output-{shard}.tsv.zst")))?;
    let mut writer = zstd::stream::write::Encoder::new(BufWriter::new(file), 1)?;
    if shard == 0 {
        writeln!(
            writer,
            "# Coordinates: zero-based forward-strand starts within each contig; input genome order"
        )?;
        writeln!(writer, "length\tstarts\tstrands\tcontigs")?;
    }
    for item in matches {
        write!(writer, "{}\t", item.len)?;
        for (i, hit) in item.hits.iter().enumerate() {
            if i != 0 {
                write!(writer, ",")?;
            }
            write!(
                writer,
                "{}",
                hit.start - genomes[i].contigs[hit.contig as usize].start
            )?;
        }
        write!(writer, "\t")?;
        for (i, hit) in item.hits.iter().enumerate() {
            if i != 0 {
                write!(writer, ",")?;
            }
            write!(writer, "{}", if hit.reverse { '-' } else { '+' })?;
        }
        write!(writer, "\t")?;
        for (i, hit) in item.hits.iter().enumerate() {
            if i != 0 {
                write!(writer, ",")?;
            }
            write!(writer, "{}", hit.contig)?;
        }
        writeln!(writer)?;
    }
    writer.finish()?.flush()?;
    Ok((count, total - count))
}

#[derive(Debug, Default)]
pub struct LocateTimings {
    pub index: f64,
    pub reference_scan: f64,
    pub scan_extend: f64,
    pub reduce: f64,
    pub output: f64,
    pub finalize: f64,
}

impl LocateTimings {
    pub fn phases(&self) -> [(&'static str, f64); 6] {
        [
            ("proto_index", self.index),
            ("reference_scan", self.reference_scan),
            ("parallel_genome_scan_extend", self.scan_extend),
            ("parallel_context_reduce", self.reduce),
            ("parallel_dedup_output", self.output),
            ("locator_finalize", self.finalize),
        ]
    }
}

#[derive(Debug, Default)]
pub struct LocateStats {
    pub proto_mums: usize,
    pub verified_mums: usize,
    pub merged_duplicates: usize,
    pub source_genomes_read: usize,
    pub alignment_computations: usize,
    pub alignment_cache_hits: usize,
    pub timings: LocateTimings,
}

/// One source pass per genome after GGCAT seed construction. Genome workers
/// produce exact position validity ranges; independent reducer and output
/// shards use all requested workers without reopening any genome sequence.
pub fn locate_mums(
    patterns: &Path,
    inputs: &[PathBuf],
    output: &Path,
    scratch_parent: Option<&Path>,
    options: LocateOptions,
) -> io::Result<LocateStats> {
    locate_mums_with_batch(patterns, inputs, output, scratch_parent, options, 256)
}

/// Locate from the complete shared de Bruijn graph (one canonical position per
/// k-mer). SSHash is a proto dictionary; source documents are scanned once.
/// Arbitrary/redundant proto sets must use `locate_mums` instead.
pub fn locate_graph_mums(
    patterns: &Path,
    inputs: &[PathBuf],
    output: &Path,
    scratch_parent: Option<&Path>,
    options: LocateOptions,
) -> io::Result<LocateStats> {
    locate_impl(patterns, inputs, output, scratch_parent, options, 256, true)
}

fn genome_scan_order(genomes: usize, first_batch: usize) -> Vec<usize> {
    if first_batch == 0 {
        return (1..genomes).collect();
    }
    let initial = first_batch.min(genomes.saturating_sub(1));
    let mut order = Vec::with_capacity(genomes.saturating_sub(1));
    let mut seen = vec![false; genomes];
    for i in 0..initial {
        let id = if initial == 1 {
            genomes - 1
        } else {
            1 + ((i as u128 * (genomes - 2) as u128) / (initial - 1) as u128) as usize
        };
        seen[id] = true;
        order.push(id);
    }
    // Visit lower-quality tail assemblies early on quality-ordered inputs.
    // Original genome IDs and output order remain unchanged.
    order.extend((1..genomes).rev().filter(|&id| !seen[id]));
    order
}

fn locate_mums_with_batch(
    patterns: &Path,
    inputs: &[PathBuf],
    output: &Path,
    scratch_parent: Option<&Path>,
    options: LocateOptions,
    first_batch: usize,
) -> io::Result<LocateStats> {
    locate_impl(
        patterns,
        inputs,
        output,
        scratch_parent,
        options,
        first_batch,
        false,
    )
}

fn locate_impl(
    patterns: &Path,
    inputs: &[PathBuf],
    output: &Path,
    scratch_parent: Option<&Path>,
    options: LocateOptions,
    first_batch: usize,
    graph: bool,
) -> io::Result<LocateStats> {
    if inputs.is_empty() {
        return Err(invalid("at least one genome is required"));
    }
    let started = Instant::now();
    let source = ProtoIndex::from_fasta(patterns, options)?;
    let mut index = if graph {
        FragmentIndex::from_graph(source, scratch_parent)?
    } else {
        FragmentIndex::new(source)?
    };
    // The compact restriction scanner uses immutable reference occurrence
    // evidence, not the mutable hash-posting pruning table. Its final reduction
    // still requires exact contexts from every genome. Avoid collecting unused
    // presence masks so verified unique alignments can skip across proto records.
    let first_batch = if index.compact.is_some() {
        0
    } else {
        first_batch
    };
    let mut timings = LocateTimings {
        index: started.elapsed().as_secs_f64(),
        ..Default::default()
    };
    let workers = index.options.threads;
    let shards = workers.min(index.patterns.len().max(1));
    eprintln!(
        "K2Rmini: indexed {} proto-MUMs ({} seed entries); {} workers; one query pass per genome",
        index.patterns.len(),
        index.seed_count(),
        workers
    );
    let scratch = if let Some(parent) = scratch_parent {
        std::fs::create_dir_all(parent)?;
        tempfile::Builder::new()
            .prefix("k2rlocate-")
            .tempdir_in(parent)?
    } else {
        tempfile::Builder::new().prefix("k2rlocate-").tempdir()?
    };
    let started = Instant::now();
    let mut reference = read_reference(&index, &inputs[0])?;
    let mut largest: Vec<_> = reference
        .fragments
        .iter()
        .enumerate()
        .map(|(id, f)| (f.len(), id))
        .collect();
    largest.sort_unstable_by(|a, b| b.cmp(a));
    eprintln!(
        "[reference_occurrences]\t{:.6}s\t{} fragments\t{} naive occurrence pairs\t{:?} largest (copies, proto)",
        started.elapsed().as_secs_f64(),
        reference.fragments.iter().map(Vec::len).sum::<usize>(),
        reference
            .fragments
            .iter()
            .map(|f| (f.len() as u128).pow(2))
            .sum::<u128>(),
        &largest[..largest.len().min(8)]
    );
    let mut reference_writer = ContextWriter::new(scratch.path(), 0, shards)?;
    let contexts_started = Instant::now();
    let restriction = if index.compact.is_some() {
        Some(restriction::RestrictionIndex::build(&index, &reference)?)
    } else {
        None
    };
    let mut scan_stats = if index.compact.is_some() {
        restriction::RestrictionIndex::identity(
            &reference,
            shards,
            index.options.k,
            &mut reference_writer,
        )?
    } else {
        reference_contexts(&reference, shards, index.options.k, &mut reference_writer)?
    };
    reference_writer.finish()?;
    eprintln!(
        "[reference_contexts]\t{:.6}s\t{} contexts\t{} alignments",
        contexts_started.elapsed().as_secs_f64(),
        scan_stats.contexts,
        scan_stats.alignments
    );
    if first_batch != 0 {
        let mut present = vec![0u64; index.active.len()];
        for &fragment in reference.fragments.iter().flatten() {
            index.mark_present(&mut present, fragment);
        }
        index.retain_present(&present);
    }
    if restriction.is_some() {
        reference.fragments.clear();
    }
    timings.reference_scan = started.elapsed().as_secs_f64();
    let started = Instant::now();
    let query_workers = workers.min(inputs.len().saturating_sub(1));
    let mut writers = (1..=query_workers)
        .map(|worker| ContextWriter::new(scratch.path(), worker, shards))
        .collect::<io::Result<Vec<_>>>()?;
    let order = genome_scan_order(inputs.len(), first_batch);
    let mut parts = Vec::with_capacity(order.len());
    let mut begin = 0;
    let mut batch_size = if first_batch == 0 {
        order.len().max(1)
    } else {
        first_batch
    };
    while begin < order.len() {
        let end = begin.saturating_add(batch_size).min(order.len());
        let batch = &order[begin..end];
        let track = first_batch != 0 && end < order.len() && !index.is_empty();
        let next = AtomicUsize::new(0);
        let results = thread::scope(
            |scope| -> io::Result<Vec<(Vec<(usize, GenomeMeta, ScanStats)>, Vec<u64>)>> {
                let mut handles = Vec::new();
                for writer in writers.iter_mut().take(batch.len()) {
                    let next = &next;
                    let index = &index;
                    let reference = &reference;
                    let restriction = restriction.as_ref();
                    handles.push(scope.spawn(move || -> io::Result<_> {
                        let mut result = Vec::new();
                        let mut surviving = if track {
                            index.active.clone()
                        } else {
                            Vec::new()
                        };
                        let mut present = vec![0; surviving.len()];
                        loop {
                            let cursor = next.fetch_add(1, Ordering::Relaxed);
                            if cursor >= batch.len() {
                                break;
                            }
                            let id = batch[cursor];
                            present.fill(0);
                            let (meta, stats) = scan_query(
                                index,
                                reference,
                                restriction,
                                &inputs[id],
                                id,
                                shards,
                                writer,
                                if track { Some(&mut present) } else { None },
                            )?;
                            if track {
                                for (keep, &seen) in surviving.iter_mut().zip(&present) {
                                    *keep &= seen;
                                }
                            }
                            result.push((id, meta, stats));
                        }
                        Ok((result, surviving))
                    }));
                }
                handles
                    .into_iter()
                    .map(|h| {
                        h.join()
                            .map_err(|_| io::Error::other("genome worker panicked"))?
                    })
                    .collect()
            },
        )?;
        let prune_started = Instant::now();
        let before = index.active_positions();
        let mut surviving = if track {
            index.active.clone()
        } else {
            Vec::new()
        };
        for (result, seen) in results {
            parts.extend(result);
            if track {
                for (keep, observed) in surviving.iter_mut().zip(seen) {
                    *keep &= observed;
                }
            }
        }
        if track {
            index.retain_present(&surviving);
            eprintln!(
                "[position_prune]\t{} genomes\t{} -> {} seed positions\t{} hashes\t{:.6}s",
                end + 1,
                before,
                index.active_positions(),
                index.seed_count(),
                prune_started.elapsed().as_secs_f64()
            );
        }
        begin = end;
        batch_size = batch_size.saturating_mul(4);
    }
    for writer in writers {
        writer.finish()?;
    }
    timings.scan_extend = started.elapsed().as_secs_f64();
    let mut metas: Vec<Option<GenomeMeta>> = (0..inputs.len()).map(|_| None).collect();
    for (id, meta, stats) in parts {
        metas[id] = Some(meta);
        scan_stats.add(stats);
    }
    let started = Instant::now();
    let k = index.options.k;
    let mut reductions = thread::scope(|scope| -> io::Result<Vec<Vec<VerifiedMatch>>> {
        let mut handles = Vec::new();
        for shard in 0..shards {
            let reference = &reference;
            let scratch = scratch.path();
            handles.push(scope.spawn(move || {
                reduce_shard(
                    reference,
                    scratch,
                    shard,
                    shards,
                    query_workers + 1,
                    inputs.len(),
                    k,
                )
            }));
        }
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .map_err(|_| io::Error::other("context reducer panicked"))?
            })
            .collect()
    })?;
    if index.compact.is_some() {
        for matches in &mut reductions {
            matches.retain(|item| {
                restriction
                    .as_ref()
                    .unwrap()
                    .unique(&index, &reference, item)
            });
        }
    }
    reference.contigs.clear();
    timings.reduce = started.elapsed().as_secs_f64();
    let reference_bases = reference
        .meta
        .contigs
        .iter()
        .map(|c| c.len as u64)
        .sum::<u64>()
        .max(1);
    metas[0] = Some(reference.meta);
    let genomes: Vec<_> = metas.into_iter().map(Option::unwrap).collect();
    let mut stats = LocateStats {
        proto_mums: index.patterns.len(),
        source_genomes_read: scan_stats.genomes,
        alignment_computations: scan_stats.alignments,
        alignment_cache_hits: scan_stats.cache_hits,
        ..Default::default()
    };
    let started = Instant::now();
    let mut output_shards: Vec<Vec<VerifiedMatch>> = (0..shards).map(|_| Vec::new()).collect();
    for matches in reductions {
        for item in matches {
            let shard =
                ((item.hits[0].start as u128 * shards as u128) / reference_bases as u128) as usize;
            output_shards[shard.min(shards - 1)].push(item);
        }
    }
    let output_stats = thread::scope(|scope| -> io::Result<Vec<(usize, usize)>> {
        let mut handles = Vec::new();
        for (shard, matches) in output_shards.into_iter().enumerate() {
            let genomes = &genomes;
            let scratch = scratch.path();
            handles.push(scope.spawn(move || write_output_shard(matches, genomes, scratch, shard)));
        }
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .map_err(|_| io::Error::other("output worker panicked"))?
            })
            .collect()
    })?;
    for (count, duplicates) in output_stats {
        stats.verified_mums += count;
        stats.merged_duplicates += duplicates;
    }
    timings.output = started.elapsed().as_secs_f64();
    let started = Instant::now();
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut result = tempfile::NamedTempFile::new_in(parent)?;
    for shard in 0..shards {
        io::copy(
            &mut File::open(scratch.path().join(format!("output-{shard}.tsv.zst")))?,
            result.as_file_mut(),
        )?;
    }
    result.as_file_mut().flush()?;
    let mut mapping_name = output.as_os_str().to_os_string();
    mapping_name.push(".contigs.tsv.zst");
    let mut mapping = tempfile::NamedTempFile::new_in(parent)?;
    {
        let mut writer =
            zstd::stream::write::Encoder::new(BufWriter::new(mapping.as_file_mut()), 1)?;
        writeln!(writer, "genome\tinput\tcontig\tname\tlength")?;
        for (id, genome) in genomes.iter().enumerate() {
            for (contig_id, contig) in genome.contigs.iter().enumerate() {
                writeln!(
                    writer,
                    "{id}\t{}\t{contig_id}\t{}\t{}",
                    genome.input.display(),
                    contig.name,
                    contig.len
                )?;
            }
        }
        writer.finish()?.flush()?;
    }
    mapping.persist(PathBuf::from(mapping_name))?;
    result.persist(output)?;
    timings.finalize = started.elapsed().as_secs_f64();
    stats.timings = timings;
    eprintln!(
        "[unique_alignment_skips]\t{} seed windows skipped across proto boundaries",
        scan_stats.unique_alignment_windows_skipped
    );
    eprintln!(
        "[reference_seed_hits]\t{} dictionary lookups reused",
        scan_stats.reference_seed_hits
    );
    eprintln!(
        "[context_witness_hits]\t{} exact context lookups reused",
        scan_stats.context_witness_hits
    );
    eprintln!(
        "[scan_summary]\t{} source genome reads\t{} contexts\t{} new exact alignments\t{} alignment cache hits\t{} extended bases\t{} first-level cache hits\t{} minimizer-absent windows skipped",
        scan_stats.genomes,
        scan_stats.contexts,
        scan_stats.alignments,
        scan_stats.cache_hits,
        scan_stats.extended_bases,
        scan_stats.hot_cache_hits,
        scan_stats.absent_windows_skipped
    );
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Deliberately slow oracle: revisit every prior interval at each start and
    // collect every covering physical occurrence, without a bounded summary.
    fn exhaustive_reduction(
        contexts: &[Context],
        genomes: usize,
        k: u64,
    ) -> Vec<(u64, Vec<(u32, u64, bool)>)> {
        let mut starts: Vec<_> = contexts.iter().map(|c| c.record.begin).collect();
        starts.sort_unstable();
        starts.dedup();
        let mut previous_end = 0;
        let mut output = Vec::new();
        for begin in starts {
            let ends: Vec<_> = (0..genomes)
                .map(|genome| {
                    contexts
                        .iter()
                        .filter(|c| c.genome == genome && c.record.begin <= begin)
                        .map(|c| c.record.end)
                        .max()
                })
                .collect();
            if ends.iter().any(Option::is_none) {
                continue;
            }
            let end = ends.into_iter().flatten().min().unwrap();
            if end <= previous_end {
                continue;
            }
            previous_end = end;
            if end < begin + k {
                continue;
            }
            let mut hits = Vec::new();
            for genome in 0..genomes {
                let mut copies = BTreeMap::new();
                for c in contexts.iter().filter(|c| {
                    c.genome == genome && c.record.begin <= begin && c.record.end >= end
                }) {
                    let h = interval_hit(c.record, begin, end);
                    copies
                        .entry((h.contig, h.start))
                        .and_modify(|r: &mut bool| *r &= h.reverse)
                        .or_insert(h.reverse);
                }
                if copies.len() != 1 {
                    break;
                }
                let ((contig, start), reverse) = copies.into_iter().next().unwrap();
                hits.push((contig, start, reverse));
            }
            if hits.len() == genomes {
                output.push((end - begin, hits));
            }
        }
        output.sort_unstable();
        output
    }

    #[test]
    fn bounded_occurrences_match_exhaustive_physical_interval_oracle() {
        let mut state = 0x5120_2026_1018u64;
        let mut random = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for genomes in [1, 2, 3, 5, 7] {
            for case in 0..160 {
                let mut contexts = vec![Context {
                    genome: 0,
                    record: Record {
                        begin: 0,
                        end: 80,
                        hit: Hit {
                            pattern: 0,
                            contig: 0,
                            start: 0,
                            reverse: false,
                        },
                    },
                }];
                for genome in 1..genomes {
                    for _ in 0..(4 + case % 48) {
                        let begin = random() % 72;
                        let end = begin + 1 + random() % (80 - begin);
                        let class = random() % (1 + case % 9);
                        let reverse = class % 2 != 0;
                        let offset = 200 + (class / 2) * 100;
                        let record = Record {
                            begin,
                            end,
                            hit: Hit {
                                pattern: 0,
                                contig: (class / 6) as u32,
                                start: offset + if reverse { 80 - end } else { begin },
                                reverse,
                            },
                        };
                        contexts.push(Context { genome, record });
                        // Overlapping records with the same affine mapping and
                        // palindromic pairs must not count as additional copies.
                        if random() % 4 == 0 {
                            contexts.push(Context { genome, record });
                        }
                    }
                }
                let expected = exhaustive_reduction(&contexts, genomes, 3);
                for shards in [1, 3] {
                    let mut actual = Vec::new();
                    for shard in 0..shards {
                        let mut input = contexts.clone();
                        actual.extend(
                            reduce_ranges(&mut input, 3, genomes, shard, shards, 80)
                                .into_iter()
                                .map(|m| {
                                    (
                                        m.len,
                                        m.hits
                                            .into_iter()
                                            .map(|h| (h.contig, h.start, h.reverse))
                                            .collect(),
                                    )
                                }),
                        );
                    }
                    actual.sort_unstable();
                    assert_eq!(
                        actual, expected,
                        "genomes={genomes}, case={case}, shards={shards}"
                    );
                }
            }
        }
    }

    #[test]
    fn occurrence_mapping_identity_handles_wide_coordinates() {
        let a = Record {
            begin: u64::MAX - 80,
            end: u64::MAX - 10,
            hit: Hit {
                pattern: 0,
                contig: 1,
                start: u64::MAX - 100,
                reverse: false,
            },
        };
        let b = Record {
            begin: a.begin + 1,
            end: a.end,
            hit: Hit {
                start: a.hit.start + 1,
                ..a.hit
            },
        };
        assert!(same_mapping(a, b));
        assert!(!same_mapping(
            a,
            Record {
                hit: Hit {
                    start: b.hit.start + 1,
                    ..b.hit
                },
                ..b
            }
        ));
        let a = Record {
            hit: Hit {
                reverse: true,
                ..a.hit
            },
            ..a
        };
        let b = Record {
            begin: a.begin + 1,
            end: a.end - 1,
            hit: Hit {
                start: a.hit.start + 1,
                ..a.hit
            },
        };
        assert!(same_mapping(a, b));
        assert!(!same_mapping(
            a,
            Record {
                hit: Hit { contig: 2, ..b.hit },
                ..b
            }
        ));
    }

    #[test]
    fn parallel_reference_batches_preserve_all_exact_serial_contexts() -> io::Result<()> {
        let temp = tempfile::tempdir()?;
        let mut state = 431u64;
        let mut dna: Vec<_> = (0..9000)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect();
        dna.extend_from_within(1000..3000);
        let options = LocateOptions {
            k: 31,
            m: 15,
            threads: 4,
            chunk_bases: 127,
            forward_only: false,
        };
        let mut source = ProtoIndex {
            options,
            patterns: Vec::new(),
            sequences: PackedSeqVec::default(),
        };
        let mut keys = std::collections::BTreeSet::new();
        for seed in dna.windows(31) {
            keys.insert(seed.to_vec().min(reverse(seed)));
        }
        for seed in keys {
            source.add_pattern(&seed)?;
        }
        let index = FragmentIndex::new(source)?;
        let input = temp.path().join("reference.fa");
        let mut bytes = b">first\n".to_vec();
        bytes.extend(&dna);
        bytes.extend(b"\n>second\n");
        bytes.extend(reverse(&dna[1000..3000]));
        bytes.push(b'\n');
        std::fs::write(&input, bytes)?;
        let reference = read_reference(&index, &input)?;
        assert!(reference.fragments.iter().map(Vec::len).sum::<usize>() > 1024);
        let mut observed = Vec::new();
        for parallel in [false, true] {
            let dir = temp
                .path()
                .join(if parallel { "parallel" } else { "serial" });
            std::fs::create_dir(&dir)?;
            let mut writer = ContextWriter::new(&dir, 0, 4)?;
            let stats = if parallel {
                parallel_reference_contexts(&reference, 4, 31, 4, &mut writer)?
            } else {
                reference_contexts(&reference, 4, 31, &mut writer)?
            };
            assert_eq!(stats.genomes, 1);
            writer.finish()?;
            let mut rows = Vec::new();
            for shard in 0..4 {
                let mut contexts = Vec::new();
                read_worker_contexts(
                    &mut File::open(worker_context_path(&dir, 0, shard))?,
                    &mut contexts,
                    &mut Vec::new(),
                    1,
                )?;
                rows.extend(contexts.into_iter().map(|c| {
                    (
                        shard,
                        c.record.begin,
                        c.record.end,
                        c.record.hit.contig,
                        c.record.hit.start,
                        c.record.hit.reverse,
                    )
                }));
            }
            rows.sort_unstable();
            rows.dedup();
            observed.push(rows);
        }
        assert_eq!(observed[0], observed[1]);
        Ok(())
    }

    #[test]
    fn presence_ranges_exactly_mark_half_open_windows_across_word_boundaries() {
        for start in 0..193 {
            for end in start..193 {
                let mut bits = vec![0; 4];
                mark_position_range(&mut bits, start, end);
                for position in 0..256 {
                    assert_eq!(
                        bits[position / 64] & (1 << (position % 64)) != 0,
                        (start..end).contains(&position),
                        "{start}..{end}: {position}"
                    );
                }
            }
        }
    }

    #[test]
    fn alignment_inline_intervals_match_ordered_map_for_disordered_insertions() {
        let mut actual = AlignmentIntervals::default();
        let mut expected = BTreeMap::new();
        let mut rng = 180019u64;
        for _ in 0..1000 {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            let begin = ((rng >> 32) % 600) as usize;
            let end = begin + (rng as usize % 150);
            actual.insert(begin, end);
            expected.insert(begin, end);
            for position in [0, begin.saturating_sub(1), begin, begin + 20, 500, 700] {
                assert_eq!(
                    actual.preceding(position),
                    expected
                        .range(..=position)
                        .next_back()
                        .map(|(&a, &b)| (a, b))
                );
            }
        }
    }

    #[test]
    fn progressive_scan_order_visits_each_genome_once() {
        for genomes in 0..129 {
            for batch in [0, 1, 2, 7, 256, usize::MAX] {
                let mut ids = genome_scan_order(genomes, batch);
                assert_eq!(ids.len(), genomes.saturating_sub(1));
                ids.sort_unstable();
                assert_eq!(ids, (1..genomes).collect::<Vec<_>>());
            }
        }
    }

    #[test]
    fn pruning_preserves_internal_seed_islands_and_exact_fragment_cache_skips() -> io::Result<()> {
        let mut rng = 180020u64;
        let core: Vec<_> = (0..500)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(rng >> 62) as usize]
            })
            .collect();
        let mut broken = core.clone();
        broken[200] = b'N';
        broken[275] = if broken[275] == b'A' { b'C' } else { b'A' };
        for k in [3, 31, 33, 65] {
            for forward_only in [false, true] {
                for chunk_bases in [1, 7, 128] {
                    let mut source = ProtoIndex {
                        options: LocateOptions {
                            k,
                            m: k.min(21),
                            threads: 3,
                            chunk_bases,
                            forward_only,
                        },
                        patterns: Vec::new(),
                        sequences: PackedSeqVec::default(),
                    };
                    source.add_pattern(&core)?;
                    source.add_pattern(&core[100..300])?;
                    let mut index = FragmentIndex::new(source)?;
                    let before = index.active_positions();
                    let mut keep = vec![0; index.active.len()];
                    for pat in 0..index.patterns.len() {
                        let offset = index.position_offsets[pat];
                        let count = index.position_offsets[pat + 1] - offset;
                        for i in 0..count {
                            if ((30..45).contains(&i) || (100..170).contains(&i)) && i % 7 != 0 {
                                mark_position_range(&mut keep, offset + i, offset + i + 1);
                            }
                        }
                    }
                    assert!(index.retain_present(&keep));
                    assert!(index.active_positions() < before);
                    let heads: Vec<_> = if let Some(dictionary) = &index.posting_dictionary {
                        dictionary
                            .heads
                            .iter()
                            .copied()
                            .filter(|&h| h != NO_POSTING)
                            .collect()
                    } else {
                        index.heads.values().copied().collect()
                    };
                    for head in heads {
                        let mut next = head;
                        while next != NO_POSTING {
                            let p = &index.postings[next as usize];
                            let bit =
                                index.position_offsets[p.pattern as usize] + p.offset as usize;
                            assert_ne!(keep[bit / 64] & (1 << (bit % 64)), 0);
                            next = p.next;
                        }
                    }
                    let queries = [&core, &broken, &reverse(&core)];
                    let actual = queries
                        .iter()
                        .map(|q| fragments(&index, q))
                        .collect::<io::Result<Vec<_>>>()?;
                    index.unique_runs.fill(0);
                    let exhaustive = queries
                        .iter()
                        .map(|q| fragments(&index, q))
                        .collect::<io::Result<Vec<_>>>()?;
                    assert_eq!(
                        actual, exhaustive,
                        "k={k}, forward={forward_only}, chunk={chunk_bases}"
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn extending_pruned_index_restores_raw_candidates_for_union_states() -> io::Result<()> {
        let mut rng = 180021u64;
        let core: Vec<_> = (0..200)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(rng >> 62) as usize]
            })
            .collect();
        let options = LocateOptions {
            k: 31,
            m: 21,
            threads: 2,
            chunk_bases: 7,
            forward_only: false,
        };
        let mut source = ProtoIndex {
            options: options.clone(),
            patterns: Vec::new(),
            sequences: PackedSeqVec::default(),
        };
        source.add_pattern(&core)?;
        let mut index = FragmentIndex::new(source)?;
        assert!(index.retain_present(&vec![0; index.active.len()]));
        assert!(index.heads.is_empty());
        let mut addition = ProtoIndex {
            options,
            patterns: Vec::new(),
            sequences: PackedSeqVec::default(),
        };
        addition.add_pattern(&core[50..150])?;
        index.extend(addition)?;
        assert_eq!(index.active_positions(), (200 - 31 + 1) + (100 - 31 + 1));
        let found = fragments(&index, &core)?;
        assert!(found.iter().any(|f| f.2 == 0 && f.0 == 0 && f.1 == 200));
        assert!(found.iter().any(|f| f.2 == 1 && f.0 == 0 && f.1 == 100));
        Ok(())
    }

    #[test]
    fn vectorizable_dna_validation_preserves_every_ascii_and_ambiguity_boundary() {
        for len in [0, 1, 7, 15, 16, 31, 32, 63, 64, 65, 129] {
            let mut sequence: Vec<_> = b"ACGTacgt".iter().copied().cycle().take(len).collect();
            assert!(all_dna(&sequence));
            for pos in 0..len {
                let previous = sequence[pos];
                for base in 0..=255 {
                    sequence[pos] = base;
                    assert_eq!(
                        all_dna(&sequence),
                        is_dna(base),
                        "len={len} pos={pos} base={base}"
                    );
                }
                sequence[pos] = previous;
            }
        }
    }

    #[test]
    fn worker_context_batches_preserve_genome_ids_dedup_strands_and_empty_shards() -> io::Result<()>
    {
        let scratch = tempfile::tempdir()?;
        let mut writer = ContextWriter::new(scratch.path(), 0, 3)?;
        let record = Record {
            begin: 2,
            end: 39,
            hit: Hit {
                pattern: 0,
                contig: 7,
                start: 53,
                reverse: true,
            },
        };
        for genome in [4, 0, 2] {
            let sink = ContextSink {
                shards: vec![vec![record, record], vec![], vec![]],
                reference_bases: 100,
                k: 3,
            };
            sink.finish_to(&mut writer, genome)?;
        }
        writer.finish()?;
        let bytes = std::fs::read(worker_context_path(scratch.path(), 0, 0))?;
        assert_eq!(bytes.len(), 3 * (16 + RECORD_BYTES));
        let mut decoded = Vec::new();
        read_worker_contexts(&mut bytes.as_slice(), &mut decoded, &mut Vec::new(), 5)?;
        assert_eq!(
            decoded.iter().map(|c| c.genome).collect::<Vec<_>>(),
            [4, 0, 2]
        );
        for c in decoded {
            assert_eq!(
                (c.record.begin, c.record.end, c.record.hit),
                (record.begin, record.end, record.hit)
            );
        }
        for shard in [1, 2] {
            assert_eq!(
                std::fs::metadata(worker_context_path(scratch.path(), 0, shard))?.len(),
                0
            );
        }
        for length in 1..16 + RECORD_BYTES {
            assert!(
                read_worker_contexts(&mut &bytes[..length], &mut Vec::new(), &mut Vec::new(), 5)
                    .is_err()
            );
        }
        assert!(
            read_worker_contexts(&mut bytes.as_slice(), &mut Vec::new(), &mut Vec::new(), 4)
                .is_err()
        );
        let mut writer = ContextWriter::new(scratch.path(), 1, 1)?;
        let records: Vec<_> = (0..4096)
            .map(|i| Record {
                begin: i,
                end: i + 37,
                hit: Hit {
                    start: i + 53,
                    ..record.hit
                },
            })
            .collect();
        let sink = ContextSink {
            shards: vec![records],
            reference_bases: 5000,
            k: 3,
        };
        sink.finish_to(&mut writer, 9)?;
        assert!(std::fs::metadata(worker_context_path(scratch.path(), 1, 0))?.len() > 64 * 1024);
        writer.finish()?;
        let bytes = std::fs::read(worker_context_path(scratch.path(), 1, 0))?;
        let mut decoded = Vec::new();
        read_worker_contexts(&mut bytes.as_slice(), &mut decoded, &mut Vec::new(), 10)?;
        assert_eq!(decoded.len(), 4096);
        assert!(decoded.iter().all(|c| c.genome == 9));
        for (i, c) in decoded.iter().enumerate() {
            assert_eq!(
                (c.record.begin, c.record.hit.start),
                (i as u64, i as u64 + 53)
            );
        }
        Ok(())
    }

    fn fragments(
        index: &FragmentIndex,
        sequence: &[u8],
    ) -> io::Result<Vec<(usize, usize, u32, u64, bool)>> {
        let query = PackedContig::new(sequence, index.options.k);
        let mut result = Vec::new();
        index.scan(&query, 0, 0, |f| {
            result.push((f.begin, f.end, f.hit.pattern, f.hit.start, f.hit.reverse));
            Ok(())
        })?;
        result.sort_unstable();
        Ok(result)
    }

    #[test]
    fn unique_seed_skips_preserve_every_fragment_across_strands_chunks_and_repeats()
    -> io::Result<()> {
        let mut state = 173101u64;
        let core: Vec<_> = (0..1400)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect();
        let mut broken = core.clone();
        broken[367] = b'N';
        broken[896] = b"ACGT"[(broken[896] as usize + 1) % 4];
        let queries = [core.clone(), reverse(&core), broken];
        for k in [3, 7, 15, 31, 32, 33, 63] {
            for forward_only in [false, true] {
                for chunk_bases in [1, 7, 32, 8192] {
                    let mut source = ProtoIndex {
                        options: LocateOptions {
                            k,
                            m: k.min(21),
                            threads: 2,
                            chunk_bases,
                            forward_only,
                        },
                        patterns: Vec::new(),
                        sequences: PackedSeqVec::default(),
                    };
                    source.add_pattern(&core)?;
                    source.add_pattern(&core[113..471])?;
                    source.add_pattern(&reverse(&core[650..900]))?;
                    source.add_pattern(&core[113..471])?;
                    let mut index = FragmentIndex::new(source)?;
                    if k >= 15 {
                        assert!(index.unique_runs.iter().any(|&n| n >= 100));
                    }
                    let actual = queries
                        .iter()
                        .map(|q| fragments(&index, q))
                        .collect::<io::Result<Vec<_>>>()?;
                    index.unique_runs.fill(0);
                    let exhaustive = queries
                        .iter()
                        .map(|q| fragments(&index, q))
                        .collect::<io::Result<Vec<_>>>()?;
                    assert_eq!(
                        actual, exhaustive,
                        "k={k}, forward={forward_only}, chunk={chunk_bases}"
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn adding_proto_invalidates_skips_that_would_hide_its_internal_matches() -> io::Result<()> {
        let mut state = 173102u64;
        let core: Vec<_> = (0..1000)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect();
        let options = LocateOptions {
            k: 31,
            m: 21,
            threads: 2,
            chunk_bases: 8192,
            forward_only: false,
        };
        let mut source = ProtoIndex {
            options: options.clone(),
            patterns: Vec::new(),
            sequences: PackedSeqVec::default(),
        };
        source.add_pattern(&core)?;
        let mut index = FragmentIndex::new(source)?;
        assert!(index.unique_runs.iter().any(|&n| n >= 900));
        let mut added = ProtoIndex {
            options,
            patterns: Vec::new(),
            sequences: PackedSeqVec::default(),
        };
        added.add_pattern(&core[300..700])?;
        index.extend(added)?;
        let actual = fragments(&index, &core)?;
        assert!(actual.iter().any(|f| f.2 == 1 && f.0 == 0 && f.1 == 400));
        index.unique_runs.fill(0);
        assert_eq!(actual, fragments(&index, &core)?);
        Ok(())
    }

    #[test]
    fn fragment_offsets_work_across_packed_word_sizes() -> io::Result<()> {
        let mut state = 5613u64;
        let core: Vec<u8> = (0..320)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect();
        let temp = tempfile::tempdir()?;
        let patterns = temp.path().join("proto.fa.zst");
        let reference = temp.path().join("reference.fa.zst");
        let text = format!(">core\n{}\n", String::from_utf8_lossy(&core));
        let compressed = zstd::stream::encode_all(text.as_bytes(), 1)?;
        std::fs::write(&patterns, &compressed)?;
        std::fs::write(&reference, &compressed)?;
        for k in [15, 31, 32, 33, 63, 64, 67] {
            let ranges = [0..120, 121 - k..240, 241 - k..320];
            let mut query_text = String::new();
            for (id, range) in ranges.iter().enumerate() {
                let dna: Vec<_> = core[range.clone()]
                    .iter()
                    .rev()
                    .copied()
                    .map(complement)
                    .collect();
                query_text.push_str(&format!(">c{id}\n{}\n", String::from_utf8_lossy(&dna)));
            }
            let query = temp.path().join("query.fa.zst");
            std::fs::write(&query, zstd::stream::encode_all(query_text.as_bytes(), 1)?)?;
            let output = temp.path().join("mums.tsv.zst");
            let stats = locate_mums(
                &patterns,
                &[reference.clone(), query],
                &output,
                Some(temp.path()),
                LocateOptions {
                    k,
                    m: k.min(21),
                    threads: 4,
                    chunk_bases: 7,
                    forward_only: false,
                },
            )?;
            let text = String::from_utf8(zstd::stream::decode_all(File::open(output)?)?).unwrap();
            assert_eq!(stats.verified_mums, 3, "k={k}: {text}");
            for (id, range) in ranges.iter().enumerate() {
                assert!(
                    text.contains(&format!(
                        "{}\t{},0\t+,-\t0,{id}\n",
                        range.len(),
                        range.start
                    )),
                    "k={k}: {text}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn incorrect_seed_posting_is_checked_against_exact_dna() -> io::Result<()> {
        let mut source = ProtoIndex {
            options: LocateOptions {
                k: 3,
                m: 2,
                threads: 2,
                chunk_bases: 2,
                forward_only: false,
            },
            patterns: Vec::new(),
            sequences: PackedSeqVec::default(),
        };
        source.add_pattern(b"ACGATC")?;
        source.add_pattern(b"TTTAAA")?;
        let mut index = FragmentIndex::new(source)?;
        let sequence = PackedContig::new(b"ACGATC", 3);
        let dictionary = index.posting_dictionary.as_mut().unwrap();
        let key = sshash_lib::Kmer::<3>::from_ascii_unchecked(b"ACG");
        let found = dictionary.compact.dictionary.query::<3>(&key);
        assert!(found.is_found());
        let slot = found.kmer_id as usize * 2 + usize::from(found.kmer_orientation < 0);
        let next = dictionary.heads[slot];
        let posting = index.postings.len() as u32;
        index.postings.push(SeedPosting {
            pattern: 1,
            offset: 0,
            next,
            reverse: false,
        });
        dictionary.heads[slot] = posting;
        let mut hits = Vec::new();
        index.scan(&sequence, 0, 0, |fragment| {
            hits.push(fragment);
            Ok(())
        })?;
        assert!(
            hits.iter()
                .any(|f| f.hit.pattern == 0 && f.begin == 0 && f.end == 6)
        );
        assert!(hits.iter().all(|f| f.hit.pattern != 1));
        Ok(())
    }

    fn reverse(sequence: &[u8]) -> Vec<u8> {
        sequence.iter().rev().copied().map(complement).collect()
    }

    #[test]
    fn exact_cached_contexts_match_basewise_oracle() -> io::Result<()> {
        let mut state = 91u64;
        let mut sequence: Vec<u8> = (0..1400)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect();
        for pos in [0, 97, 431, 997, 1399] {
            sequence[pos] = b'N';
        }
        for reverse_query in [false, true] {
            let mut query = b"AACCA".to_vec();
            query.extend_from_slice(&sequence);
            query.extend_from_slice(b"TGT");
            for pos in [34 + 5, 451 + 5, 1103 + 5] {
                query[pos] = if query[pos] == b'A' { b'C' } else { b'A' };
            }
            let physical = if reverse_query {
                reverse(&query)
            } else {
                query.clone()
            };
            let packed = PackedContig::new(&physical, 7);
            let reference = Reference {
                contigs: vec![PackedContig::new(&sequence, 7)],
                meta: GenomeMeta {
                    input: "reference".into(),
                    contigs: vec![ContigMeta {
                        name: "chr".into(),
                        start: 37,
                        len: sequence.len(),
                    }],
                },
                fragments: Vec::new(),
            };
            for descending in [false, true] {
                let mut cache = AlignmentCache::default();
                let mut monotonic_cache = MonotonicAlignmentCache::default();
                let mut monotonic_stats = ScanStats::default();
                let mut stats = ScanStats::default();
                let mut checked = 0;
                let mut positions: Vec<_> = (0..sequence.len() - 7).step_by(3).collect();
                if descending {
                    positions.reverse();
                }
                for r in positions {
                    if !sequence[r..r + 7].iter().copied().all(is_dna)
                        || sequence[r..r + 7] != query[r + 5..r + 12]
                    {
                        continue;
                    }
                    let mut left = 0usize;
                    while r > left
                        && is_dna(sequence[r - left - 1])
                        && sequence[r - left - 1] == query[r + 5 - left - 1]
                    {
                        left += 1;
                    }
                    let mut right = 0usize;
                    while r + 7 + right < sequence.len()
                        && is_dna(sequence[r + 7 + right])
                        && sequence[r + 7 + right] == query[r + 12 + right]
                    {
                        right += 1;
                    }
                    for reverse_anchor in [false, true] {
                        let anchor = Anchor {
                            hit: Hit {
                                pattern: 0,
                                contig: 0,
                                start: 37 + r as u64,
                                reverse: reverse_anchor,
                            },
                            len: 7,
                        };
                        let q = if reverse_query {
                            physical.len() - (r + 5) - 7
                        } else {
                            r + 5
                        };
                        let hit = Hit {
                            pattern: 0,
                            contig: 0,
                            start: 71 + q as u64,
                            reverse: reverse_query != reverse_anchor,
                        };
                        let expected = if reverse_anchor {
                            (right as u64, left as u64)
                        } else {
                            (left as u64, right as u64)
                        };
                        assert_eq!(
                            alignment_context(
                                &reference, &anchor, &packed, hit, 71, &mut cache, &mut stats
                            )?,
                            expected
                        );
                        assert_eq!(
                            alignment_context_monotonic(
                                &reference,
                                &anchor,
                                &packed,
                                hit,
                                71,
                                &mut monotonic_cache,
                                &mut monotonic_stats,
                            )?,
                            expected
                        );
                        assert_eq!(
                            alignment_context_uncached(
                                &reference,
                                &anchor,
                                &packed,
                                hit,
                                71,
                                &mut ScanStats::default()
                            )?,
                            expected
                        );
                        checked += 1;
                    }
                }
                assert!(checked > 500);
                assert!(
                    stats.alignments < 20,
                    "{} alignments for {} seeds",
                    stats.alignments,
                    checked
                );
                assert!(stats.cache_hits > 500);
                assert_eq!(monotonic_stats.alignments, stats.alignments);
                assert_eq!(monotonic_stats.cache_hits, stats.cache_hits);
            }
        }
        Ok(())
    }

    #[test]
    fn long_identical_region_is_extended_once_for_many_seeds() -> io::Result<()> {
        let sequence: Vec<_> = b"ACGATCTG".iter().copied().cycle().take(500_003).collect();
        let query = PackedContig::new(&sequence, 31);
        let reference = Reference {
            contigs: vec![PackedContig::new(&sequence, 31)],
            meta: GenomeMeta {
                input: "reference".into(),
                contigs: vec![ContigMeta {
                    name: "chr".into(),
                    start: 0,
                    len: sequence.len(),
                }],
            },
            fragments: Vec::new(),
        };
        let mut stats = ScanStats::default();
        let mut cache = AlignmentCache::default();
        for pos in (13..sequence.len() - 31).step_by(1237) {
            let hit = Hit {
                pattern: 0,
                contig: 0,
                start: pos as u64,
                reverse: false,
            };
            let anchor = Anchor { hit, len: 31 };
            assert_eq!(
                alignment_context(&reference, &anchor, &query, hit, 0, &mut cache, &mut stats)?,
                (pos as u64, (sequence.len() - pos - 31) as u64)
            );
        }
        assert_eq!(stats.alignments, 1);
        assert!(stats.cache_hits > 400);
        assert_eq!(stats.extended_bases, (sequence.len() - 31) as u64);
        Ok(())
    }
}
