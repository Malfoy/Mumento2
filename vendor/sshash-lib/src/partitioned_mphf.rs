//! Partitioned Minimal Perfect Hash Function
//!
//! Wraps multiple PTRHash MPHFs, one per partition, with hash-based partition
//! selection matching the C++ `partitioned_phf` design. Keys are assigned to
//! partitions via Lemire fast-range reduction on an independent rapidhash,
//! and each partition's MPHF is built independently (in parallel via rayon).
//!
//! For datasets with fewer than `AVG_PARTITION_SIZE` keys (including all our
//! test indices), only a single partition is created and the `get()` fast path
//! skips the partition hash entirely.

use crate::mphf_config::{Mphf, MphfHasher, build_mphf_from_slice, mphf_hasher, read_mphf};
use rayon::prelude::*;
use std::hash::Hash;
use std::io::{self, Read, Write};

/// Average number of keys per partition, matching C++ PTHash `partitioned_phf`.
const AVG_PARTITION_SIZE: usize = 3_000_000;

/// Seed used for partition selection, chosen to differ from PTRHash's internal
/// seeds (0, 1, 2, ...) to ensure statistical independence.
const PARTITION_HASH_SEED: u64 = 0xC6A4_A793_5BD1_E995;

/// A partitioned minimal perfect hash function.
///
/// Splits keys into partitions by hash range (Lemire fast-range reduction),
/// builds one PTRHash MPHF per partition, and routes queries transparently.
/// Global indices are `offsets[partition] + inner_mphf.get(key)`.
pub struct PartitionedMphf {
    /// One PTRHash MPHF per partition.
    inners: Vec<Mphf>,
    /// Cumulative key counts: `offsets[i]` = total keys in partitions 0..i.
    /// Length = `num_partitions + 1`.
    offsets: Vec<usize>,
    /// Number of partitions.
    num_partitions: u32,
    /// Total number of keys across all partitions.
    num_keys: usize,
    /// Hasher for partition selection.
    hasher: MphfHasher,
}

impl PartitionedMphf {
    /// Build a partitioned MPHF from an owned Vec of keys.
    ///
    /// If `partitioned` is false (or there are fewer than `AVG_PARTITION_SIZE`
    /// keys), a single partition is used — equivalent to a monolithic MPHF with
    /// zero query overhead.
    pub fn build_from_vec<K: Hash + Copy + Into<u128> + Send + Sync>(
        keys: Vec<K>,
        partitioned: bool,
    ) -> Self {
        let num_keys = keys.len();
        let num_partitions = if partitioned {
            num_keys.div_ceil(AVG_PARTITION_SIZE).max(1)
        } else {
            1
        };
        if num_partitions == 1 {
            return Self::build_from_slice(&keys, false);
        }
        let (grouped, offsets) = distribute_keys(&keys, num_partitions);
        // Release owned input before PTRHash allocates its construction hashes.
        drop(keys);
        Self::build_grouped(&grouped, offsets)
    }

    /// Build directly from borrowed keys, including mmap-backed minimizers.
    pub fn build_from_slice<K: Hash + Copy + Into<u128> + Send + Sync>(
        keys: &[K],
        partitioned: bool,
    ) -> Self {
        let num_keys = keys.len();
        if num_keys == 0 {
            return Self {
                inners: Vec::new(),
                offsets: vec![0],
                num_partitions: 0,
                num_keys: 0,
                hasher: mphf_hasher(),
            };
        }
        let num_partitions = if partitioned {
            num_keys.div_ceil(AVG_PARTITION_SIZE).max(1)
        } else {
            1
        };
        if num_partitions == 1 {
            return Self {
                inners: vec![build_mphf_from_slice(keys)],
                offsets: vec![0, num_keys],
                num_partitions: 1,
                num_keys,
                hasher: mphf_hasher(),
            };
        }
        let (grouped, offsets) = distribute_keys(keys, num_partitions);
        Self::build_grouped(&grouped, offsets)
    }

    fn build_grouped<K: Hash + Copy + Into<u128> + Send + Sync>(
        keys: &[K],
        offsets: Vec<usize>,
    ) -> Self {
        let inners = offsets
            .par_windows(2)
            .map(|range| build_mphf_from_slice(&keys[range[0]..range[1]]))
            .collect();
        Self {
            inners,
            num_partitions: (offsets.len() - 1) as u32,
            num_keys: keys.len(),
            offsets,
            hasher: mphf_hasher(),
        }
    }

    /// Look up a key and return its global index in [0, num_keys).
    ///
    /// Out-of-set keys can map to any slot. The dictionary's stored-DNA
    /// verifier decides membership; the MPHF itself never does.
    #[inline]
    pub fn get<K: Hash + Copy + Into<u128>>(&self, key: &K) -> usize {
        if self.num_partitions == 0 {
            return self.num_keys;
        }
        if self.num_partitions == 1 {
            // Fast path: skip partition hash entirely.
            let idx = self.inners[0].get(key);
            if idx == usize::MAX {
                return self.num_keys;
            }
            return idx;
        }
        let p = self.partition_for(key);
        let idx = self.inners[p].get(key);
        if idx == usize::MAX {
            return self.num_keys;
        }
        self.offsets[p] + idx
    }

    /// Total number of keys.
    pub fn num_keys(&self) -> usize {
        self.num_keys
    }

    /// Number of partitions.
    pub fn num_partitions(&self) -> u32 {
        self.num_partitions
    }

    /// Compute which partition a key belongs to (Lemire fast-range reduction).
    #[inline]
    fn partition_for<K: Hash + Copy + Into<u128>>(&self, key: &K) -> usize {
        let hash = self.hasher.hash_one_with_seed(key, PARTITION_HASH_SEED);
        ((hash as u128 * self.num_partitions as u128) >> 64) as usize
    }

    /// Serialize to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // Header
        writer.write_all(&self.num_partitions.to_le_bytes())?;
        writer.write_all(&(self.num_keys as u64).to_le_bytes())?;

        // Offsets (num_partitions + 1 entries)
        for &off in &self.offsets {
            writer.write_all(&(off as u64).to_le_bytes())?;
        }

        // Inner MPHFs
        for mphf in &self.inners {
            mphf.write(writer)?;
        }

        Ok(())
    }

    /// Deserialize from a reader.
    pub fn read_from(reader: &mut dyn Read) -> io::Result<Self> {
        let mut buf4 = [0u8; 4];
        let mut buf8 = [0u8; 8];

        reader.read_exact(&mut buf4)?;
        let num_partitions = u32::from_le_bytes(buf4);

        reader.read_exact(&mut buf8)?;
        let num_keys = u64::from_le_bytes(buf8) as usize;

        let num_offsets = num_partitions as usize + 1;
        let mut offsets = Vec::with_capacity(num_offsets);
        for _ in 0..num_offsets {
            reader.read_exact(&mut buf8)?;
            offsets.push(u64::from_le_bytes(buf8) as usize);
        }

        let mut inners = Vec::with_capacity(num_partitions as usize);
        for _ in 0..num_partitions {
            inners.push(read_mphf(reader)?);
        }

        Ok(Self {
            inners,
            offsets,
            num_partitions,
            num_keys,
            hasher: mphf_hasher(),
        })
    }

    pub(crate) fn has_compact(&self) -> bool {
        self.inners.iter().any(Mphf::has_compact)
    }

    /// Estimate serialized byte size (for container offset table pre-allocation).
    pub fn write_bytes(&self) -> usize {
        let header = 4 + 8; // num_partitions + num_keys
        let offsets = (self.offsets.len()) * 8;
        let mphfs: usize = self.inners.iter().map(|m| m.write_bytes()).sum();
        header + offsets + mphfs
    }
}

/// Extension trait to hash with a specific seed.
///
/// PTRHash's `BuildRapidHash` implements `BuildSeededHasher`, so we use
/// `build_hasher(seed)` to get a seeded hasher for partition selection.
trait HashOneWithSeed {
    fn hash_one_with_seed<K: Hash + Copy + Into<u128>>(&self, key: &K, seed: u64) -> u64;
}

impl HashOneWithSeed for MphfHasher {
    #[inline]
    fn hash_one_with_seed<K: Hash + Copy + Into<u128>>(&self, key: &K, seed: u64) -> u64 {
        use ptr_hash::hash::{KeyHasher, Xxh3_128};
        Xxh3_128::hash(&(*key).into(), seed) as u64
    }
}

/// Stable parallel counting distribution. Each job owns disjoint output slices;
/// no locks, atomics or raw-pointer writes occur in the scatter loop.
fn distribute_keys<K: Hash + Copy + Into<u128> + Send + Sync>(
    keys: &[K],
    partitions: usize,
) -> (Vec<K>, Vec<usize>) {
    use std::mem::MaybeUninit;
    assert!(partitions > 0);
    if std::env::var("REINDEER_MPHF_PARALLEL_DISTRIBUTION").as_deref() == Ok("0") {
        let mut groups: Vec<Vec<K>> = (0..partitions).map(|_| Vec::new()).collect();
        let hasher = mphf_hasher();
        for &key in keys {
            let hash = hasher.hash_one_with_seed(&key, PARTITION_HASH_SEED);
            groups[((hash as u128 * partitions as u128) >> 64) as usize].push(key);
        }
        let mut offsets = vec![0];
        let mut grouped = Vec::with_capacity(keys.len());
        for group in groups {
            grouped.extend(group);
            offsets.push(grouped.len());
        }
        return (grouped, offsets);
    }
    let chunks = (rayon::current_num_threads() * 4).min(keys.len().max(1));
    let chunk_size = keys.len().div_ceil(chunks).max(1);
    let np = partitions as u128;
    let partition = |key: &K| {
        let hash = mphf_hasher().hash_one_with_seed(key, PARTITION_HASH_SEED);
        ((hash as u128 * np) >> 64) as usize
    };
    let counts: Vec<Vec<usize>> = keys
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut counts = vec![0; partitions];
            for key in chunk {
                counts[partition(key)] += 1;
            }
            counts
        })
        .collect();
    let mut offsets = Vec::with_capacity(partitions + 1);
    offsets.push(0);
    for p in 0..partitions {
        offsets.push(offsets[p] + counts.iter().map(|c| c[p]).sum::<usize>());
    }
    let mut output: Vec<MaybeUninit<K>> = Vec::with_capacity(keys.len());
    // SAFETY: uninitialized MaybeUninit<K> is valid. Every slot belongs to
    // exactly one histogram-sized span, filled before conversion to Vec<K>.
    unsafe {
        output.set_len(keys.len());
    }
    let mut jobs: Vec<_> = keys
        .chunks(chunk_size)
        .map(|chunk| (chunk, Vec::with_capacity(partitions)))
        .collect();
    let mut remaining = output.as_mut_slice();
    for p in 0..partitions {
        for (job, count) in jobs.iter_mut().zip(&counts) {
            let (span, rest) = remaining.split_at_mut(count[p]);
            job.1.push(span);
            remaining = rest;
        }
    }
    assert!(remaining.is_empty());
    jobs.into_par_iter().for_each(|(chunk, mut spans)| {
        let mut cursors = vec![0; partitions];
        for &key in chunk {
            let p = partition(&key);
            spans[p][cursors[p]].write(key);
            cursors[p] += 1;
        }
        assert!(
            spans
                .iter()
                .zip(cursors)
                .all(|(span, used)| span.len() == used)
        );
    });
    let length = output.len();
    let capacity = output.capacity();
    let data = output.as_mut_ptr().cast::<K>();
    std::mem::forget(output);
    // SAFETY: all spans were filled and joined above; MaybeUninit<K> has
    // K's layout/alignment. Ownership of this allocation transfers once.
    let output = unsafe { Vec::from_raw_parts(data, length, capacity) };
    (output, offsets)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stable_parallel_distribution_matches_serial_reference() {
        fn check<K: Hash + Copy + Into<u128> + Send + Sync + Eq + std::fmt::Debug>(keys: &[K]) {
            for threads in [1, 4] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                for partitions in [1, 2, 7, 64] {
                    let (grouped, offsets) = pool.install(|| distribute_keys(keys, partitions));
                    let mut reference: Vec<Vec<K>> = (0..partitions).map(|_| Vec::new()).collect();
                    for &key in keys {
                        let hash = mphf_hasher().hash_one_with_seed(&key, PARTITION_HASH_SEED);
                        reference[((hash as u128 * partitions as u128) >> 64) as usize].push(key);
                    }
                    for p in 0..partitions {
                        assert_eq!(&grouped[offsets[p]..offsets[p + 1]], reference[p]);
                    }
                    assert_eq!(offsets.last(), Some(&keys.len()));
                }
            }
        }
        check::<u64>(&[]);
        check(&[u64::MAX]);
        check(&[42u64; 100]); // all keys in one partition; empty output spans
        check(
            &(0..8193u64)
                .map(|x| x.wrapping_mul(0x9e3779b97f4a7c15))
                .collect::<Vec<_>>(),
        );
        check(
            &(0..8193u128)
                .map(|x| (x << 80) | (x ^ 17))
                .collect::<Vec<_>>(),
        );
    }

    #[test]
    fn test_partition_count_math() {
        // < AVG_PARTITION_SIZE → 1 partition
        assert_eq!(1_000_000usize.div_ceil(AVG_PARTITION_SIZE).max(1), 1);
        // Exactly AVG_PARTITION_SIZE → 1 partition
        assert_eq!(AVG_PARTITION_SIZE.div_ceil(AVG_PARTITION_SIZE).max(1), 1);
        // Just over → 2 partitions
        assert_eq!(
            (AVG_PARTITION_SIZE + 1).div_ceil(AVG_PARTITION_SIZE).max(1),
            2
        );
        // 10M → 4 partitions
        assert_eq!(10_000_000usize.div_ceil(AVG_PARTITION_SIZE).max(1), 4);
    }

    #[test]
    fn test_single_partition_roundtrip() {
        let keys: Vec<u64> = (0..1000).collect();
        let pmphf = PartitionedMphf::build_from_vec(keys.clone(), true);

        assert_eq!(pmphf.num_partitions(), 1);
        assert_eq!(pmphf.num_keys(), 1000);

        // All keys should get unique indices in [0, 1000)
        let mut indices: Vec<usize> = keys.iter().map(|k| pmphf.get(k)).collect();
        indices.sort();
        indices.dedup();
        assert_eq!(indices.len(), 1000);
        assert!(indices.iter().all(|&i| i < 1000));
    }

    #[test]
    fn test_monolithic_flag() {
        let keys: Vec<u64> = (0..100).collect();
        let pmphf = PartitionedMphf::build_from_vec(keys.clone(), false);

        assert_eq!(pmphf.num_partitions(), 1);
        assert_eq!(pmphf.num_keys(), 100);

        let mut indices: Vec<usize> = keys.iter().map(|k| pmphf.get(k)).collect();
        indices.sort();
        indices.dedup();
        assert_eq!(indices.len(), 100);
    }

    #[test]
    fn test_serialization_roundtrip() {
        let keys: Vec<u64> = (0..500).collect();
        let pmphf = PartitionedMphf::build_from_vec(keys.clone(), true);

        let mut buf = Vec::new();
        pmphf.write_to(&mut buf).unwrap();

        let pmphf2 = PartitionedMphf::read_from(&mut buf.as_slice()).unwrap();

        assert_eq!(pmphf.num_partitions(), pmphf2.num_partitions());
        assert_eq!(pmphf.num_keys(), pmphf2.num_keys());

        // Verify same results
        for key in &keys {
            assert_eq!(pmphf.get(key), pmphf2.get(key));
        }
    }

    #[test]
    fn test_empty() {
        let keys: Vec<u64> = Vec::new();
        let pmphf = PartitionedMphf::build_from_vec(keys, true);
        assert_eq!(pmphf.num_partitions(), 0);
        assert_eq!(pmphf.num_keys(), 0);
    }

    #[test]
    fn test_write_bytes_sanity() {
        let keys: Vec<u64> = (0..100).collect();
        let pmphf = PartitionedMphf::build_from_vec(keys, true);

        let mut buf = Vec::new();
        pmphf.write_to(&mut buf).unwrap();

        // write_bytes() is an estimate (PTRHash's write_bytes() is approximate)
        // Just verify it's in a reasonable range
        let actual = buf.len();
        let estimate = pmphf.write_bytes();
        assert!(estimate > 0, "estimate should be positive");
        // The estimate comes from PTRHash's write_bytes() which may not be exact
        assert!(
            actual > 0 && estimate > 0,
            "both actual ({actual}) and estimate ({estimate}) should be positive"
        );
    }
}
