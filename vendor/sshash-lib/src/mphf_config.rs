//! PTRHash backend for SSHash's minimizer control map and skew partitions.
//! Both integer widths map injectively to u128; SSHash verifies stored DNA.
use ptr_epserde::prelude::{Deserialize, Serialize};
use ptr_hash::{
    PtrHash, PtrHashParams,
    bucket_fn::{CubicEps, Linear},
    hash::Xxh3_128,
};
use rayon::prelude::*;
use std::{
    hash::{BuildHasherDefault, Hash},
    io::{self, Read, Write},
};
type Inner = PtrHash<u128, Linear, Vec<u32>, Xxh3_128, Vec<u8>, false>;
type CompactInner = PtrHash<
    u128,
    CubicEps,
    cacheline_ef::CachelineEfVec<Vec<cacheline_ef::CachelineEf>>,
    Xxh3_128,
    Vec<u8>,
    false,
>;
const COMPACT_TAG: u64 = 1 << 63;
enum Backend {
    Fast(Inner),
    Compact(CompactInner),
}
impl Backend {
    fn pilot_bits_per_key(&self) -> f64 {
        match self {
            Self::Fast(mph) => mph.bits_per_element().0,
            Self::Compact(mph) => mph.bits_per_element().0,
        }
    }
    #[inline]
    fn index(&self, key: &u128) -> usize {
        match self {
            Self::Fast(mph) => mph.index(key),
            Self::Compact(mph) => mph.index(key),
        }
    }
    fn write(&self, writer: &mut impl Write) -> io::Result<()> {
        match self {
            Self::Fast(mph) => unsafe { mph.serialize(writer) }
                .map(|_| ())
                .map_err(io::Error::other),
            Self::Compact(mph) => unsafe { mph.serialize(writer) }
                .map(|_| ())
                .map_err(io::Error::other),
        }
    }
}
/// Deterministic partition hasher, independent of CPU instruction selection.
pub type MphfHasher = BuildHasherDefault<xxhash_rust::xxh3::Xxh3Default>;
/// Internal SSHash MPHF: PTRHash, also used by skew partitions.
pub struct Mphf {
    inner: Option<Backend>,
    serialized_bytes: usize,
}
/// Construct the partition hasher.
pub fn mphf_hasher() -> MphfHasher {
    MphfHasher::default()
}
struct ByteCounter(usize);
impl Write for ByteCounter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 += bytes.len();
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
impl Mphf {
    /// Return an index; SSHash verifies membership against stored DNA.
    #[inline]
    pub fn get<T: Copy + Into<u128>>(&self, key: &T) -> usize {
        self.inner
            .as_ref()
            .map_or(usize::MAX, |mph| mph.index(&(*key).into()))
    }
    pub(crate) fn has_compact(&self) -> bool {
        matches!(self.inner, Some(Backend::Compact(_)))
    }
    /// Size of the length-delimited serialized PTRHash.
    pub fn write_bytes(&self) -> usize {
        8 + self.serialized_bytes
    }
    /// Write PTRHash for SSHash's MPHF container.
    pub fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        let length = self.serialized_bytes as u64;
        assert!(length < COMPACT_TAG);
        writer.write_all(
            &(length | if self.has_compact() { COMPACT_TAG } else { 0 }).to_le_bytes(),
        )?;
        if let Some(inner) = &self.inner {
            inner.write(writer)?;
        }
        Ok(())
    }
}
/// Build minimizer or canonical k-mer PTRHash functions.
pub fn build_mphf_from_slice<T: Hash + Copy + Into<u128> + Send + Sync>(keys: &[T]) -> Mphf {
    if keys.is_empty() {
        return Mphf {
            inner: None,
            serialized_bytes: 0,
        };
    }
    // The iterator preserves input order and produces the same u128 keys as
    // the former Vec. PTRHash hashes them directly, without a second key array.
    let compact =
        std::env::var("REINDEER_MPHF_COMPACT").as_deref() == Ok("1") && keys.len() >= 4096;
    let key_copy = std::env::var("REINDEER_MPHF_KEY_COPY").as_deref() == Ok("1");
    let inner = if compact {
        if key_copy {
            let words: Vec<u128> = keys.iter().map(|&key| key.into()).collect();
            Backend::Compact(CompactInner::new(&words, PtrHashParams::default_compact()))
        } else {
            Backend::Compact(CompactInner::new_from_par_iter(
                keys.len(),
                keys.par_iter().map(|&key| Into::<u128>::into(key)),
                PtrHashParams::default_compact(),
            ))
        }
    } else if key_copy {
        let words: Vec<u128> = keys.iter().map(|&key| key.into()).collect();
        Backend::Fast(Inner::new(&words, PtrHashParams::default()))
    } else {
        Backend::Fast(Inner::new_from_par_iter(
            keys.len(),
            keys.par_iter().map(|&key| Into::<u128>::into(key)),
            PtrHashParams::default(),
        ))
    };
    let mut count = ByteCounter(0);
    // PTRHash 2.1.2's iterator constructor discards compute_pilots' Option.
    // A failed construction leaves pilots empty; never return that dictionary.
    assert!(
        inner.pilot_bits_per_key() > 0.0,
        "PTRHash iterator construction failed after its ten seed attempts"
    );
    inner
        .write(&mut count)
        .expect("Count PTRHash serialization bytes");
    Mphf {
        inner: Some(inner),
        serialized_bytes: count.0,
    }
}
/// Build from an owned vector.
pub fn build_mphf_from_vec<T: Hash + Copy + Into<u128> + Send + Sync>(keys: Vec<T>) -> Mphf {
    build_mphf_from_slice(&keys)
}
/// Build parts in parallel using the requested Rayon pool.
pub fn build_mphf_from_slice_mt<T: Hash + Copy + Into<u128> + Send + Sync>(
    keys: &[T],
    threads: usize,
) -> Mphf {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("PTRHash thread pool");
    pool.install(|| build_mphf_from_slice(keys))
}
/// Read the PTRHash backend used in this SSHash fork.
pub fn read_mphf(reader: &mut dyn Read) -> io::Result<Mphf> {
    let mut bytes = [0; 8];
    reader.read_exact(&mut bytes)?;
    let tagged_length = u64::from_le_bytes(bytes);
    let compact = tagged_length & COMPACT_TAG != 0;
    let len = usize::try_from(tagged_length & !COMPACT_TAG).map_err(io::Error::other)?;
    if len == 0 {
        return Ok(Mphf {
            inner: None,
            serialized_bytes: 0,
        });
    }
    let mut blob = vec![0; len];
    reader.read_exact(&mut blob)?;
    // Dictionary::load checked the fork's format header; the types and
    // epserde version match the data emitted by write().
    let inner = if compact {
        Backend::Compact(
            unsafe { CompactInner::deserialize_full(&mut io::Cursor::new(blob)) }
                .map_err(io::Error::other)?,
        )
    } else {
        Backend::Fast(
            unsafe { Inner::deserialize_full(&mut io::Cursor::new(blob)) }
                .map_err(io::Error::other)?,
        )
    };
    Ok(Mphf {
        inner: Some(inner),
        serialized_bytes: len,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore = "isolated before/after construction benchmark"]
    fn paper_mphf_construction_trial() {
        let n = std::env::var("REINDEER_TRIAL_KEYS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(10_000_000usize);
        let threads = std::env::var("REINDEER_TRIAL_THREADS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(64usize);
        let keys: Vec<u128> = (0..n as u128)
            .map(|i| (i << 80) | i.wrapping_mul(0x9e3779b97f4a7c15))
            .collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let started = std::time::Instant::now();
        let mphf = pool
            .install(|| crate::partitioned_mphf::PartitionedMphf::build_from_slice(&keys, true));
        let seconds = started.elapsed().as_secs_f64();
        let mut used = vec![false; n];
        for key in &keys {
            let slot = mphf.get(key);
            assert!(slot < n && !used[slot]);
            used[slot] = true;
        }
        assert!(used.iter().all(|&occupied| occupied));
        eprintln!(
            "PAPER_TRIAL {{\"kind\":\"mphf\",\"keys\":{n},\"threads\":{threads},\"construction_seconds\":{seconds},\"serialized_bytes\":{},\"verified\":true}}",
            mphf.write_bytes()
        );
    }

    #[test]
    fn iterator_builder_preserves_permutation_and_serialization() {
        let keys: Vec<u128> = (0..20000u128)
            .map(|i| (i << 80) | i.wrapping_mul(0x9e3779b97f4a7c15))
            .collect();
        let reference = Inner::new(&keys, PtrHashParams::default());
        let candidate = Inner::new_from_par_iter(
            keys.len(),
            keys.par_iter().copied(),
            PtrHashParams::default(),
        );
        // Parallel PTRHash may choose different valid pilots on equal inputs.
        // Each mapping must be a bijection and survive its own serialization.
        for mph in [reference, candidate] {
            let mut slots: Vec<_> = keys.iter().map(|key| mph.index(key)).collect();
            slots.sort_unstable();
            assert_eq!(slots, (0..keys.len()).collect::<Vec<_>>());
            let mut bytes = Vec::new();
            unsafe { mph.serialize(&mut bytes) }.unwrap();
            let loaded = unsafe { Inner::deserialize_full(&mut bytes.as_slice()) }.unwrap();
            for key in &keys {
                assert_eq!(mph.index(key), loaded.index(key));
            }
        }
    }

    #[test]
    fn compact_backend_is_a_permutation_and_roundtrips_with_tag() {
        let keys: Vec<u128> = (0..20000u128).map(|i| (i << 80) | (i ^ 97)).collect();
        let backend = Backend::Compact(CompactInner::new_from_par_iter(
            keys.len(),
            keys.par_iter().copied(),
            PtrHashParams::default_compact(),
        ));
        let mut count = ByteCounter(0);
        backend.write(&mut count).unwrap();
        let mphf = Mphf {
            inner: Some(backend),
            serialized_bytes: count.0,
        };
        let mut slots: Vec<_> = keys.iter().map(|key| mphf.get(key)).collect();
        slots.sort_unstable();
        assert_eq!(slots, (0..keys.len()).collect::<Vec<_>>());
        let mut blob = Vec::new();
        mphf.write(&mut blob).unwrap();
        assert_ne!(
            u64::from_le_bytes(blob[..8].try_into().unwrap()) & COMPACT_TAG,
            0
        );
        let loaded = read_mphf(&mut blob.as_slice()).unwrap();
        assert!(loaded.has_compact());
        for key in &keys {
            assert_eq!(mphf.get(key), loaded.get(key));
        }
    }
}
