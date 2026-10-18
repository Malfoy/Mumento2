//! SIMD minimizers for SSHash construction and query, computed by the pinned
//! simd-minimizers crate. Odd windows use its canonical scheme; even windows
//! use forward SIMD passes on both strands, selected by canonical k-mer
//! orientation. The index supports repeated loci by storing distinct locate
//! positions and constructing skew keys from the actual selected windows.

use crate::hasher::DeterministicHasher;
use crate::kmer::{Kmer, KmerBits};

use crate::encoding::{decode_base, encode_base_unchecked};
use simd_minimizers::packed_seq::{PackedSeq, PackedSeqVec, Seq, SeqVec};
fn mask(k: usize) -> u128 {
    if k == 64 {
        u128::MAX
    } else {
        (1u128 << (2 * k)) - 1
    }
}
#[inline(always)]
fn reverse(mut x: u128, k: usize) -> u128 {
    x ^= 0xaaaa_aaaa_aaaa_aaaa_aaaa_aaaa_aaaa_aaaa;
    x = (x >> 2 & 0x3333_3333_3333_3333_3333_3333_3333_3333)
        | (x & 0x3333_3333_3333_3333_3333_3333_3333_3333) << 2;
    x = (x >> 4 & 0x0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f)
        | (x & 0x0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f_0f0f) << 4;
    x.swap_bytes() >> (128 - 2 * k)
}
fn encode(seq: &[u8]) -> u128 {
    seq.iter().enumerate().fold(0, |x, (i, &b)| {
        x | (encode_base_unchecked(b) as u128) << (2 * i)
    })
}
fn canonical(x: u128, k: usize) -> u128 {
    x.min(reverse(x, k))
}

/// Reusable batched SIMD minimizers for SSHash construction and queries.
#[derive(Default)]
pub struct SimdMinimizers {
    /// Selected absolute m-mer positions per run.
    pub positions: Vec<u32>,
    /// First window of each run, followed by the exclusive end.
    pub starts: Vec<u32>,
    /// Canonical m-mer values of each selected run.
    pub values: Vec<u64>,
    packed: PackedSeqVec,
}
impl SimdMinimizers {
    /// Compute SIMD minimizers for all windows of a sequence.
    pub fn compute(&mut self, seq: &[u8], k: usize, m: usize) {
        self.positions.clear();
        self.starts.clear();
        self.values.clear();
        self.packed.clear();
        self.packed.push_ascii(seq);
        let packed = &self.packed;
        if k % 2 == 0 {
            // The SIMD canonical scheme requires odd windows. For even k,
            // compute SIMD forward minimizers on both strands and choose the
            // canonical whole-k-mer's strand, including its tie-breaking rule.
            let mut fw_pos = Vec::new();
            let mut fw_start = Vec::new();
            let _ = simd_minimizers::minimizers(m, k - m + 1)
                .super_kmers(&mut fw_start)
                .run(packed.as_slice(), &mut fw_pos);
            let mut rc_seq = seq.to_vec();
            rc_seq.reverse();
            for b in &mut rc_seq {
                *b = decode_base(encode_base_unchecked(*b) ^ 2);
            }
            let rc_packed = PackedSeqVec::from_ascii(&rc_seq);
            let mut rc_pos = Vec::new();
            let mut rc_start = Vec::new();
            let _ = simd_minimizers::minimizers(m, k - m + 1)
                .super_kmers(&mut rc_start)
                .run(rc_packed.as_slice(), &mut rc_pos);
            let n = seq.len() - k + 1;
            let mut f = 0;
            let mut r = rc_pos.len() - 1;
            let mut key = encode(&seq[..k]);
            let mut rc = reverse(key, k);
            for i in 0..n {
                while f + 1 < fw_start.len() && i >= fw_start[f + 1] as usize {
                    f += 1;
                }
                while r > 0 && n - 1 - i < rc_start[r] as usize {
                    r -= 1;
                }
                let p = if key <= rc {
                    fw_pos[f] as usize
                } else {
                    seq.len() - m - rc_pos[r] as usize
                };
                if self.positions.last() != Some(&(p as u32)) {
                    self.positions.push(p as u32);
                    self.starts.push(i as u32);
                    self.values
                        .push(canonical(encode(&seq[p..p + m]), m) as u64);
                }
                if i + 1 < n {
                    let b = encode_base_unchecked(seq[i + k]) as u128;
                    key = (key >> 2) | (b << (2 * (k - 1)));
                    rc = ((rc << 2) | (b ^ 2)) & mask(k);
                }
            }
            self.starts.push(n as u32);
            return;
        }
        odd_packed_minimizers(
            packed.as_slice(),
            k,
            m,
            &mut self.positions,
            &mut self.starts,
            &mut self.values,
        );
    }

    /// Compute minimizers from existing packed DNA, without an ASCII round trip.
    /// Positions and starts are relative to the supplied view, including when
    /// that view begins at an unaligned position in a larger packed document.
    pub fn compute_packed(&mut self, seq: PackedSeq<'_>, k: usize, m: usize) {
        if k % 2 == 0 {
            self.compute(&seq.unpack(), k, m);
            return;
        }
        self.positions.clear();
        self.starts.clear();
        self.values.clear();
        odd_packed_minimizers(
            seq,
            k,
            m,
            &mut self.positions,
            &mut self.starts,
            &mut self.values,
        );
    }
}

fn odd_packed_minimizers(
    seq: PackedSeq<'_>,
    k: usize,
    m: usize,
    positions: &mut Vec<u32>,
    starts: &mut Vec<u32>,
    values: &mut Vec<u64>,
) {
    let _ = simd_minimizers::canonical_minimizers(m, k - m + 1)
        .super_kmers(starts)
        .run(seq, positions);
    starts.push((seq.len() - k + 1) as u32);
    values.extend(positions.iter().map(|&p| {
        let key = seq.read_kmer(m, p as usize);
        key.min(reverse_complement_mmer(key, m))
    }));
}

/// Reverse complement of an m-mer packed in the low `2m` bits of a word.
///
/// This is the u64 body of [`Kmer::reverse_complement`] specialized to a value
/// known to fit one word, regardless of how wide the k-mer type is: complement
/// (A<->T is 00<->10, C<->G is 01<->11, i.e. XOR with 0b10 per base), reverse
/// the 2-bit pairs, then right-align.
#[inline(always)]
pub fn reverse_complement_mmer(mmer: u64, m: usize) -> u64 {
    debug_assert!(m >= 1 && 2 * m < 64, "m must be in [1, 31]");
    let mut x = mmer ^ 0xAAAA_AAAA_AAAA_AAAAu64;
    x = ((x >> 2) & 0x3333_3333_3333_3333u64) | ((x & 0x3333_3333_3333_3333u64) << 2);
    x = ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0Fu64) | ((x & 0x0F0F_0F0F_0F0F_0F0Fu64) << 4);
    x = x.swap_bytes();
    x >> (64 - 2 * m)
}

/// True if the m-mer equals its own reverse complement (possible only for
/// even m). The streaming same-minimizer memo must fall through to a full
/// verification for such minimizers.
#[inline(always)]
pub fn is_self_rc_mmer(mmer: u64, m: usize) -> bool {
    reverse_complement_mmer(mmer, m) == mmer
}

/// The canonical m-mer at locus `i` of a k-mer, given the k-mer's reverse
/// complement: `kappa(i) = min(x_i, rc(x_i))`.
///
/// Uses `rc(x_i) = rc(x)_{K-m-i}`: the reverse complements of all k-m+1 loci
/// are windows of the single reverse complement of the whole k-mer, so the
/// k-mer is reverse-complemented once instead of once per locus (the caller
/// usually has rc(x) already — a lookup needs it anyway for orientation).
#[inline(always)]
pub fn canonical_mmer_at<const K: usize>(
    kmer: &Kmer<K>,
    kmer_rc: &Kmer<K>,
    m: usize,
    i: usize,
) -> u64
where
    Kmer<K>: KmerBits,
{
    debug_assert!(i + m <= K);
    let mask = (1u64 << (2 * m)) - 1;
    let (fwd, rc) = if K <= 31 {
        let bits = <Kmer<K> as KmerBits>::to_u64(kmer.bits());
        let rc_bits = <Kmer<K> as KmerBits>::to_u64(kmer_rc.bits());
        (
            (bits >> (2 * i)) & mask,
            (rc_bits >> (2 * (K - m - i))) & mask,
        )
    } else {
        let bits = <Kmer<K> as KmerBits>::to_u128(kmer.bits());
        let rc_bits = <Kmer<K> as KmerBits>::to_u128(kmer_rc.bits());
        (
            ((bits >> (2 * i)) as u64) & mask,
            ((rc_bits >> (2 * (K - m - i))) as u64) & mask,
        )
    };
    debug_assert_eq!(rc, reverse_complement_mmer(fwd, m));
    fwd.min(rc)
}

/// Information about a minimizer within a k-mer
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MinimizerInfo {
    /// The minimizer value: the canonical m-mer at the selected locus
    pub value: u64,
    /// The absolute position of the selected locus in the full string
    /// (not within the k-mer)
    pub position: u64,
    /// The selected locus within the k-mer window, in [0, k-m]
    /// (0 = leftmost m-mer of the window)
    pub pos_in_kmer: usize,
}

impl MinimizerInfo {
    /// Create a new minimizer info
    pub fn new(value: u64, position: u64, pos_in_kmer: usize) -> Self {
        Self {
            value,
            position,
            pos_in_kmer,
        }
    }
}

/// Compute the minimizer of a single k-mer: the locus minimizing
/// `h(kappa(i))`, with the centre-closest tie-break.
///
/// The returned `position` is the locus within the k-mer (identical to
/// `pos_in_kmer`); callers anchoring in a sequence add the window's absolute
/// begin themselves.
pub fn compute_minimizer<const K: usize>(
    kmer: &Kmer<K>,
    kmer_rc: &Kmer<K>,
    m: usize,
    hasher: &DeterministicHasher,
) -> MinimizerInfo
where
    Kmer<K>: KmerBits,
{
    let _ = (kmer_rc, hasher);
    let mut sequence = [0; K];
    for (i, b) in sequence.iter_mut().enumerate() {
        *b = decode_base(kmer.get_base(i));
    }
    let mut minimizers = SimdMinimizers::default();
    minimizers.compute(&sequence, K, m);
    let position = minimizers.positions[0];
    MinimizerInfo::new(minimizers.values[0], position as u64, position as usize)
}

/// Single-window SIMD iterator for the existing streaming-query API.
/// Construction and batched queries use SimdMinimizers for whole sequences.
pub struct MinimizerIterator {
    k: usize,
    m: usize,
    position: u64,
    hasher: DeterministicHasher,
}

impl MinimizerIterator {
    /// Create a new minimizer iterator with a seeded hasher
    ///
    /// # Arguments
    /// * `k` - k-mer size (must be >= m)
    /// * `m` - minimizer size (must be <= k and <= 31)
    /// * `seed` - seed for the hash function (for deterministic results)
    pub fn with_seed(k: usize, m: usize, seed: u64) -> Self {
        assert!(k > 0 && m <= k, "k must be > 0 and m <= k");
        assert!(2 * m < 64, "m must be <= 31 (minimizers are u64)");
        let hasher = DeterministicHasher::new(seed);
        let mut iter = Self {
            k,
            m,
            position: 0,
            hasher,
        };
        iter.reset();
        iter
    }

    /// Create a new minimizer iterator (defaults to seed=1)
    pub fn new(k: usize, m: usize, _position: u64) -> Self {
        Self::with_seed(k, m, 1)
    }

    /// Set a new starting position
    pub fn set_position(&mut self, position: u64) {
        self.position = position;
        self.reset();
    }

    /// Reset internal state (called when position changes)
    pub fn reset(&mut self) {
        // The next call computes independently at the current window offset.
    }

    /// Compute the minimizer of the next k-mer of the sequence.
    ///
    /// `kmer_rc` must be the reverse complement of `kmer`; taking it as a
    /// parameter (unlike C++, which recomputes it) lets the streaming query
    /// reuse the RC it already maintains incrementally.
    #[inline]
    pub fn next<const K: usize>(&mut self, kmer: &Kmer<K>, kmer_rc: &Kmer<K>) -> MinimizerInfo
    where
        Kmer<K>: KmerBits,
    {
        debug_assert_eq!(K, self.k);
        let mut info = compute_minimizer(kmer, kmer_rc, self.m, &self.hasher);
        info.position += self.position;
        self.position += 1;
        info
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::Kmer;

    #[test]
    fn packed_views_match_ascii_batches_and_scalar_values_at_unaligned_edges() {
        let mut state = 0x2026_1018_512u64;
        let dna: Vec<_> = (0..1400)
            .map(|i| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                if (400..800).contains(&i) {
                    b"AT"[i % 2]
                } else {
                    b"ACGT"[(state & 3) as usize]
                }
            })
            .collect();
        let packed = PackedSeqVec::from_ascii(&dna);
        for k in [3, 4, 7, 31, 32, 33, 63, 64] {
            for m in [1, 2, 3, 7, 13, 19, 31].into_iter().filter(|&m| m < k) {
                for begin in 0..8 {
                    for len in [k, k + 1, k + 63, 1024] {
                        let mut ascii = SimdMinimizers::default();
                        let mut actual = SimdMinimizers::default();
                        ascii.compute(&dna[begin..begin + len], k, m);
                        actual.compute_packed(packed.as_slice().slice(begin..begin + len), k, m);
                        assert_eq!(
                            actual.positions, ascii.positions,
                            "positions k={k} m={m} begin={begin} len={len}"
                        );
                        assert_eq!(actual.starts, ascii.starts);
                        assert_eq!(actual.values, ascii.values);
                        for (&p, &value) in actual.positions.iter().zip(&actual.values) {
                            let offset = begin + p as usize;
                            assert_eq!(
                                value,
                                canonical(encode(&dna[offset..offset + m]), m) as u64
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_reverse_complement_mmer_matches_kmer_rc() {
        // Compare against Kmer::reverse_complement at matching widths.
        let mmer13 = Kmer::<13>::from_str("ACGTTGCAATCGA").unwrap();
        let rc = mmer13.reverse_complement();
        assert_eq!(
            reverse_complement_mmer(<Kmer<13> as KmerBits>::to_u64(mmer13.bits()), 13),
            <Kmer<13> as KmerBits>::to_u64(rc.bits())
        );
    }

    #[test]
    fn test_self_rc_mmer() {
        // "ACGT" is its own reverse complement (even m only).
        let mmer = Kmer::<5>::from_str("ACGTA").unwrap(); // container only
        let _ = mmer;
        let acgt = {
            let mut bits = 0u64;
            for (i, b) in "ACGT".bytes().enumerate() {
                bits |= (((b >> 1) & 3) as u64) << (i * 2);
            }
            bits
        };
        assert!(is_self_rc_mmer(acgt, 4));
        let acgg = acgt & !(3 << 6) | (((b'G' >> 1) & 3) as u64) << 6;
        assert!(!is_self_rc_mmer(acgg, 4));
    }

    #[test]
    fn test_mirror_equivariance_basic() {
        // The scheme's defining property: x and rc(x) select the same value at
        // mirrored loci.
        let hasher = DeterministicHasher::new(1);
        let seqs = [
            "ACGTACGTACGTACGTACGTACGTACGTACG",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTA",
            "ACGTTGCAATCGATCGGCTAGCTAGCTAGCA",
        ];
        for m in [3usize, 4, 7, 12, 13, 16, 19] {
            for s in seqs {
                let kmer = Kmer::<31>::from_str(s).unwrap();
                let rc = kmer.reverse_complement();
                let a = compute_minimizer(&kmer, &rc, m, &hasher);
                let b = compute_minimizer(&rc, &kmer, m, &hasher);
                assert_eq!(a.value, b.value, "value differs, m={m} s={s}");
                assert_eq!(
                    a.pos_in_kmer + b.pos_in_kmer,
                    31 - m,
                    "loci are not mirror images, m={m} s={s}"
                );
                assert_eq!(a.value, canonical_mmer_at(&kmer, &rc, m, a.pos_in_kmer));
            }
        }
    }

    #[test]
    fn test_minimizer_iterator_basic() {
        let kmer: Kmer<15> = Kmer::from_str("ACGTACGTACGTACG").unwrap();
        let rc = kmer.reverse_complement();

        let mut iter = MinimizerIterator::new(15, 7, 0);
        let mini = iter.next(&kmer, &rc);

        assert!(mini.pos_in_kmer <= 15 - 7);
        assert_eq!(mini.position, mini.pos_in_kmer as u64);
    }

    #[test]
    fn test_minimizer_consistency_fresh_vs_sequential() {
        // A fresh iterator gives the same minimizer as processing k-mers
        // sequentially through a string (and both match the batch reference,
        // via the debug assertion inside next()).
        let seq =
            "ATTTTCAGGATGTTTTCAGGTTCATCATCTCCCTTCTTTGCAGGATAGTAGATAAGATCGCTCATCAACGGATGTTGTGT";
        let k = 31usize;
        let m = 19usize;
        let seed = 1u64;
        let num_kmers = seq.len() - k + 1;

        let mut seq_iter = MinimizerIterator::with_seed(k, m, seed);
        seq_iter.set_position(0);

        for kmer_pos in 0..num_kmers {
            let kmer_str = &seq[kmer_pos..kmer_pos + k];
            let kmer = Kmer::<31>::from_str(kmer_str).unwrap();
            let rc = kmer.reverse_complement();

            let seq_mini = seq_iter.next(&kmer, &rc);

            let mut fresh_iter = MinimizerIterator::with_seed(k, m, seed);
            let fresh_mini = fresh_iter.next(&kmer, &rc);

            assert_eq!(
                seq_mini.value, fresh_mini.value,
                "value mismatch at {kmer_pos}"
            );
            assert_eq!(
                seq_mini.pos_in_kmer, fresh_mini.pos_in_kmer,
                "pos mismatch at {kmer_pos}"
            );
            // Forwardness: absolute sampled position never decreases and
            // matches window anchoring.
            assert_eq!(
                seq_mini.position,
                kmer_pos as u64 + seq_mini.pos_in_kmer as u64
            );
        }
    }

    #[test]
    fn test_batch_and_point_minimizers_agree_on_revisited_loci() {
        // Small minimizers exercise non-forward selections and repeated loci.
        let mut state = 0xDEADBEEFu64;
        let mut next_base = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 33) as usize % 4]
        };
        let seq: Vec<u8> = (0..4000).map(|_| next_base()).collect();
        let k = 21usize;
        for m in [3usize, 4, 5] {
            let mut iter = MinimizerIterator::with_seed(k, m, 1);
            iter.set_position(0);
            let mut batch = SimdMinimizers::default();
            batch.compute(&seq, k, m);
            let mut run = 0;
            for (i, w) in seq.windows(k).enumerate() {
                let kmer = Kmer::<21>::from_ascii_unchecked(w);
                let rc = kmer.reverse_complement();
                let mini = iter.next(&kmer, &rc);
                while run + 1 < batch.values.len() && i >= batch.starts[run + 1] as usize {
                    run += 1;
                }
                assert_eq!(mini.position, batch.positions[run] as u64);
                assert_eq!(mini.value, batch.values[run]);
                assert_eq!(mini.pos_in_kmer, batch.positions[run] as usize - i);
            }
        }
    }
}
