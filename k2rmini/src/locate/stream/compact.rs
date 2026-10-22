//! SSHash indexes the graph's distinct canonical k-mers, not document suffixes.
//! Document occurrences remain exact proto-coordinate intervals.
use super::*;
use sshash_lib::{BuildConfiguration, Dictionary, DictionaryBuilder, Kmer, KmerBits};

pub(super) struct CompactProtoIndex {
    pub(super) dictionary: Dictionary,
}

/// SSHash minimizer length. Independent of `LocateOptions.m` (which is capped
/// at 29 for legacy state). Pass `Some(m)` from `--sshash-m` to override;
/// `None` uses the k-adaptive default.
///
/// Defaults keep the validated k=63 / m=31 pair, and for smaller k prefer a
/// longer minimizer (stronger absent-window skips) than floor(k/2).
pub(super) fn minimizer_length(k: usize, override_m: Option<usize>) -> usize {
    let clamp = |m: usize| {
        let hi = k.saturating_sub(1).min(31).max(1);
        m.clamp(1, hi)
    };
    if let Some(m) = override_m {
        return clamp(m);
    }
    if k == 0 {
        return 1;
    }
    if k == 63 {
        return 31;
    }
    // Larger m improves miss skips at small/mid k. Cap at 23 for k!=63 so
    // dictionary build stays within the common SSHash word path used in tests.
    let hi = k.saturating_sub(1).min(23).max(1);
    let target = if k <= 7 {
        k / 2
    } else if k <= 31 {
        // ~2/3 k, at least floor(k/2), at most hi
        ((k * 2) / 3).max(k / 2)
    } else {
        k / 2
    };
    target.clamp(1, hi)
}

pub(super) fn build_dictionary(
    k: usize,
    threads: usize,
    parent: Option<&Path>,
    strings: Vec<String>,
    override_m: Option<usize>,
) -> io::Result<Dictionary> {
    let scratch = if let Some(parent) = parent {
        std::fs::create_dir_all(parent)?;
        tempfile::Builder::new()
            .prefix("sshash-")
            .tempdir_in(parent)?
    } else {
        tempfile::Builder::new().prefix("sshash-").tempdir()?
    };
    let m = minimizer_length(k, override_m);
    let mut config = BuildConfiguration::new(k, m).map_err(invalid)?;
    config.num_threads = threads;
    config.ram_limit_gib = 4;
    config.tmp_dirname = scratch.path().to_path_buf();
    config.verbose = false;
    eprintln!("[sshash_parameters]\tk={k}\tm={m}");
    DictionaryBuilder::new(config)
        .map_err(invalid)?
        .build_from_sequences(strings)
        .map_err(io::Error::other)
}

/// Bounded SIMD batches retain the graph dictionary's exact lookup and allow
/// proven alignment skips without computing a minimizer separately per seed.
pub(super) struct CompactQuery<'a, const K: usize>
where
    Kmer<K>: KmerBits,
{
    dictionary: &'a Dictionary,
    packed: &'a PackedSeqVec,
    stream: sshash_lib::streaming_query::StreamingQueryEngine<'a, K>,
    minimizers: sshash_lib::minimizer::SimdMinimizers,
    begin: usize,
    end: usize,
    cursor: usize,
    next: Option<usize>,
}

impl<const K: usize> CompactQuery<'_, K>
where
    Kmer<K>: KmerBits,
{
    pub(super) fn reset(&mut self) {
        self.stream.reset_for_packed_jump();
        self.next = None;
    }

    #[inline(always)]
    pub(super) fn lookup(&mut self, position: usize, run_end: usize) -> sshash_lib::LookupResult {
        if self.next.is_some_and(|expected| expected != position) {
            self.stream.reset_for_packed_jump();
        }
        if position < self.begin || position >= self.end {
            self.refill_minimizers(position, run_end);
        }
        let local = position - self.begin;
        while self.cursor + 1 < self.minimizers.values.len()
            && local >= self.minimizers.starts[self.cursor + 1] as usize
        {
            self.cursor += 1;
        }
        let selected = self.begin + self.minimizers.positions[self.cursor] as usize;
        let minimizer = sshash_lib::MinimizerInfo::new(
            self.minimizers.values[self.cursor],
            selected as u64,
            selected - position,
        );
        self.next = Some(position + 1);
        self.stream.lookup_packed_with_minimizer(
            self.packed.as_slice().slice(position..position + K),
            minimizer,
        )
    }
    /// An absent minimizer excludes every window in its current SIMD run.
    /// Stop at the first window whose selected minimizer changes, or at the
    /// batch end. A negative k-mer with a present minimizer permits no skip.
    #[inline(always)]
    pub(super) fn next_after_miss(
        &self,
        position: usize,
        result: &sshash_lib::LookupResult,
    ) -> usize {
        debug_assert!(!result.is_found());
        if result.minimizer_found {
            position + 1
        } else {
            let next = self.begin + self.minimizers.starts[self.cursor + 1] as usize;
            debug_assert!(next > position && next <= self.end);
            next
        }
    }

    // Keep batch construction out of the per-seed instruction path. Inlining
    // the remaining lookup lets the caller consume only the needed result fields.
    #[inline(never)]
    fn refill_minimizers(&mut self, position: usize, run_end: usize) {
        self.begin = position;
        let end = run_end.min(position + 16 * 1024 + K - 1);
        self.minimizers.compute_packed(
            self.packed.as_slice().slice(position..end),
            K,
            self.dictionary.m(),
        );
        self.end = end - K + 1;
        self.cursor = 0;
    }
}

impl FragmentIndex {
    pub(super) fn from_graph(source: ProtoIndex, parent: Option<&Path>) -> io::Result<Self> {
        // Graph path: one canonical position per k-mer via compact SSHash.
        // Even k is supported here. Forward-only / k>63 / empty use Self::new,
        // which builds packed occurrence lists behind SSHash (postings.rs).
        if source.is_empty()
            || source.options.forward_only
            || source.options.k < 3
            || source.options.k > 63
        {
            return Self::new(source);
        }
        let options = source.options.clone();
        let k = options.k;
        let mut strings = Vec::with_capacity(source.len());
        let mut patterns = Vec::with_capacity(source.len());
        let mut position_offsets = vec![0];
        for pattern in &source.patterns {
            let dna: Vec<_> = (pattern.start..pattern.start + pattern.len)
                .map(|i| source.sequences.as_slice().get_ascii(i))
                .collect();
            patterns.push(PackedContig::new(&dna, k));
            position_offsets.push(position_offsets.last().unwrap() + dna.len() - k + 1);
            strings.push(String::from_utf8(dna).map_err(io::Error::other)?);
        }
        let override_m = options.sshash_m;
        drop(source);
        let dictionary = build_dictionary(k, options.threads, parent, strings, override_m)?;
        if dictionary.num_strings() as usize != patterns.len() {
            return Err(invalid("SSHash changed proto-MUM record order/count"));
        }
        let count = *position_offsets.last().unwrap();
        let mut active = vec![u64::MAX; count.div_ceil(64)];
        if count % 64 != 0 {
            *active.last_mut().unwrap() = (1u64 << (count % 64)) - 1;
        }
        eprintln!(
            "[proto_dictionary]\tSSHash Rust\t{} positions\t{} bytes",
            count,
            dictionary.num_bits().div_ceil(8)
        );
        Ok(Self {
            options,
            patterns,
            heads: FxHashMap::default(),
            postings: Vec::new(),
            posting_hashes: Vec::new(),
            position_offsets,
            active,
            unique_runs: Vec::new(),
            probe: Vec::new(),
            compact: Some(CompactProtoIndex { dictionary }),
            posting_dictionary: None,
        })
    }
}

impl CompactProtoIndex {
    pub(super) fn query<'a, const K: usize>(
        &'a self,
        source: &'a PackedSeqVec,
    ) -> CompactQuery<'a, K>
    where
        Kmer<K>: KmerBits,
    {
        CompactQuery {
            dictionary: &self.dictionary,
            packed: source,
            stream: self.dictionary.create_streaming_query::<K>(),
            minimizers: Default::default(),
            begin: 0,
            end: 0,
            cursor: 0,
            next: None,
        }
    }

    pub(super) fn scan(
        &self,
        index: &FragmentIndex,
        query: &PackedContig,
        contig: u32,
        origin: u64,
        accept: impl FnMut(Fragment) -> io::Result<()>,
    ) -> io::Result<()> {
        sshash_lib::dispatch_on_k!(index.options.k, K => self.scan_k::<K>(index, query, contig, origin, accept))
    }

    fn scan_k<const K: usize>(
        &self,
        index: &FragmentIndex,
        query: &PackedContig,
        contig: u32,
        origin: u64,
        mut accept: impl FnMut(Fragment) -> io::Result<()>,
    ) -> io::Result<()>
    where
        Kmer<K>: KmerBits,
    {
        if index.is_empty() {
            return Ok(());
        }
        let mut stream = self.query::<K>(&query.forward);
        for run in &query.runs {
            stream.reset();
            let mut pos = run.start;
            while pos + K <= run.end {
                let found = stream.lookup(pos, run.end);
                if !found.is_found() {
                    pos = stream.next_after_miss(pos, &found);
                    continue;
                }
                let id = found.string_id as usize;
                let r = found.kmer_id_in_string as usize;
                let bit = index.position_offsets[id] + r;
                if index.active[bit / 64] & (1 << (bit % 64)) == 0 {
                    pos += 1;
                    continue;
                }
                let reverse = found.kmer_orientation < 0;
                let proto = &index.patterns[id];
                let q = if reverse { query.len() - pos - K } else { pos };
                let (qforward, qbackward, qb) = if reverse {
                    (
                        &query.reverse,
                        &query.forward,
                        query.len() - run.end..query.len() - run.start,
                    )
                } else {
                    (&query.forward, &query.reverse, run.clone())
                };
                if proto
                    .forward
                    .slice(r..r + K)
                    .cmp_lcp(&qforward.slice(q..q + K))
                    .1
                    != K
                {
                    return Err(invalid("SSHash returned an inexact proto seed"));
                }
                let right = proto
                    .forward
                    .slice(r + K..proto.len())
                    .cmp_lcp(&qforward.slice(q + K..qb.end))
                    .1;
                let left = proto
                    .reverse
                    .slice(proto.len() - r..proto.len())
                    .cmp_lcp(&qbackward.slice(query.len() - q..query.len() - qb.start))
                    .1;
                accept(Fragment {
                    begin: r - left,
                    end: r + K + right,
                    hit: Hit {
                        pattern: id as u32,
                        contig,
                        start: origin + if reverse { pos - right } else { pos - left } as u64,
                        reverse,
                    },
                })?;
                // Every canonical k-mer has one graph position. Interior seeds
                // cannot discover another proto or alignment. A jump resets
                // SSHash's rolling state before the next lookup.
                pos += 1 + if reverse { left } else { right };
                stream.reset();
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_k_reduces_minimizer_length_without_changing_dictionary() {
        for k in 3..=64 {
            let m = minimizer_length(k, None);
            assert!((1..k).contains(&m), "k={k} m={m}");
            assert!(m <= 31);
            assert!(m <= k.saturating_sub(1).min(23) || k == 63);
        }
        assert_eq!(minimizer_length(3, None), 1);
        assert_eq!(minimizer_length(7, None), 3); // tiny-k keeps floor(k/2)
        assert_eq!(minimizer_length(21, None), 14); // 2*21/3 = 14
        assert_eq!(minimizer_length(63, None), 31);
        assert_eq!(minimizer_length(21, Some(10)), 10);
        assert_eq!(minimizer_length(21, Some(100)), 20); // clamped to k-1
    }

    fn batch_matches_point<const K: usize>() -> io::Result<()>
    where
        Kmer<K>: KmerBits,
    {
        let mut state = 0x8123_4567_89ab_cdefu64;
        let dna: Vec<_> = (0..150_000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        let mut config = BuildConfiguration::new(K, 19).map_err(invalid)?;
        config.num_threads = 2;
        config.verbose = false;
        let index = CompactProtoIndex {
            dictionary: DictionaryBuilder::new(config)
                .map_err(invalid)?
                .build_from_sequences(vec![String::from_utf8(dna.clone()).unwrap()])
                .map_err(io::Error::other)?,
        };
        let mut query: Vec<_> = dna
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        query.push(b'N');
        query.extend_from_slice(&dna);
        query.extend(
            std::iter::repeat_n(b"GATTACA".as_slice(), 600)
                .flatten()
                .copied(),
        );
        for p in (4096..query.len()).step_by(16_381) {
            if query[p] != b'N' {
                query[p] = if query[p] == b'A' { b'C' } else { b'A' };
            }
        }
        let packed = PackedContig::new(&query, K);
        let mut stream = index.query::<K>(&packed.forward);
        for run in &packed.runs {
            stream.reset();
            let mut p = run.start;
            let mut checks = 0;
            while p + K <= run.end {
                let found = stream.lookup(p, run.end);
                let kmer = Kmer::<K>::from_str(std::str::from_utf8(&query[p..p + K]).unwrap())
                    .map_err(|error| invalid(error.to_string()))?;
                let expected = index.dictionary.query::<K>(&kmer);
                assert_eq!(found.is_found(), expected.is_found(), "K={K}, pos={p}");
                if found.is_found() {
                    assert_eq!(
                        (
                            found.string_id,
                            found.kmer_id_in_string,
                            found.kmer_orientation
                        ),
                        (
                            expected.string_id,
                            expected.kmer_id_in_string,
                            expected.kmer_orientation
                        ),
                        "K={K}, pos={p}"
                    );
                }
                checks += 1;
                if !found.is_found() {
                    let next = stream.next_after_miss(p, &found);
                    for skipped in p + 1..next {
                        let key = Kmer::<K>::from_ascii_unchecked(&query[skipped..skipped + K]);
                        assert!(
                            !index.dictionary.query::<K>(&key).is_found(),
                            "skipped present seed K={K}, pos={skipped}, span={p}..{next}"
                        );
                    }
                    p = next;
                    continue;
                }
                if checks % 3 == 0 {
                    stream.reset();
                    let revisited = stream.lookup(p, run.end);
                    assert_eq!(revisited.is_found(), expected.is_found());
                    if revisited.is_found() {
                        assert_eq!(
                            (
                                revisited.string_id,
                                revisited.kmer_id_in_string,
                                revisited.kmer_orientation
                            ),
                            (
                                expected.string_id,
                                expected.kmer_id_in_string,
                                expected.kmer_orientation
                            )
                        );
                    }
                }
                p += [1, 1, 2, 3, 7, 31, 63, 127, 1021][checks % 9];
            }
        }
        Ok(())
    }

    #[test]
    fn simd_batches_match_exact_point_queries_across_skips_strands_and_boundaries() -> io::Result<()>
    {
        batch_matches_point::<20>()?;
        batch_matches_point::<21>()?;
        batch_matches_point::<31>()?;
        batch_matches_point::<33>()?;
        batch_matches_point::<61>()?;
        batch_matches_point::<63>()
    }

    #[test]
    fn even_k_uses_sshash_compact_index() -> io::Result<()> {
        let mut source = ProtoIndex {
            options: LocateOptions {
                k: 20,
                min_length: 20,

                m: 15,
                sshash_m: None,
                threads: 2,
                chunk_bases: 127,
                forward_only: false,
                mumemto_format: false,
                binary: false,
            },
            patterns: Vec::new(),
            sequences: PackedSeqVec::default(),
        };
        let dna = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        source.add_pattern(dna)?;
        let index = FragmentIndex::from_graph(source, None)?;
        assert!(
            index.compact.is_some(),
            "even k=20 must stay on the SSHash compact path"
        );
        Ok(())
    }
}
