//! SSHash indexes the graph's distinct canonical k-mers, not document suffixes.
//! Document occurrences remain exact proto-coordinate intervals.
use super::*;
use sshash_lib::{BuildConfiguration, Dictionary, DictionaryBuilder, Kmer, KmerBits};

pub(super) struct CompactProtoIndex {
    pub(super) dictionary: Dictionary,
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
        let end = run_end.min(position + 128 * 1024 + K - 1);
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
        if source.is_empty()
            || source.options.forward_only
            || source.options.k < 3
            || source.options.k > 63
            || source.options.k % 2 == 0
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
        drop(source);
        let scratch = if let Some(parent) = parent {
            std::fs::create_dir_all(parent)?;
            tempfile::Builder::new()
                .prefix("sshash-")
                .tempdir_in(parent)?
        } else {
            tempfile::Builder::new().prefix("sshash-").tempdir()?
        };
        // Longer minimizers reduce failed dictionary work at the measured k=63.
        // Preserve the existing smaller-k configuration.
        let minimizer = if k == 63 { 31 } else { 23.min(k - 1) };
        let mut config = BuildConfiguration::new(k, minimizer).map_err(invalid)?;
        config.num_threads = options.threads;
        config.ram_limit_gib = 4;
        config.tmp_dirname = scratch.path().to_path_buf();
        config.verbose = false;
        let dictionary = DictionaryBuilder::new(config)
            .map_err(invalid)?
            .build_from_sequences(strings)
            .map_err(io::Error::other)?;
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
        })
    }
}

impl CompactProtoIndex {
    pub(super) fn query<'a, const K: usize>(&'a self, source: &'a PackedSeqVec) -> CompactQuery<'a, K>
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
        batch_matches_point::<31>()?;
        batch_matches_point::<33>()?;
        batch_matches_point::<61>()?;
        batch_matches_point::<63>()
    }
}
