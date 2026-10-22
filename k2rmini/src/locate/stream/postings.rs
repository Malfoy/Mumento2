//! SSHash plus packed occurrence lists for overlapping/union proto sets.
//! Canonical seeds are deduplicated by sorting, then stored as contiguous
//! proto runs. Small k changes the minimizer length, never the dictionary kind.
use super::*;
use sshash_lib::{Kmer, KmerBits};

pub(super) struct PostingDictionary {
    pub(super) compact: compact::CompactProtoIndex,
    pub(super) heads: Vec<u32>,
    pub(super) slots: Vec<u32>,
}

impl PostingDictionary {
    pub(super) fn build(index: &FragmentIndex) -> io::Result<Option<Self>> {
        if index.postings.is_empty() {
            return Ok(None);
        }
        sshash_lib::dispatch_on_k!(index.options.k, K => Self::build_k::<K>(index).map(Some))
    }

    fn build_k<const K: usize>(index: &FragmentIndex) -> io::Result<Self>
    where
        Kmer<K>: KmerBits,
    {
        let count = *index.position_offsets.last().unwrap();
        let mut entries = Vec::with_capacity(count);
        let mut slots = vec![0u32; index.postings.len()];
        let strands = if index.options.forward_only { 1 } else { 2 };
        for (id, proto) in index.patterns.iter().enumerate() {
            let positions = index.position_offsets[id + 1] - index.position_offsets[id];
            let base = index.position_offsets[id] * strands;
            for offset in 0..positions {
                let window = proto.forward.slice(offset..offset + K);
                let bits = if K <= 32 {
                    window.as_u64() as u128
                } else {
                    window.slice(0..32).as_u64() as u128
                        | ((window.slice(32..K).as_u64() as u128) << 64)
                };
                let key = Kmer::<K>::from_bits(bits);
                let rc = key.reverse_complement();
                let reverse = key > rc;
                let canonical = if reverse { rc } else { key };
                let forward = base + offset;
                slots[forward] = u32::from(reverse);
                if strands == 2 {
                    // A palindromic seed belongs to both oriented lists at
                    // the same canonical key. Neither occurrence may be lost.
                    slots[base + 2 * positions - 1 - offset] = u32::from(key != rc && !reverse);
                }
                entries.push((
                    canonical.bits(),
                    u32::try_from(forward).map_err(|_| invalid("too many proto positions"))?,
                ));
            }
        }
        entries.sort_unstable_by_key(|entry| entry.0);
        let mut selected = vec![0u64; count.div_ceil(64)];
        let mut previous = None;
        let mut distinct = 0usize;
        for &(key, posting) in &entries {
            if previous != Some(key) {
                distinct += 1;
                let p = &index.postings[posting as usize];
                let bit = index.position_offsets[p.pattern as usize] + p.offset as usize;
                selected[bit / 64] |= 1 << (bit % 64);
                previous = Some(key);
            }
        }
        // Preserve contiguous proto sequence rather than allocating a string
        // per k-mer. Every canonical seed occurs exactly once in this SPSS.
        let mut strings = Vec::new();
        for (id, proto) in index.patterns.iter().enumerate() {
            let base = index.position_offsets[id];
            let positions = index.position_offsets[id + 1] - base;
            let included = |p: usize| selected[(base + p) / 64] & (1 << ((base + p) % 64)) != 0;
            let mut pos = 0;
            while pos < positions {
                if !included(pos) {
                    pos += 1;
                    continue;
                }
                let begin = pos;
                while pos < positions && included(pos) {
                    pos += 1;
                }
                let dna: Vec<_> = (begin..pos + K - 1)
                    .map(|i| proto.forward.as_slice().get_ascii(i))
                    .collect();
                strings.push(String::from_utf8(dna).map_err(io::Error::other)?);
            }
        }
        drop(selected);
        let dictionary = compact::build_dictionary(
            K,
            index.options.threads,
            None,
            strings,
            index.options.sshash_m,
        )?;
        let head_count = distinct
            .checked_mul(2)
            .filter(|&n| n <= u32::MAX as usize)
            .ok_or_else(|| invalid("too many distinct proto seeds"))?;
        previous = None;
        let mut base = 0u32;
        let mut flip = false;
        for (key, posting_id) in entries {
            if previous != Some(key) {
                let found = dictionary.query::<K>(&Kmer::<K>::new(key));
                if !found.is_found() {
                    return Err(invalid("SSHash lost an indexed proto seed"));
                }
                base = u32::try_from(found.kmer_id * 2)
                    .map_err(|_| invalid("too many distinct proto seeds"))?;
                flip = found.kmer_orientation < 0;
                previous = Some(key);
            }
            let posting_id = posting_id as usize;
            let p = &index.postings[posting_id];
            slots[posting_id] = base + (slots[posting_id] ^ u32::from(flip));
            if strands == 2 {
                let pattern = p.pattern as usize;
                let positions =
                    index.position_offsets[pattern + 1] - index.position_offsets[pattern];
                let reverse =
                    index.position_offsets[pattern] * 2 + 2 * positions - 1 - p.offset as usize;
                slots[reverse] = base + (slots[reverse] ^ u32::from(flip));
            }
        }
        eprintln!(
            "[proto_dictionary]\tSSHash Rust\t{} distinct seeds\t{} bytes\tpacked occurrence lists",
            distinct,
            dictionary.num_bits().div_ceil(8)
        );
        Ok(Self {
            compact: compact::CompactProtoIndex { dictionary },
            heads: vec![NO_POSTING; head_count],
            slots,
        })
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
        let mut stream = self.compact.query::<K>(&query.forward);
        let mut covered: FxHashMap<(u32, bool, i64), (usize, usize)> = FxHashMap::default();
        for run in &query.runs {
            stream.reset();
            let mut pos = run.start;
            while pos + K <= run.end {
                let found = stream.lookup(pos, run.end);
                if !found.is_found() {
                    pos = stream.next_after_miss(pos, &found);
                    continue;
                }
                let head = self.heads
                    [found.kmer_id as usize * 2 + usize::from(found.kmer_orientation < 0)];
                if head == NO_POSTING {
                    pos += 1;
                    continue;
                }
                let unique_run = if index.postings[head as usize].next == NO_POSTING {
                    index.unique_runs[head as usize] as usize
                } else {
                    0
                };
                let mut advance = 1;
                let mut next = head;
                while next != NO_POSTING {
                    let posting = &index.postings[next as usize];
                    next = posting.next;
                    let r = posting.offset as usize;
                    let q = if posting.reverse {
                        query.len() - pos - K
                    } else {
                        pos
                    };
                    let key = (posting.pattern, posting.reverse, q as i64 - r as i64);
                    if let Some(&(a, b)) = covered.get(&key) {
                        if a <= r && r + K <= b {
                            if unique_run > 1 {
                                let remaining = if posting.reverse {
                                    r - a + 1
                                } else {
                                    b - r - K + 1
                                };
                                advance = advance.max(unique_run.min(remaining));
                            }
                            continue;
                        }
                    }
                    let proto = &index.patterns[posting.pattern as usize];
                    let (qforward, qbackward, qb) = if posting.reverse {
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
                        continue;
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
                    let (begin, end) = (r - left, r + K + right);
                    covered.insert(key, (begin, end));
                    let physical_start = if posting.reverse {
                        pos - right
                    } else {
                        pos - left
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
                            end - r - K + 1
                        };
                        advance = advance.max(unique_run.min(remaining));
                    }
                }
                pos += advance;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source(k: usize, forward_only: bool, dna: &[Vec<u8>]) -> io::Result<ProtoIndex> {
        let mut source = ProtoIndex {
            options: LocateOptions {
                k,
                min_length: k,
                m: k.min(21),
                sshash_m: None,
                threads: 1,
                chunk_bases: 7,
                forward_only,
                mumemto_format: false,
                binary: false,
            },
            patterns: Vec::new(),
            sequences: PackedSeqVec::default(),
        };
        for pattern in dna {
            source.add_pattern(pattern)?;
        }
        Ok(source)
    }

    #[test]
    fn all_supported_k_keep_sshash_for_graph_forward_and_overlapping_protos() -> io::Result<()> {
        for k in 3..=64 {
            for forward in [false, true] {
                let dna = vec![b"ACGT".repeat(32), b"ACGT".repeat(24)];
                let index = FragmentIndex::new(source(k, forward, &dna)?)?;
                assert!(
                    index.heads.is_empty()
                        && index.posting_hashes.is_empty()
                        && index.probe.is_empty(),
                    "hash index at k={k}"
                );
                let dictionary = index
                    .posting_dictionary
                    .as_ref()
                    .expect("SSHash occurrence dictionary");
                assert_eq!(
                    dictionary.compact.dictionary.m(),
                    compact::minimizer_length(k, None)
                );
                for (id, &slot) in dictionary.slots.iter().enumerate() {
                    let mut next = dictionary.heads[slot as usize];
                    let mut found = false;
                    while next != NO_POSTING {
                        found |= next as usize == id;
                        next = index.postings[next as usize].next;
                    }
                    assert!(
                        found,
                        "lost proto occurrence k={k} forward={forward} id={id}"
                    );
                }
                // A single k-mer is a valid distinct graph, including even k
                // and forward-only mode, which formerly selected hash postings.
                let graph =
                    FragmentIndex::from_graph(source(k, forward, &[dna[0][..k].to_vec()])?, None)?;
                assert!(graph.compact.is_some() || graph.posting_dictionary.is_some());
                assert!(graph.heads.is_empty() && graph.posting_hashes.is_empty());
            }
        }
        Ok(())
    }

    #[test]
    fn sshash_occurrences_match_basewise_fragments_with_even_k_palindromes_and_ambiguity()
    -> io::Result<()> {
        let rc = |s: &[u8]| {
            s.iter()
                .rev()
                .map(|&b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    _ => b'N',
                })
                .collect::<Vec<_>>()
        };
        let mut rng = 573188u64;
        let core: Vec<_> = (0..130)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(rng >> 62) as usize]
            })
            .collect();
        for k in [3, 4, 7, 20, 31, 32, 63, 64] {
            let patterns = vec![
                core.clone(),
                core[20..110].to_vec(),
                rc(&core),
                b"ACGT".repeat(32),
            ];
            let mut broken = core.clone();
            broken[50] = b'N';
            for forward in [false, true] {
                let index = FragmentIndex::new(source(k, forward, &patterns)?)?;
                for query in [&core, &rc(&core), &broken, &patterns[3]] {
                    let mut actual = std::collections::BTreeSet::new();
                    index.scan(&PackedContig::new(query, k), 0, 0, |f| {
                        actual.insert((f.begin, f.end, f.hit.pattern, f.hit.start, f.hit.reverse));
                        Ok(())
                    })?;
                    let mut expected = std::collections::BTreeSet::new();
                    for (id, pattern) in patterns.iter().enumerate() {
                        for reverse in [false, true] {
                            if reverse && forward {
                                continue;
                            }
                            let oriented = if reverse { rc(query) } else { query.to_vec() };
                            for r in 0..=pattern.len() - k {
                                for q in 0..=oriented.len() - k {
                                    if pattern[r..r + k] != oriented[q..q + k] {
                                        continue;
                                    }
                                    let mut left = 0;
                                    while left < r
                                        && left < q
                                        && pattern[r - left - 1] == oriented[q - left - 1]
                                    {
                                        left += 1;
                                    }
                                    let mut right = 0;
                                    while r + k + right < pattern.len()
                                        && q + k + right < oriented.len()
                                        && pattern[r + k + right] == oriented[q + k + right]
                                    {
                                        right += 1;
                                    }
                                    let physical = if reverse {
                                        query.len() - q - k - right
                                    } else {
                                        q - left
                                    };
                                    expected.insert((
                                        r - left,
                                        r + k + right,
                                        id as u32,
                                        physical as u64,
                                        reverse,
                                    ));
                                }
                            }
                        }
                    }
                    assert_eq!(actual, expected, "k={k} forward={forward}");
                }
            }
        }
        Ok(())
    }
}
