#[cfg(test)]
use std::fs::File;
use std::{
    io,
    path::{Path, PathBuf},
};

use helicase::{input::*, *};
use rustc_hash::{FxHashMap, FxHashSet};
use simd_minimizers::{
    minimizers,
    packed_seq::{PackedSeqVec, Seq, SeqVec},
    seq_hash::{KmerHasher, NtHasher},
};

const CONFIG: Config = ParserOptions::default().config();
const NO_POSTING: u32 = u32::MAX;

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message.into())
}

fn is_dna(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

#[cfg(test)]
fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

#[derive(Clone)]
pub struct LocateOptions {
    pub k: usize,
    pub m: usize,
    pub threads: usize,
    pub chunk_bases: usize,
    pub forward_only: bool,
}

impl LocateOptions {
    pub fn validate(&self) -> io::Result<()> {
        if self.k == 0 || self.m == 0 || self.m > self.k || self.m > 29 {
            return Err(invalid("require k >= m >= 1 and m <= 29"));
        }
        if self.k - self.m + 1 >= 65536 {
            return Err(invalid("SIMD minimizer window must be smaller than 65536"));
        }
        if self.threads == 0 || self.chunk_bases == 0 {
            return Err(invalid("threads and chunk-bases must be positive"));
        }
        Ok(())
    }
}

struct Pattern {
    start: usize,
    len: usize,
}

/// Packed proto-MUM sequences. The streaming index below indexes every k-mer
/// position so a rejected extension cannot hide valid internal fragments.
pub struct ProtoIndex {
    options: LocateOptions,
    patterns: Vec<Pattern>,
    sequences: PackedSeqVec,
}

impl ProtoIndex {
    pub fn from_fasta(path: &Path, options: LocateOptions) -> io::Result<Self> {
        options.validate()?;
        let mut index = Self {
            options,
            patterns: Vec::new(),
            sequences: PackedSeqVec::default(),
        };
        let data = FileInput::new(path)?;
        if data.first_byte() == 0 {
            return Ok(index);
        }
        let mut parser = FastxParser::<CONFIG>::from_input(data)?;
        while parser.next().is_some() {
            index.add_pattern(parser.get_dna_string())?;
        }
        Ok(index)
    }

    pub fn len(&self) -> usize {
        self.patterns.len()
    }
    pub fn is_empty(&self) -> bool {
        self.patterns.is_empty()
    }

    fn add_pattern(&mut self, sequence: &[u8]) -> io::Result<()> {
        if sequence.len() < self.options.k || !sequence.iter().copied().all(is_dna) {
            return Err(invalid(
                "proto-MUMs must contain only A/C/G/T and have length >= k",
            ));
        }
        u32::try_from(self.patterns.len()).map_err(|_| invalid("too many proto-MUMs"))?;
        let start = self.sequences.len();
        self.sequences.push_ascii(sequence);
        self.patterns.push(Pattern {
            start,
            len: sequence.len(),
        });
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Hit {
    pattern: u32,
    contig: u32,
    start: u64,
    reverse: bool,
}

mod stream;
pub use stream::{LocateStats, LocateTimings, locate_mums};

#[cfg(test)]
mod tests {
    use super::*;

    fn run(patterns: &str, genomes: &[&str], chunk_bases: usize) -> (LocateStats, String) {
        run_mode(patterns, genomes, chunk_bases, false)
    }

    fn run_mode(
        patterns: &str,
        genomes: &[&str],
        chunk_bases: usize,
        forward_only: bool,
    ) -> (LocateStats, String) {
        let temp = tempfile::tempdir().unwrap();
        let pattern_path = temp.path().join("protos.fasta.zst");
        std::fs::write(
            &pattern_path,
            zstd::stream::encode_all(patterns.as_bytes(), 1).unwrap(),
        )
        .unwrap();
        let inputs: Vec<_> = genomes
            .iter()
            .enumerate()
            .map(|(i, text)| {
                let path = temp.path().join(format!("genome-{i}.fasta.zst"));
                std::fs::write(&path, zstd::stream::encode_all(text.as_bytes(), 1).unwrap())
                    .unwrap();
                path
            })
            .collect();
        let output = temp.path().join("mums.tsv");
        let stats = locate_mums(
            &pattern_path,
            &inputs,
            &output,
            Some(temp.path()),
            LocateOptions {
                k: 3,
                m: 2,
                threads: 3,
                chunk_bases,
                forward_only,
            },
        )
        .unwrap();
        (
            stats,
            String::from_utf8(zstd::stream::decode_all(File::open(output).unwrap()).unwrap())
                .unwrap(),
        )
    }

    #[test]
    fn all_genomes_chimeras_chunk_boundaries_and_duplicate_extensions() {
        for chunk in [1, 4, 1024] {
            let (stats, output) = run(
                ">valid\nACGATC\n>chimera\nACGAAA\n>same_mum\nCGATC\n",
                &[">a\nTACGATCA\n>b\nACG\n>c\nAAA\n", ">z\nGGATCGTC\n"],
                chunk,
            );
            assert_eq!(stats.verified_mums, 1);
            assert_eq!(stats.source_genomes_read, 2);
            assert!(output.ends_with("6\t1,1\t+,-\t0,0\n"), "{output}");
        }
    }

    #[test]
    fn repeated_full_pattern_is_rejected() {
        let (stats, _) = run(">a\nACGATC\n", &[">a\nACGATCNNACGATC\n", ">b\nACGATC\n"], 4);
        assert_eq!(stats.verified_mums, 0);
        assert_eq!(stats.source_genomes_read, 2);
    }

    #[test]
    fn repeated_seeds_extend_to_unique_complete_matches() {
        for chunk in [1, 4, 1024] {
            let (stats, output) = run(
                ">repeated_seed\nACG\n",
                &[">a\nACGATCNNACGAAA\n", ">b\nACGATCNNACGCCC\n"],
                chunk,
            );
            assert_eq!(stats.verified_mums, 1, "{output}");
            assert!(output.ends_with("6\t0,0\t+,+\t0,0\n"), "{output}");
            let (stats, output) = run(
                ">repeated_seed\nACG\n",
                &[">a\nACGATCNNACGAAA\n", ">b\nACGAAANNACGATC\n"],
                chunk,
            );
            assert_eq!(stats.verified_mums, 2, "{output}");
            assert!(output.contains("6\t0,8\t+,+\t0,0\n"), "{output}");
            assert!(output.contains("6\t8,0\t+,+\t0,0\n"), "{output}");
        }
    }

    #[test]
    fn repeated_seed_disambiguation_uses_both_flanks_and_strands() {
        let (stats, output) = run(
            ">seed\nACG\n",
            &[">a\nACGATCNNACGAAA\n", ">b\nGATCGTNNACGCCC\n"],
            1,
        );
        assert_eq!(stats.verified_mums, 1, "{output}");
        assert!(output.ends_with("6\t0,0\t+,-\t0,0\n"), "{output}");
        let (stats, output) = run(
            ">seed\nACG\n",
            &[
                ">a\nTACGATNNTACGACNNGACGAT\n",
                ">b\nTACGATNNTACGAGNNCACGAT\n",
            ],
            2,
        );
        assert_eq!(stats.verified_mums, 1, "{output}");
        assert!(output.ends_with("6\t0,0\t+,+\t0,0\n"), "{output}");
    }

    #[test]
    fn palindromic_seed_disambiguates_after_asymmetric_extension() {
        let (stats, output) = run(
            ">palindrome\nACGT\n",
            &[">a\nAACGTCC\n", ">b\nGGACGTT\n"],
            1,
        );
        assert_eq!(stats.verified_mums, 1, "{output}");
        assert!(output.ends_with("7\t0,0\t+,-\t0,0\n"), "{output}");
    }

    #[test]
    fn repeated_seeds_match_independent_exhaustive_mum_oracle() {
        // All length-three seed strings are provided, including repeats. The
        // oracle enumerates actual substrings and counts their complete DNA,
        // without using the locator's context/extension algorithm.
        let mut state = 917u64;
        let mut dna = |len| -> String {
            (0..len)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(state >> 62) as usize] as char
                })
                .collect()
        };
        for case in 0..32 {
            let a = dna(7);
            let b = dna(5);
            let extra = dna(4);
            let sequences = [
                format!("{a}N{b}N{}", &a[..3]),
                format!("{b}N{a}N{extra}"),
                format!("{a}N{b}N{}", if case % 2 == 0 { &a[..4] } else { &extra }),
            ];
            let mut expected = FxHashSet::default();
            let reference = sequences[0].as_bytes();
            for start in 0..reference.len() {
                for end in start + 3..=reference.len() {
                    let pattern = &reference[start..end];
                    if !pattern.iter().copied().all(is_dna) {
                        continue;
                    }
                    let positions: Vec<Vec<_>> = sequences
                        .iter()
                        .map(|s| {
                            s.as_bytes()
                                .windows(pattern.len())
                                .enumerate()
                                .filter_map(|(i, p)| (p == pattern).then_some(i))
                                .collect()
                        })
                        .collect();
                    if !positions.iter().all(|p| p.len() == 1) {
                        continue;
                    }
                    let starts: Vec<_> = positions.iter().map(|p| p[0]).collect();
                    let extends = |left: bool| {
                        let flanks: Vec<_> = sequences
                            .iter()
                            .zip(&starts)
                            .map(|(s, &p)| {
                                let i = if left {
                                    p.checked_sub(1)
                                } else {
                                    Some(p + pattern.len())
                                };
                                i.and_then(|i| s.as_bytes().get(i)).copied()
                            })
                            .collect();
                        flanks[0].is_some_and(is_dna) && flanks.iter().all(|&f| f == flanks[0])
                    };
                    if !extends(true) && !extends(false) {
                        expected.insert((pattern.len(), starts));
                    }
                }
            }
            let seeds: std::collections::BTreeSet<_> = reference
                .windows(3)
                .filter(|w| w.iter().copied().all(is_dna))
                .map(|w| String::from_utf8(w.to_vec()).unwrap())
                .collect();
            let patterns: String = seeds
                .iter()
                .enumerate()
                .map(|(i, s)| format!(">s{i}\n{s}\n"))
                .collect();
            let fastas: Vec<_> = sequences
                .iter()
                .map(|s| format!(">contig\n{s}\n"))
                .collect();
            let inputs: Vec<_> = fastas.iter().map(String::as_str).collect();
            let (_, output) = run_mode(&patterns, &inputs, 2, true);
            let observed: FxHashSet<_> = output
                .lines()
                .skip(2)
                .map(|line| {
                    let fields: Vec<_> = line.split('\t').collect();
                    (
                        fields[0].parse::<usize>().unwrap(),
                        fields[1]
                            .split(',')
                            .map(|s| s.parse::<usize>().unwrap())
                            .collect::<Vec<_>>(),
                    )
                })
                .collect();
            assert_eq!(observed, expected, "case {case}: {sequences:?}\n{output}");
        }
    }

    #[test]
    fn ambiguity_breaks_matches_and_preserves_offsets() {
        let (stats, output) = run(
            ">a\nACGATC\n",
            &[">a\nACGNNATC\n>b\nNNtacgatca\n", ">b\nNNggatcgtc\n"],
            2,
        );
        assert_eq!(stats.verified_mums, 1);
        assert!(output.ends_with("6\t3,3\t+,-\t1,0\n"), "{output}");
    }

    #[test]
    fn empty_genome_rejects_every_pattern() {
        let (stats, _) = run(">a\nACGATC\n", &[">a\nACGATC\n", ">empty\n"], 1);
        assert_eq!(stats.verified_mums, 0);
        assert_eq!(stats.source_genomes_read, 2);
    }

    #[test]
    fn compressed_empty_patterns_and_genomes_are_valid() {
        let (stats, output) = run("", &["", ""], 1);
        assert_eq!(stats.proto_mums, 0);
        assert_eq!(stats.verified_mums, 0);
        assert_eq!(output.lines().count(), 2);
    }

    #[test]
    fn palindromic_physical_occurrence_is_counted_once() {
        let (stats, output) = run(">a\nACGT\n", &[">a\nACGT\n", ">b\nACGT\n"], 1);
        assert_eq!(stats.verified_mums, 1);
        assert!(output.ends_with("4\t0,0\t+,+\t0,0\n"), "{output}");
    }

    #[test]
    fn invalid_whole_proto_retains_prefix_interior_and_suffix_ranges() {
        let proto = "ACGATCTGCATGACCTAGTC";
        for chunk in [1, 4, 1024] {
            for reversed in [false, true] {
                let pieces = [&proto[..9], &proto[7..15], &proto[13..]];
                let query: String = pieces
                    .iter()
                    .enumerate()
                    .map(|(i, s)| {
                        let dna: String = if reversed {
                            s.bytes().rev().map(|b| complement(b) as char).collect()
                        } else {
                            s.to_string()
                        };
                        format!(">c{i}\n{dna}\n")
                    })
                    .collect();
                let reference = format!(">reference\n{proto}\n");
                let patterns = format!(">proto\n{proto}\n");
                let (stats, output) = run(&patterns, &[&reference, &query], chunk);
                assert_eq!(stats.verified_mums, 3, "{output}");
                let strand = if reversed { '-' } else { '+' };
                for (start, len, contig) in [(0, 9, 0), (7, 8, 1), (13, 7, 2)] {
                    assert!(
                        output.contains(&format!("{len}\t{start},0\t+,{strand}\t0,{contig}\n")),
                        "{output}"
                    );
                }
            }
        }
    }

    #[test]
    fn fragments_match_exhaustive_oracle_with_both_strands_and_contig_breaks() {
        let mut state = 20391u64;
        for case in 0..48 {
            let mut dna = |len| -> Vec<u8> {
                (0..len)
                    .map(|_| {
                        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                        b"ACGT"[(state >> 62) as usize]
                    })
                    .collect()
            };
            let core = dna(13);
            let extra = dna(6);
            let reverse = |s: &[u8]| s.iter().rev().copied().map(complement).collect::<Vec<_>>();
            let genomes = [
                vec![core.clone(), extra.clone(), core[..3].to_vec()],
                vec![core[..8].to_vec(), reverse(&core[6..]), extra.clone()],
                vec![
                    if case % 2 == 0 {
                        reverse(&core)
                    } else {
                        core.clone()
                    },
                    extra,
                ],
            ];
            let mut expected = FxHashSet::default();
            for (cid, reference) in genomes[0].iter().enumerate() {
                for start in 0..reference.len() {
                    for end in start + 3..=reference.len() {
                        let pattern = &reference[start..end];
                        let rc = reverse(pattern);
                        let positions: Vec<Vec<_>> = genomes
                            .iter()
                            .map(|genome| {
                                let mut hits = Vec::new();
                                for (c, seq) in genome.iter().enumerate() {
                                    if seq.len() < pattern.len() {
                                        continue;
                                    }
                                    for (p, w) in seq.windows(pattern.len()).enumerate() {
                                        if w == pattern || w == rc {
                                            hits.push((c, p, w != pattern));
                                        }
                                    }
                                }
                                hits
                            })
                            .collect();
                        if !positions.iter().all(|p| p.len() == 1) {
                            continue;
                        }
                        assert_eq!((positions[0][0].0, positions[0][0].1), (cid, start));
                        let extends = |left: bool| {
                            let flanks: Vec<FxHashSet<u8>> = genomes
                                .iter()
                                .zip(&positions)
                                .enumerate()
                                .map(|(gidx, (g, p))| {
                                    let (c, start, rev) = p[0];
                                    let orientations: &[bool] = if pattern == rc && gidx != 0 {
                                        &[false, true]
                                    } else if rev {
                                        &[true]
                                    } else {
                                        &[false]
                                    };
                                    orientations
                                        .iter()
                                        .filter_map(|&rev| {
                                            let index = if left != rev {
                                                start.checked_sub(1)
                                            } else {
                                                Some(start + pattern.len())
                                            };
                                            index
                                                .and_then(|i| g[c].get(i))
                                                .copied()
                                                .map(|b| if rev { complement(b) } else { b })
                                        })
                                        .collect()
                                })
                                .collect();
                            flanks[0]
                                .iter()
                                .any(|b| flanks.iter().all(|set| set.contains(b)))
                        };
                        if !extends(true) && !extends(false) {
                            expected.insert((
                                pattern.len(),
                                positions.iter().map(|p| p[0]).collect::<Vec<_>>(),
                            ));
                        }
                    }
                }
            }
            // A long proto containing the core plus a separate extra proto.
            // The core is absent as a complete string from genome 1.
            let text = |s: &[u8]| String::from_utf8(s.to_vec()).unwrap();
            let patterns = format!(">long\n{}\n>extra\n{}\n", text(&core), text(&genomes[0][1]));
            let fastas: Vec<_> = genomes
                .iter()
                .map(|g| {
                    g.iter()
                        .enumerate()
                        .map(|(c, s)| format!(">c{c}\n{}\n", text(s)))
                        .collect::<String>()
                })
                .collect();
            let inputs: Vec<_> = fastas.iter().map(String::as_str).collect();
            let (_, output) = run(&patterns, &inputs, 1 + case % 5);
            let observed: FxHashSet<_> = output
                .lines()
                .skip(2)
                .map(|row| {
                    let f: Vec<_> = row.split('\t').collect();
                    let starts: Vec<usize> = f[1].split(',').map(|v| v.parse().unwrap()).collect();
                    let strands: Vec<_> = f[2].split(',').collect();
                    let contigs: Vec<usize> = f[3].split(',').map(|v| v.parse().unwrap()).collect();
                    (
                        f[0].parse::<usize>().unwrap(),
                        (0..genomes.len())
                            .map(|g| (contigs[g], starts[g], strands[g] == "-"))
                            .collect::<Vec<_>>(),
                    )
                })
                .collect();
            assert_eq!(observed, expected, "case {case}: {genomes:?}\n{output}");
        }
    }
}
