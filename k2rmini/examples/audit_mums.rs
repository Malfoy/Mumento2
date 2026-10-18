//! Independent exact uniqueness audit: Aho-Corasick, no locator seed/range code.
use aho_corasick::AhoCorasick;
use clap::Parser;
use helicase::{input::*, *};
use serde::Serialize;
use std::{
    io,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
    thread,
};

const CONFIG: Config = ParserOptions::default().config();
#[derive(Parser)]
struct Args {
    #[arg(long)]
    patterns: PathBuf,
    #[arg(long)]
    input_list: PathBuf,
    #[arg(short = 'j', long, default_value_t = 32)]
    threads: usize,
    #[arg(short, long)]
    output: PathBuf,
}
#[derive(Default, Serialize)]
struct Count {
    absent: usize,
    repeated: usize,
    examples: Vec<(usize, u8)>,
}
fn reverse(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}
fn build(patterns: &[Vec<u8>]) -> io::Result<(AhoCorasick, Vec<usize>)> {
    let mut strings = Vec::new();
    let mut owners = Vec::new();
    for (id, s) in patterns.iter().enumerate() {
        strings.push(s.clone());
        owners.push(id);
        let rc = reverse(s);
        if &rc != s {
            strings.push(rc);
            owners.push(id);
        }
    }
    let ac = AhoCorasick::builder()
        .ascii_case_insensitive(true)
        .build(strings)
        .map_err(io::Error::other)?;
    Ok((ac, owners))
}
fn accumulate(ac: &AhoCorasick, owners: &[usize], sequence: &[u8], counts: &mut [u8]) {
    for hit in ac.find_overlapping_iter(sequence) {
        let n = &mut counts[owners[hit.pattern().as_usize()]];
        *n = (*n + 1).min(2);
    }
}
fn main() -> io::Result<()> {
    let args = Args::parse();
    if args.threads == 0 {
        return Err(io::Error::other("threads must be positive"));
    }
    let paths: Vec<PathBuf> = std::fs::read_to_string(&args.input_list)?
        .lines()
        .map(PathBuf::from)
        .collect();
    let input = FileInput::new(&args.patterns)?;
    let mut names = Vec::new();
    let mut patterns = Vec::new();
    if input.first_byte() != 0 {
        let mut parser = FastxParser::<CONFIG>::from_input(input)?;
        while parser.next().is_some() {
            let pattern = parser.get_dna_string().to_ascii_uppercase();
            if pattern.is_empty() || !pattern.iter().all(|b| b"ACGT".contains(b)) {
                return Err(io::Error::other(
                    "audit patterns must be nonempty unambiguous DNA",
                ));
            }
            names.push(String::from_utf8_lossy(parser.get_header()).into_owned());
            patterns.push(pattern);
        }
    }
    let (ac, owners) = build(&patterns)?;
    let next = AtomicUsize::new(0);
    let finished = AtomicUsize::new(0);
    let results = thread::scope(|scope| -> io::Result<Vec<Vec<Count>>> {
        let mut handles = Vec::new();
        for _ in 0..args.threads.min(paths.len()) {
            let ac = &ac;
            let owners = &owners;
            let paths = &paths;
            let next = &next;
            let finished = &finished;
            let n = patterns.len();
            handles.push(scope.spawn(move || -> io::Result<Vec<Count>> {
                let mut totals: Vec<Count> = (0..n).map(|_| Count::default()).collect();
                let mut counts = vec![0u8; n];
                loop {
                    let id = next.fetch_add(1, Ordering::Relaxed);
                    let Some(path) = paths.get(id) else {
                        break;
                    };
                    counts.fill(0);
                    let input = FileInput::new(path)?;
                    if input.first_byte() != 0 {
                        let mut parser = FastxParser::<CONFIG>::from_input(input)?;
                        while parser.next().is_some() {
                            accumulate(ac, owners, parser.get_dna_string(), &mut counts);
                        }
                    }
                    for (i, &count) in counts.iter().enumerate() {
                        if count != 1 {
                            totals[i].absent += usize::from(count == 0);
                            totals[i].repeated += usize::from(count > 1);
                            if totals[i].examples.len() < 5 {
                                totals[i].examples.push((id, count));
                            }
                        }
                    }
                    let done = finished.fetch_add(1, Ordering::Relaxed) + 1;
                    if done % 500 == 0 || done == paths.len() {
                        eprintln!("[uniqueness_audit] {done}/{} genomes", paths.len());
                    }
                }
                Ok(totals)
            }));
        }
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .map_err(|_| io::Error::other("audit worker panicked"))?
            })
            .collect()
    })?;
    let mut totals: Vec<Count> = (0..patterns.len()).map(|_| Count::default()).collect();
    for partial in results {
        for (dest, source) in totals.iter_mut().zip(partial) {
            dest.absent += source.absent;
            dest.repeated += source.repeated;
            dest.examples.extend(
                source
                    .examples
                    .into_iter()
                    .take(5usize.saturating_sub(dest.examples.len())),
            );
        }
    }
    let report = serde_json::json!({"genomes":paths.len(),"patterns":patterns.len(),"source_passes":1,
        "algorithm":"independent exact overlapping Aho-Corasick; both strands; contigs searched separately; palindromes counted once",
        "all_unique":totals.iter().all(|c|c.absent==0&&c.repeated==0),"names":names,"counts":totals});
    serde_json::to_writer_pretty(std::fs::File::create(args.output)?, &report)?;
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    fn count(pattern: &[u8], contigs: &[&[u8]]) -> u8 {
        let (ac, owners) = build(&[pattern.to_vec()]).unwrap();
        let mut counts = [0];
        for contig in contigs {
            accumulate(&ac, &owners, contig, &mut counts);
        }
        counts[0]
    }
    #[test]
    fn palindrome_is_one_physical_hit() {
        assert_eq!(count(b"ACGT", &[b"ACGT"]), 1);
    }
    #[test]
    fn overlapping_reverse_hits_are_distinct() {
        assert_eq!(count(b"ACG", &[b"ACGT"]), 2);
    }
    #[test]
    fn contigs_never_join() {
        assert_eq!(count(b"ACGT", &[b"AC", b"GT"]), 0);
    }
    #[test]
    fn repeated_contigs_count_twice() {
        assert_eq!(count(b"ACGT", &[b"ACGT", b"ACGT"]), 2);
    }
    #[test]
    fn lowercase_and_ambiguity() {
        assert_eq!(count(b"ACGT", &[b"nnacgtnn"]), 1);
    }
}
