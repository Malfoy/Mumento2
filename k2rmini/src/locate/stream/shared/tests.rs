use super::*;
use std::collections::BTreeSet;

fn options(forward_only: bool) -> LocateOptions {
    LocateOptions {
        k: 3,
        min_length: 3,

        m: 2,
        sshash_m: None,
        threads: 3,
        chunk_bases: 2,
        forward_only,
        mumemto_format: false,
        binary: false,
    }
}
fn rc(s: &[u8]) -> Vec<u8> {
    s.iter().rev().copied().map(complement).collect()
}

type MatchKey = (Vec<u8>, Vec<Option<(usize, usize)>>);

// Independent exhaustive substring oracle; no seed, alignment or range logic.
fn oracle(genomes: &[Vec<String>], minimum: usize, forward: bool) -> BTreeSet<MatchKey> {
    oracle_with_k(genomes, minimum, forward, 3)
}

pub(super) fn oracle_with_k(
    genomes: &[Vec<String>],
    minimum: usize,
    forward: bool,
    k: usize,
) -> BTreeSet<MatchKey> {
    let mut candidates = BTreeSet::new();
    for contig in genomes.iter().flatten() {
        for start in 0..contig.len() {
            for end in start + k..=contig.len() {
                let p = &contig.as_bytes()[start..end];
                if !p.iter().copied().all(is_dna) {
                    continue;
                }
                let rev = rc(p);
                candidates.insert(if forward || p <= rev.as_slice() {
                    p.to_vec()
                } else {
                    rev
                });
            }
        }
    }
    let mut output = BTreeSet::new();
    for pattern in candidates {
        let rev = rc(&pattern);
        let mut hits = Vec::new();
        let mut flanks = Vec::new();
        let mut repeated = false;
        for genome in genomes {
            let mut positions = Vec::new();
            for (c, contig) in genome.iter().enumerate() {
                let s = contig.as_bytes();
                if s.len() < pattern.len() {
                    continue;
                }
                for (start, p) in s.windows(pattern.len()).enumerate() {
                    let reverse = p != pattern && !forward && p == rev;
                    if p != pattern && !reverse {
                        continue;
                    }
                    positions.push((c, start));
                    let left = start
                        .checked_sub(1)
                        .and_then(|i| s.get(i))
                        .copied()
                        .filter(|&b| is_dna(b));
                    let right = s.get(start + pattern.len()).copied().filter(|&b| is_dna(b));
                    let (l, r) = if reverse {
                        (right.map(complement), left.map(complement))
                    } else {
                        (left, right)
                    };
                    let mut left_set: BTreeSet<_> = l.into_iter().collect();
                    let mut right_set: BTreeSet<_> = r.into_iter().collect();
                    if !forward && pattern == rev {
                        // A palindrome is one physical occurrence, with two
                        // possible orientations for checking maximality.
                        left_set.extend(right.map(complement));
                        right_set.extend(left.map(complement));
                    }
                    flanks.push((left_set, right_set));
                }
            }
            repeated |= positions.len() > 1;
            hits.push(positions.first().copied());
        }
        if repeated || hits.iter().flatten().count() < minimum {
            continue;
        }
        if flanks[0]
            .0
            .iter()
            .any(|b| flanks.iter().all(|p| p.0.contains(b)))
        {
            continue;
        }
        if flanks[0]
            .1
            .iter()
            .any(|b| flanks.iter().all(|p| p.1.contains(b)))
        {
            continue;
        }
        output.insert((pattern, hits));
    }
    output
}

fn observed(path: &Path, genomes: &[Vec<String>], forward: bool) -> BTreeSet<MatchKey> {
    let text =
        String::from_utf8(zstd::stream::decode_all(File::open(path).unwrap()).unwrap()).unwrap();
    let mut output = BTreeSet::new();
    for line in text.lines().filter(|row| !row.is_empty() && !row.starts_with('#')) {
        let fields: Vec<_> = line.split('\t').collect();
        let len = fields[0].parse::<usize>().unwrap();
        let starts: Vec<_> = fields[1].split(',').collect();
        let strands: Vec<_> = fields[2].split(',').collect();
        let contigs: Vec<_> = fields[3].split(',').collect();
        let mut hits = Vec::new();
        let mut sequence = None;
        for g in 0..genomes.len() {
            if starts[g] == "." {
                assert_eq!(strands[g], ".");
                assert_eq!(contigs[g], ".");
                hits.push(None);
                continue;
            }
            let start = starts[g].parse::<usize>().unwrap();
            let c = contigs[g].parse::<usize>().unwrap();
            let p = &genomes[g][c].as_bytes()[start..start + len];
            let p = if strands[g] == "-" {
                assert!(!forward);
                rc(p)
            } else {
                p.to_vec()
            };
            if let Some(s) = &sequence {
                assert_eq!(&p, s, "{line}");
            } else {
                assert_eq!(strands[g], "+");
                sequence = Some(p);
            }
            hits.push(Some((c, start)));
        }
        let p = sequence.unwrap();
        let rev = rc(&p);
        let canonical = if forward || p <= rev { p } else { rev };
        assert!(output.insert((canonical, hits)), "duplicate row: {line}");
    }
    output
}

fn prepare(dir: &Path, genomes: &[Vec<String>], offset: usize) -> (PathBuf, Vec<PathBuf>) {
    let mut patterns = String::new();
    let inputs = genomes
        .iter()
        .enumerate()
        .map(|(i, contigs)| {
            let text: String = contigs
                .iter()
                .enumerate()
                .map(|(c, s)| format!(">c{c}\n{s}\n"))
                .collect();
            let path = dir.join(format!("g{}.fa", i + offset));
            fs::write(&path, text).unwrap();
            for contig in contigs {
                for run in contig.split('N').filter(|s| s.len() >= 3) {
                    patterns.push_str(&format!(">seed\n{run}\n"));
                }
            }
            path
        })
        .collect();
    let path = dir.join(format!("patterns-{offset}.fa.zst"));
    fs::write(
        &path,
        zstd::stream::encode_all(patterns.as_bytes(), 1).unwrap(),
    )
    .unwrap();
    (path, inputs)
}

#[test]
fn min_length_floor_drops_short_survivors_and_state_stores_it() {
    // Unique blocks of length 5 and 8; seed k=3 so both are discoverable.
    let genomes = vec![
        vec!["AAACCCGG".into(), "TTTTT".into()],
        vec!["AAACCCGG".into(), "TTTTT".into()],
    ];
    let temp = tempfile::tempdir().unwrap();
    let (patterns, inputs) = prepare(temp.path(), &genomes, 0);
    let mut opts = options(true);
    opts.k = 3;
    opts.min_length = 3;
    opts.m = 2;
    let out = temp.path().join("all.zst");
    let state = temp.path().join("state");
    locate_shared_mums(
        &patterns,
        &inputs,
        &out,
        Some(temp.path()),
        opts.clone(),
        None,
        Some(&state),
        None,
    )
    .unwrap();
    let all = observed(&out, &genomes, true);
    assert!(all.iter().any(|(s, _)| s.len() == 5), "{all:?}");
    assert!(all.iter().any(|(s, _)| s.len() == 8), "{all:?}");
    let info = read_state_info(&state).unwrap();
    assert_eq!((info.k, info.min_length, info.version), (3, 3, 3));

    opts.min_length = 6;
    let filtered = temp.path().join("filtered.zst");
    locate_shared_mums(
        &patterns,
        &inputs,
        &filtered,
        Some(temp.path()),
        opts,
        None,
        None,
        None,
    )
    .unwrap();
    let kept = observed(&filtered, &genomes, true);
    assert!(kept.iter().all(|(s, _)| s.len() >= 6), "{kept:?}");
    assert!(kept.iter().any(|(s, _)| s.len() == 8), "{kept:?}");
    assert!(!kept.iter().any(|(s, _)| s.len() == 5), "{kept:?}");
    assert!(kept.len() < all.len());
}

#[test]
fn thresholds_match_exhaustive_substrings_including_missing_first_genome() {
    let mut state = 314159u64;
    let mut dna = |len| {
        (0..len)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize] as char
            })
            .collect::<String>()
    };
    for case in 0..16 {
        let a = dna(7);
        let b = dna(6);
        let extra = dna(5);
        let genomes = vec![
            vec![format!("{a}N{extra}")],
            vec![format!("{b}N{a}"), extra.clone()],
            vec![format!(
                "{b}N{}",
                String::from_utf8(rc(a.as_bytes())).unwrap()
            )],
            vec![format!(
                "{b}N{}",
                if case % 2 == 0 { &b[..4] } else { &a[..5] }
            )],
        ];
        let temp = tempfile::tempdir().unwrap();
        let (patterns, inputs) = prepare(temp.path(), &genomes, 0);
        for forward in [true, false] {
            for minimum in 1..=genomes.len() {
                let out = temp.path().join("matches.zst");
                let stats = locate_shared_mums(
                    &patterns,
                    &inputs,
                    &out,
                    Some(temp.path()),
                    options(forward),
                    Some(minimum),
                    None,
                    None,
                )
                .unwrap();
                let expected = oracle(&genomes, minimum, forward);
                let actual = observed(&out, &genomes, forward);
                assert_eq!(
                    actual, expected,
                    "case {case}, N={minimum}, forward={forward}, {genomes:?}"
                );
                assert_eq!(stats.verified_mums, actual.len());
                assert_eq!(stats.source_genomes_read, genomes.len());
            }
        }
    }
}

#[test]
fn updates_equal_rebuild_without_old_sources_and_discover_new_support() {
    let initial = vec![vec!["AACGATCNGGACCA".into()], vec!["GATCGTNCTTAGG".into()]];
    // GGACCA was below N=2; TGCCTG exists only in the new documents.
    let new = vec![
        vec!["GGACCANTGCCTG".into()],
        vec!["CAGGCANACGATNACCA".into()],
    ];
    for forward in [true, false] {
        for minimum in [None, Some(2)] {
            let temp = tempfile::tempdir().unwrap();
            let state = temp.path().join("state");
            let (patterns, inputs) = prepare(temp.path(), &initial, 0);
            let out = temp.path().join("initial.zst");
            locate_shared_mums(
                &patterns,
                &inputs,
                &out,
                Some(temp.path()),
                options(forward),
                minimum,
                Some(&state),
                None,
            )
            .unwrap();
            assert_eq!(
                observed(&out, &initial, forward),
                oracle(&initial, minimum.unwrap_or(initial.len()), forward)
            );
            for path in inputs {
                fs::remove_file(path).unwrap();
            }
            let mut all = initial.clone();
            for (i, genome) in new.iter().enumerate() {
                let (patterns, inputs) = prepare(temp.path(), std::slice::from_ref(genome), 2 + i);
                all.push(genome.clone());
                let out = temp.path().join(format!("updated-{i}.zst"));
                let stats = locate_shared_mums(
                    &patterns,
                    &inputs,
                    &out,
                    Some(temp.path()),
                    options(forward),
                    minimum,
                    None,
                    Some(&state),
                )
                .unwrap();
                assert_eq!(stats.source_genomes_read, 1);
                assert_eq!(
                    observed(&out, &all, forward),
                    oracle(&all, minimum.unwrap_or(all.len()), forward)
                );
                let info = read_state_info(&state).unwrap();
                assert_eq!(info.genomes, all.len());
                assert_eq!(info.minimum_genomes, minimum);
                fs::remove_file(&inputs[0]).unwrap();
            }
            let fresh = temp.path().join("fresh");
            fs::create_dir(&fresh).unwrap();
            let (patterns, inputs) = prepare(&fresh, &all, 0);
            let out = fresh.join("mums.zst");
            locate_shared_mums(
                &patterns,
                &inputs,
                &out,
                Some(temp.path()),
                options(forward),
                minimum,
                None,
                None,
            )
            .unwrap();
            assert_eq!(
                observed(&out, &all, forward),
                observed(&temp.path().join("updated-1.zst"), &all, forward)
            );
        }
    }
}

#[test]
fn invalid_update_preserves_current_state_and_releases_lock() {
    let temp = tempfile::tempdir().unwrap();
    let genomes = vec![vec!["ACGATC".into()], vec!["ACGATC".into()]];
    let (patterns, inputs) = prepare(temp.path(), &genomes, 0);
    let state = temp.path().join("state");
    let out = temp.path().join("mums.zst");
    locate_shared_mums(
        &patterns,
        &inputs,
        &out,
        None,
        options(false),
        None,
        Some(&state),
        None,
    )
    .unwrap();
    let current = fs::read(state.join("CURRENT")).unwrap();
    let lock = StateLock::acquire(&state).unwrap();
    let error = locate_shared_mums(
        &patterns,
        &inputs,
        &out,
        None,
        options(false),
        None,
        None,
        Some(&state),
    )
    .unwrap_err();
    assert!(error.to_string().contains("locked"));
    drop(lock);
    let mut changed = options(false);
    changed.k = 4;
    assert!(
        locate_shared_mums(
            &patterns,
            &inputs,
            &out,
            None,
            changed,
            None,
            None,
            Some(&state)
        )
        .is_err()
    );
    assert_eq!(fs::read(state.join("CURRENT")).unwrap(), current);
    drop(StateLock::acquire(&state).unwrap());
    // Duplicate documents must not silently increase support.
    assert!(
        locate_shared_mums(
            &patterns,
            &inputs[..1],
            &out,
            None,
            options(false),
            None,
            None,
            Some(&state)
        )
        .is_err()
    );
    assert_eq!(fs::read(state.join("CURRENT")).unwrap(), current);
}

#[test]
fn forked_update_reuses_seeds_without_changing_original_state() {
    let temp = tempfile::tempdir().unwrap();
    let genomes = vec![vec!["ACGATC".into()], vec!["ACGATC".into()]];
    let (patterns, inputs) = prepare(temp.path(), &genomes, 0);
    let state = temp.path().join("state");
    let out = temp.path().join("mums.zst");
    locate_shared_mums(
        &patterns,
        &inputs,
        &out,
        None,
        options(false),
        None,
        Some(&state),
        None,
    )
    .unwrap();
    let old_info = read_state_info(&state).unwrap();
    let current = fs::read(state.join("CURRENT")).unwrap();
    let (patterns, inputs) = prepare(temp.path(), &genomes[..1], 2);
    // Incoming FASTA may be plain while stored proto frames are compressed.
    let plain = temp.path().join("plain.fa");
    fs::write(
        &plain,
        zstd::stream::decode_all(File::open(patterns).unwrap()).unwrap(),
    )
    .unwrap();
    let fork = temp.path().join("fork");
    locate_shared_mums(
        &plain,
        &inputs,
        &out,
        None,
        options(false),
        None,
        Some(&fork),
        Some(&state),
    )
    .unwrap();
    assert_eq!(fs::read(state.join("CURRENT")).unwrap(), current);
    assert_eq!(read_state_info(&state).unwrap().genomes, 2);
    let info = read_state_info(&fork).unwrap();
    assert_eq!(info.genomes, 3);
    assert_eq!(
        info.patterns, old_info.patterns,
        "already covered seeds must not accumulate"
    );
    let text =
        String::from_utf8(zstd::stream::decode_all(File::open(out).unwrap()).unwrap()).unwrap();
    assert!(text.ends_with("6\t0,0,0\t+,+,+\t0,0,0\n"));
}

#[test]
fn empty_genome_and_lowered_threshold_can_recover_previously_absent_matches() {
    let temp = tempfile::tempdir().unwrap();
    let old = vec![vec![], vec!["ACGATC".into()]];
    let (patterns, inputs) = prepare(temp.path(), &old, 0);
    let state = temp.path().join("state");
    let out = temp.path().join("mums.zst");
    let stats = locate_shared_mums(
        &patterns,
        &inputs,
        &out,
        None,
        options(false),
        None,
        Some(&state),
        None,
    )
    .unwrap();
    assert_eq!(stats.verified_mums, 0);
    for path in inputs {
        fs::remove_file(path).unwrap();
    }
    let new = vec![vec!["GATCGT".into()]];
    let (patterns, inputs) = prepare(temp.path(), &new, 2);
    let stats = locate_shared_mums(
        &patterns,
        &inputs,
        &out,
        None,
        options(false),
        Some(2),
        None,
        Some(&state),
    )
    .unwrap();
    assert_eq!(stats.verified_mums, 1);
    let text =
        String::from_utf8(zstd::stream::decode_all(File::open(out).unwrap()).unwrap()).unwrap();
    assert!(text.ends_with("6\t.,0,0\t.,+,-\t.,0,0\n"));
}
