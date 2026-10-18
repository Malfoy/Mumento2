//! Direct library tests using the same reviewed biological scenarios as the CLI
//! suite. Expected substrings come from an independent occurrence-table oracle.
use super::*;
use std::{collections::BTreeSet, sync::OnceLock};

#[derive(Deserialize)]
struct Expected {
    sequence: String,
    hits: Vec<Option<[usize; 2]>>,
}
#[derive(Deserialize)]
struct Step {
    genomes: Vec<Vec<String>>,
    minimum: Option<usize>,
    #[serde(default)]
    fork: bool,
    expected: Vec<Expected>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    kind: String,
    description: String,
    genomes: Vec<Vec<String>>,
    k: usize,
    minimum: Option<usize>,
    forward_only: bool,
    threads: usize,
    chunk_bases: usize,
    m: Option<usize>,
    #[serde(default)]
    expected: Vec<Expected>,
    #[serde(default)]
    updates: Vec<Step>,
    #[serde(default)]
    compression: String,
    #[serde(default)]
    fastx: String,
    #[serde(default)]
    header: String,
    #[serde(default)]
    quality: String,
    wrap: Option<usize>,
    #[serde(default)]
    crlf: bool,
    #[serde(default = "yes")]
    final_newline: bool,
    #[serde(default)]
    seed_mode: String,
    patterns: Option<Vec<String>>,
    error: Option<String>,
    #[serde(default)]
    error_contains: String,
    #[serde(default)]
    oracle_mode: String,
    #[serde(default)]
    header_padding: usize,
    #[serde(default)]
    stable_patterns: bool,
}
fn yes() -> bool {
    true
}
fn all_cases() -> &'static [Case] {
    static CASES: OnceLock<Vec<Case>> = OnceLock::new();
    CASES.get_or_init(|| {
        serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../tests/corpus/cases.json"
        )))
        .unwrap()
    })
}
fn case_options(c: &Case) -> LocateOptions {
    LocateOptions {
        k: c.k,
        m: c.m.unwrap_or(c.k.min(21)),
        threads: c.threads,
        chunk_bases: c.chunk_bases,
        forward_only: c.forward_only,
    }
}
fn reverse(s: &[u8]) -> Vec<u8> {
    s.iter().rev().copied().map(complement).collect()
}
fn name(c: &Case, id: usize) -> String {
    if c.header_padding > 0 {
        return format!("c{id} {}", "x".repeat(c.header_padding));
    }
    match c.header.as_str() {
        "same" => "same".into(),
        "description" => format!("c{id} annotation with spaces"),
        "unicode" => format!("c{id}_é_菌"),
        _ => format!("c{id}"),
    }
}
fn records(c: &Case, genome: &[String]) -> Vec<u8> {
    let mut text = String::new();
    for (id, s) in genome.iter().enumerate() {
        if c.fastx == "fastq" {
            let q = if c.quality.is_empty() {
                "I"
            } else {
                &c.quality
            };
            let quality: String = q.chars().cycle().take(s.len()).collect();
            text.push_str(&format!("@{}\n{s}\n+\n{quality}\n", name(c, id)));
        } else {
            text.push_str(&format!(">{}\n", name(c, id)));
            if s.is_empty() {
                text.push('\n');
            } else {
                for part in s.as_bytes().chunks(c.wrap.unwrap_or(s.len()).max(1)) {
                    text.push_str(std::str::from_utf8(part).unwrap());
                    text.push('\n');
                }
            }
        }
    }
    if !c.final_newline {
        text = text.trim_end_matches('\n').into();
    }
    if c.crlf {
        text = text.replace('\n', "\r\n");
    }
    text.into_bytes()
}
fn compressed(data: &[u8], codec: &str) -> Vec<u8> {
    match codec {
        "gzip" => {
            let mut out = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            out.write_all(data).unwrap();
            out.finish().unwrap()
        }
        "xz" => {
            let mut out = liblzma::write::XzEncoder::new(Vec::new(), 1);
            out.write_all(data).unwrap();
            out.finish().unwrap()
        }
        "zstd" => zstd::stream::encode_all(data, 1).unwrap(),
        "gzip_members" | "xz_members" | "zstd_frames" => {
            let base = match codec {
                "gzip_members" => "gzip",
                "xz_members" => "xz",
                _ => "zstd",
            };
            let split = data.len().min(9);
            let mut output = compressed(&data[..split], base);
            output.extend(compressed(&data[split..], base));
            output
        }
        _ => data.to_vec(),
    }
}
fn write_inputs(dir: &Path, c: &Case, genomes: &[Vec<String>], prefix: &str) -> Vec<PathBuf> {
    genomes
        .iter()
        .enumerate()
        .map(|(g, genome)| {
            let codec = if c.compression == "mixed" {
                ["gzip", "xz", "zstd", "plain"][g % 4]
            } else {
                &c.compression
            };
            let path = dir.join(format!("{prefix}-{g}.fastx"));
            fs::write(&path, compressed(&records(c, genome), codec)).unwrap();
            path
        })
        .collect()
}
fn write_patterns(
    dir: &Path,
    c: &Case,
    genomes: &[Vec<String>],
    prefix: &str,
    union: bool,
) -> PathBuf {
    let mut text = Vec::new();
    let mut add = |s: &[u8]| {
        text.extend_from_slice(b">seed\n");
        text.extend_from_slice(s);
        text.push(b'\n');
    };
    if let Some(patterns) = c.patterns.as_ref().filter(|_| !union) {
        for p in patterns {
            add(p.as_bytes());
        }
    } else {
        for s in genomes.iter().flatten() {
            for run in s
                .as_bytes()
                .split(|&b| !is_dna(b))
                .filter(|s| s.len() >= c.k)
            {
                match if union { "" } else { c.seed_mode.as_str() } {
                    "kmers" => {
                        for kmer in run.windows(c.k) {
                            add(kmer);
                        }
                    }
                    "duplicates" => {
                        for _ in 0..3 {
                            add(run);
                        }
                    }
                    "reverse" => add(&reverse(run)),
                    _ => add(run),
                }
            }
        }
    }
    let path = dir.join(format!("{prefix}-patterns.fa.zst"));
    fs::write(&path, zstd::stream::encode_all(text.as_slice(), 1).unwrap()).unwrap();
    path
}

type Key = (String, Vec<Option<[usize; 2]>>);
fn identical_contigs_oracle(c: &Case, genomes: &[Vec<String>]) -> BTreeSet<Key> {
    let canonical = |s: &str| {
        let reverse = String::from_utf8(reverse(s.as_bytes())).unwrap();
        if c.forward_only || s <= reverse.as_str() {
            s.to_string()
        } else {
            reverse
        }
    };
    let mut first: Vec<_> = genomes[0].iter().map(|s| canonical(s)).collect();
    first.sort();
    for genome in genomes {
        assert!(genome.iter().all(|s| s.bytes().all(is_dna)));
        let mut set: Vec<_> = genome.iter().map(|s| canonical(s)).collect();
        set.sort();
        assert_eq!(set, first, "long-fixture identity premise");
    }
    let mut output = BTreeSet::new();
    for pattern in first.into_iter().collect::<BTreeSet<_>>() {
        if pattern.len() < c.k {
            continue;
        }
        let rc = reverse(pattern.as_bytes());
        let mut hits = Vec::new();
        for genome in genomes {
            let mut positions = BTreeSet::new();
            for (contig, s) in genome.iter().enumerate() {
                if s.len() < pattern.len() {
                    continue;
                }
                for (start, p) in s.as_bytes().windows(pattern.len()).enumerate() {
                    if p == pattern.as_bytes() || (!c.forward_only && p == rc) {
                        positions.insert([contig, start]);
                    }
                }
            }
            hits.push(if positions.len() == 1 {
                positions.first().copied()
            } else {
                None
            });
        }
        if hits.iter().all(Option::is_some) {
            output.insert((pattern, hits));
        }
    }
    output
}
fn verify(
    c: &Case,
    path: &Path,
    genomes: &[Vec<String>],
    minimum: Option<usize>,
    expected: &[Expected],
) -> BTreeSet<Key> {
    let text =
        String::from_utf8(zstd::stream::decode_all(File::open(path).unwrap()).unwrap()).unwrap();
    let mut lines = text.lines();
    assert!(lines.next().unwrap().starts_with("# Coordinates:"));
    assert_eq!(lines.next().unwrap(), "length\tstarts\tstrands\tcontigs");
    let mut observed = BTreeSet::new();
    let mut ordering = Vec::new();
    for row in lines {
        let fields: Vec<_> = row.split('\t').collect();
        assert_eq!(fields.len(), 4, "{row}");
        let len = fields[0].parse::<usize>().unwrap();
        assert!(len >= c.k, "{row}");
        let starts: Vec<_> = fields[1].split(',').collect();
        let strands: Vec<_> = fields[2].split(',').collect();
        let contigs: Vec<_> = fields[3].split(',').collect();
        assert_eq!(
            (starts.len(), strands.len(), contigs.len()),
            (genomes.len(), genomes.len(), genomes.len())
        );
        let mut hits = Vec::new();
        let mut sequence = None;
        for g in 0..genomes.len() {
            if starts[g] == "." {
                assert_eq!((strands[g], contigs[g]), (".", "."));
                hits.push(None);
                continue;
            }
            let start = starts[g].parse::<usize>().unwrap();
            let contig = contigs[g].parse::<usize>().unwrap();
            assert!(contig < genomes[g].len(), "{row}");
            let s = genomes[g][contig].as_bytes();
            assert!(start + len <= s.len(), "{row}");
            let physical: Vec<_> = s[start..start + len]
                .iter()
                .map(u8::to_ascii_uppercase)
                .collect();
            assert!(physical.iter().copied().all(is_dna));
            assert!(["+", "-"].contains(&strands[g]));
            if c.forward_only {
                assert_eq!(strands[g], "+");
            }
            let p = if strands[g] == "-" {
                reverse(&physical)
            } else {
                physical
            };
            if let Some(sequence) = &sequence {
                assert_eq!(&p, sequence, "{row}");
            } else {
                assert_eq!(strands[g], "+");
                sequence = Some(p);
            }
            hits.push(Some([contig, start]));
        }
        assert!(hits.iter().flatten().count() >= minimum.unwrap_or(genomes.len()));
        let sequence = sequence.unwrap();
        let rc = reverse(&sequence);
        let sequence = if c.forward_only || sequence <= rc {
            sequence
        } else {
            rc
        };
        let owner = hits.iter().position(Option::is_some).unwrap();
        let [contig, start] = hits[owner].unwrap();
        let absolute = genomes[owner][..contig]
            .iter()
            .map(String::len)
            .sum::<usize>()
            + start;
        ordering.push((owner, absolute, len));
        assert!(
            observed.insert((String::from_utf8(sequence).unwrap(), hits)),
            "duplicate row: {row}"
        );
    }
    assert!(
        ordering.windows(2).all(|w| w[0] <= w[1]),
        "output coordinate order"
    );
    let wanted: BTreeSet<_> = expected
        .iter()
        .map(|m| (m.sequence.clone(), m.hits.clone()))
        .collect();
    // Check the committed Python occurrence-table oracle against a separate
    // Rust oracle that searches each candidate in all original substrings.
    let normalized: Vec<Vec<String>> = genomes
        .iter()
        .map(|g| g.iter().map(|s| s.to_ascii_uppercase()).collect())
        .collect();
    let independent: BTreeSet<Key> = if c.oracle_mode == "identical_contigs" {
        identical_contigs_oracle(c, &normalized)
    } else {
        super::tests::oracle_with_k(
            &normalized,
            minimum.unwrap_or(genomes.len()),
            c.forward_only,
            c.k,
        )
        .into_iter()
        .map(|(s, hits)| {
            (
                String::from_utf8(s).unwrap(),
                hits.into_iter().map(|p| p.map(|(c, p)| [c, p])).collect(),
            )
        })
        .collect()
    };
    assert_eq!(
        wanted, independent,
        "independent oracles disagree for {}",
        c.name
    );
    assert_eq!(observed, wanted, "{}: {}", c.name, c.description);
    let mut mapping_name = path.as_os_str().to_os_string();
    mapping_name.push(".contigs.tsv.zst");
    let mapping = String::from_utf8(
        zstd::stream::decode_all(File::open(PathBuf::from(mapping_name)).unwrap()).unwrap(),
    )
    .unwrap();
    let mut mapping_lines = mapping.lines();
    assert_eq!(
        mapping_lines.next().unwrap(),
        "genome\tinput\tcontig\tname\tlength"
    );
    let mut actual = BTreeMap::new();
    for row in mapping_lines {
        let f: Vec<_> = row.split('\t').collect();
        assert_eq!(f.len(), 5);
        let (g, id) = (
            f[0].parse::<usize>().unwrap(),
            f[2].parse::<usize>().unwrap(),
        );
        assert!(!f[1].is_empty());
        assert!(
            actual
                .insert((g, id), (f[3].to_string(), f[4].parse::<usize>().unwrap()))
                .is_none()
        );
    }
    let wanted: BTreeMap<_, _> = genomes
        .iter()
        .enumerate()
        .flat_map(|(g, contigs)| {
            contigs.iter().enumerate().map(move |(id, s)| {
                (
                    (g, id),
                    (
                        name(c, id).split_whitespace().next().unwrap().to_string(),
                        s.len(),
                    ),
                )
            })
        })
        .collect();
    assert_eq!(actual, wanted, "contig metadata");
    observed
}
fn discovery(c: &Case, dir: &Path) {
    let inputs = write_inputs(dir, c, &c.genomes, "g");
    let patterns = write_patterns(dir, c, &c.genomes, "g", false);
    let out = dir.join("result.zst");
    let stats = locate_shared_mums(
        &patterns,
        &inputs,
        &out,
        Some(dir),
        case_options(c),
        c.minimum,
        None,
        None,
    )
    .unwrap();
    verify(c, &out, &c.genomes, c.minimum, &c.expected);
    assert_eq!(stats.verified_mums, c.expected.len());
    assert_eq!(stats.source_genomes_read, inputs.len());
    if c.minimum.is_none_or(|n| n == c.genomes.len()) {
        // A graph dictionary has one canonical position per k-mer. Supply a
        // complete spectrum from genome zero, independently of arbitrary proto
        // overlaps/duplicates used by the legacy-path fixture.
        if c.genomes.iter().flatten().map(String::len).sum::<usize>() < 20_000 {
            let mut keys = BTreeSet::new();
            for sequence in &c.genomes[0] {
                for run in sequence
                    .as_bytes()
                    .split(|&b| !is_dna(b))
                    .filter(|r| r.len() >= c.k)
                {
                    for seed in run.windows(c.k) {
                        let mut key = seed.to_vec();
                        if !c.forward_only {
                            key = key.min(reverse(seed));
                        }
                        keys.insert(key);
                    }
                }
            }
            let graph = dir.join("graph-distinct-kmers.fa");
            let mut text = Vec::new();
            for key in keys {
                text.extend_from_slice(b">k\n");
                text.extend(key);
                text.push(b'\n');
            }
            fs::write(&graph, text).unwrap();
            let compact = dir.join("sshash.zst");
            let stats = super::super::locate_graph_mums(
                &graph,
                &inputs,
                &compact,
                Some(dir),
                case_options(c),
            )
            .unwrap();
            assert_eq!(stats.verified_mums, c.expected.len(), "SSHash match count");
            assert_eq!(stats.source_genomes_read, inputs.len());
            verify(c, &compact, &c.genomes, c.minimum, &c.expected);
        }
        let classic = dir.join("classic.zst");
        let stats = locate_mums(&patterns, &inputs, &classic, Some(dir), case_options(c)).unwrap();
        assert_eq!(stats.verified_mums, c.expected.len());
        verify(c, &classic, &c.genomes, c.minimum, &c.expected);
        // Force pruning between tiny batches, including reordered inputs;
        // compare every position/strand against the independent oracle.
        for batch in [0, 1, 2] {
            let batched = dir.join(format!("pruned-{batch}.zst"));
            let stats = super::super::locate_mums_with_batch(
                &patterns,
                &inputs,
                &batched,
                Some(dir),
                case_options(c),
                batch,
            )
            .unwrap();
            assert_eq!(stats.verified_mums, c.expected.len());
            assert_eq!(stats.source_genomes_read, inputs.len());
            verify(c, &batched, &c.genomes, c.minimum, &c.expected);
        }
    }
}
fn update(c: &Case, dir: &Path) {
    let mut genomes = c.genomes.clone();
    let inputs = write_inputs(dir, c, &genomes, "old");
    let patterns = write_patterns(dir, c, &genomes, "old", true);
    let mut state = dir.join("state");
    let out = dir.join("initial.zst");
    locate_shared_mums(
        &patterns,
        &inputs,
        &out,
        Some(dir),
        case_options(c),
        c.minimum,
        Some(&state),
        None,
    )
    .unwrap();
    verify(c, &out, &genomes, c.minimum, &c.expected);
    let first_patterns = read_state_info(&state).unwrap().patterns;
    for path in inputs {
        fs::remove_file(path).unwrap();
    }
    let mut threshold = c.minimum;
    let mut last_keys = BTreeSet::new();
    for (id, step) in c.updates.iter().enumerate() {
        let previous = state.clone();
        let pointer = fs::read(previous.join("CURRENT")).unwrap();
        let inputs = write_inputs(dir, c, &step.genomes, &format!("new-{id}"));
        let patterns = write_patterns(dir, c, &step.genomes, &format!("new-{id}"), true);
        if step.fork {
            state = dir.join(format!("fork-{id}"));
        }
        threshold = step.minimum.or(threshold);
        let out = dir.join(format!("updated-{id}.zst"));
        let stats = locate_shared_mums(
            &patterns,
            &inputs,
            &out,
            Some(dir),
            case_options(c),
            threshold,
            step.fork.then_some(state.as_path()),
            Some(&previous),
        )
        .unwrap();
        genomes.extend(step.genomes.clone());
        last_keys = verify(c, &out, &genomes, threshold, &step.expected);
        assert_eq!(stats.source_genomes_read, inputs.len());
        let info = read_state_info(&state).unwrap();
        assert_eq!(
            (info.genomes, info.minimum_genomes),
            (genomes.len(), threshold)
        );
        assert_eq!(
            (info.k, info.m, info.forward_only),
            (c.k, case_options(c).m, c.forward_only)
        );
        if c.stable_patterns {
            assert_eq!(info.patterns, first_patterns);
        }
        if step.fork {
            assert_eq!(fs::read(previous.join("CURRENT")).unwrap(), pointer);
        }
        for path in inputs {
            fs::remove_file(path).unwrap();
        }
    }
    let fresh = dir.join("fresh");
    fs::create_dir(&fresh).unwrap();
    let inputs = write_inputs(&fresh, c, &genomes, "g");
    let patterns = write_patterns(&fresh, c, &genomes, "g", true);
    let out = fresh.join("rebuilt.zst");
    locate_shared_mums(
        &patterns,
        &inputs,
        &out,
        Some(dir),
        case_options(c),
        threshold,
        None,
        None,
    )
    .unwrap();
    assert_eq!(
        verify(
            c,
            &out,
            &genomes,
            threshold,
            &c.updates.last().unwrap().expected
        ),
        last_keys
    );
}
fn failure(c: &Case, dir: &Path) {
    let error = c.error.as_deref().unwrap();
    let mut inputs = write_inputs(dir, c, &c.genomes, "g");
    let mut patterns = write_patterns(dir, c, &c.genomes, "g", true);
    let mut options = case_options(c);
    let mut minimum = c.minimum;
    let output = dir.join("failure.zst");
    let sentinel = b"preserve previous output";
    fs::write(&output, sentinel).unwrap();
    let mut chosen_output = output.clone();
    let mut source = None;
    let mut dest = None;
    let mut pointer = None;
    let mut lock = None;
    match error {
        "minimum_zero" => minimum = Some(0),
        "minimum_exceeds_documents" => minimum = Some(3),
        "threads_zero" => options.threads = 0,
        "chunk_zero" => options.chunk_bases = 0,
        "minimizer_zero" => options.m = 0,
        "minimizer_exceeds_k" => options.m = options.k + 1,
        "minimizer_exceeds_word" => {
            options.k = 40;
            options.m = 30;
        }
        "missing_input" => inputs = vec![dir.join("absent.fa")],
        "duplicate_input" => inputs = vec![inputs[0].clone(), inputs[0].clone()],
        "missing_state" => source = Some(dir.join("absent")),
        "missing_output_parent" => chosen_output = dir.join("absent/out.zst"),
        _ => {
            let state = dir.join("state");
            let initial = dir.join("initial.zst");
            locate_shared_mums(
                &patterns,
                &inputs,
                &initial,
                Some(dir),
                options.clone(),
                minimum,
                Some(&state),
                None,
            )
            .unwrap();
            let current = generation(&state).unwrap();
            let new = vec![vec!["ACGTTGCAACCTGGA".into()]];
            let original = inputs.clone();
            inputs = write_inputs(dir, c, &new, "new");
            patterns = write_patterns(dir, c, &new, "new", true);
            source = Some(state.clone());
            match error {
                "existing_initial_state" => {
                    source = None;
                    dest = Some(state.clone());
                }
                "update_changed_k" => options.k = 9,
                "update_changed_minimizer" => options.m = 1,
                "update_changed_strand" => options.forward_only = true,
                "update_duplicate_document" => inputs = vec![original[0].clone()],
                "unsupported_state_version" | "zero_state_genomes" => {
                    let path = current.join("state.json");
                    let mut info: serde_json::Value =
                        serde_json::from_reader(File::open(&path).unwrap()).unwrap();
                    if error == "unsupported_state_version" {
                        info["version"] = 999.into();
                    } else {
                        info["genomes"] = 0.into();
                    }
                    fs::write(path, serde_json::to_vec(&info).unwrap()).unwrap();
                }
                "missing_current" => fs::remove_file(state.join("CURRENT")).unwrap(),
                "invalid_current_path" => fs::write(state.join("CURRENT"), "../escape").unwrap(),
                "missing_state_manifest" => fs::remove_file(current.join("state.json")).unwrap(),
                "invalid_state_json" => fs::write(current.join("state.json"), "{broken").unwrap(),
                "invalid_proto_dna" => fs::write(
                    current.join("protos.fa.zst"),
                    zstd::stream::encode_all(b">invalid\nACNNNNT\n".as_slice(), 1).unwrap(),
                )
                .unwrap(),
                "invalid_packed_magic" => fs::write(
                    current.join("genomes/0.zst"),
                    zstd::stream::encode_all(b"badmagic".as_slice(), 1).unwrap(),
                )
                .unwrap(),
                "truncated_context_record" => fs::write(
                    current.join("contexts/0/contexts-0-0.zst"),
                    zstd::stream::encode_all(b"bad-record".as_slice(), 1).unwrap(),
                )
                .unwrap(),
                "missing_proto_archive" => fs::remove_file(current.join("protos.fa.zst")).unwrap(),
                "locked_state" => lock = Some(StateLock::acquire(&state).unwrap()),
                "existing_fork_destination" => {
                    let fork = dir.join("fork");
                    let out = dir.join("fork-initial.zst");
                    let original_patterns = write_patterns(dir, c, &c.genomes, "fork", true);
                    locate_shared_mums(
                        &original_patterns,
                        &original,
                        &out,
                        Some(dir),
                        options.clone(),
                        minimum,
                        Some(&fork),
                        None,
                    )
                    .unwrap();
                    dest = Some(fork);
                }
                _ => panic!("unknown error fixture {error}"),
            }
            pointer = Some((state.clone(), fs::read(state.join("CURRENT")).ok()));
        }
    }
    let result = locate_shared_mums(
        &patterns,
        &inputs,
        &chosen_output,
        Some(dir),
        options,
        minimum,
        dest.as_deref(),
        source.as_deref(),
    );
    let error = result.expect_err(&format!("{} should fail", c.name));
    assert!(
        error.to_string().contains(&c.error_contains),
        "{} failed for the wrong reason: {error}",
        c.name
    );
    drop(lock);
    if let Some((state, before)) = pointer {
        assert_eq!(
            fs::read(state.join("CURRENT")).ok(),
            before,
            "failed update changed CURRENT"
        );
    }
    assert_eq!(
        fs::read(output).unwrap(),
        sentinel,
        "failed operation replaced output"
    );
}
fn run_case(case_name: &str) {
    let c = all_cases().iter().find(|c| c.name == case_name).unwrap();
    let temp = tempfile::tempdir().unwrap();
    match c.kind.as_str() {
        "discovery" => discovery(c, temp.path()),
        "update" => update(c, temp.path()),
        "error" => failure(c, temp.path()),
        _ => panic!("unknown fixture kind"),
    }
}
macro_rules! corpus_test {
    ($name:ident) => {
        #[test]
        fn $name() {
            run_case(stringify!($name));
        }
    };
}
include!("corpus_cases.rs");
