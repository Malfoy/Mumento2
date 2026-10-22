use std::{
    io::Write,
    path::{Path, PathBuf},
    process::Command,
    thread,
    time::Instant,
};

use anyhow::{Context, Result, ensure};
use clap::Parser;
use k2rmini::locate::{LocateOptions, locate_graph_mums, locate_shared_mums, read_state_info};

#[derive(Parser)]
#[command(
    version,
    about = "Build shared proto-MUMs inside GGCAT, then locate and verify them with K2Rmini"
)]
struct Args {
    /// One FASTA/FASTQ file per genome; all records in a file belong to that genome
    #[arg(required_unless_present = "input_list")]
    genomes: Vec<PathBuf>,
    /// UTF-8 file with one genome path per line; relative paths use the list's directory
    #[arg(long)]
    input_list: Option<PathBuf>,
    /// Minimum MUM length (primary knob). Sets discovery seed k = l unless --seed-k
    /// overrides. Bounds: 3..=64. Completeness is only down to this length.
    #[arg(short = 'l', long = "min-match-len", required_unless_present_any = ["update", "k"])]
    min_match_len: Option<usize>,
    /// Advanced: seed / dictionary k-mer length. Default: same as -l.
    /// Must satisfy k <= l. Smaller k than -l only discovers then discards short
    /// survivors (more RAM/CPU, never recovers MUMs shorter than k). Prefer -l alone.
    #[arg(short, long = "seed-k", visible_alias = "k")]
    k: Option<usize>,
    /// Minimum genomes sharing a MUM (default: all). Negative = relative to total
    /// (e.g. -N -1 means present in all but one), like mumemto -k.
    #[arg(short = 'N', long = "minimum-genomes", allow_hyphen_values = true)]
    minimum_genomes: Option<i64>,
    /// Write updateable merge metadata (packed genomes, protos, evidence) here.
    /// A later --update against this directory adds assemblies as if re-run on the union.
    /// Mumento2-only; not compatible with classic mumemto merge files.
    #[arg(long, alias = "mergeable")]
    state: Option<PathBuf>,
    /// Add NEW genome FASTAs to an existing --state/--mergeable directory (in place,
    /// or fork with --state DIR while updating from --update OTHER).
    #[arg(long)]
    update: Option<PathBuf>,
    /// Graph seed inputs (default: up to 16 for stateless all-genome searches, otherwise all)
    #[arg(long)]
    seed_genomes: Option<usize>,
    /// MUM output path; also writes <prefix>.lengths. With -b, write .bumbl.
    #[arg(short, long)]
    output: PathBuf,
    /// Write mumemto-compatible .bumbl (global coords). Incompatible with --contig-coords.
    #[arg(short = 'b', long = "binary")]
    binary: bool,
    /// Contig-relative 4-col coords (default: mumemto global length/starts/strands)
    #[arg(long)]
    contig_coords: bool,
    #[arg(short = 'j', long)]
    threads: Option<usize>,
    /// GGCAT's suggested scratch-memory budget in GiB
    #[arg(short = 'm', long, default_value_t = 2.0)]
    memory: f64,
    /// Modified GGCAT executable (default: the local build)
    #[arg(long)]
    ggcat: Option<PathBuf>,
    #[arg(long)]
    temp_dir: Option<PathBuf>,
    /// Keep intermediate proto-MUMs and GGCAT scratch files
    #[arg(long)]
    keep_temp: bool,
    #[arg(long)]
    forward_only: bool,
    #[arg(long)]
    ggcat_minimizer_length: Option<usize>,
    /// Legacy locator minimizer length, retained for checkpoint compatibility
    #[arg(long)]
    locate_minimizer_length: Option<usize>,
    #[arg(long)]
    buckets_count_log: Option<usize>,
    #[arg(long, default_value_t = 65_536)]
    chunk_bases: usize,
    /// Write phase wall times as TSV (also printed to stderr)
    #[arg(long)]
    timings: Option<PathBuf>,
    /// Save the shared proto-MUM FASTA for validation without keeping GGCAT scratch
    #[arg(long)]
    proto_output: Option<PathBuf>,
}

fn local_ggcat() -> PathBuf {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    for profile in ["mumento", "release", "debug"] {
        let path = root.join("ggcat/target").join(profile).join("ggcat");
        if path.is_file() {
            return path;
        }
    }
    PathBuf::from("ggcat")
}

/// Resolve `-N`: positive absolute count; negative means present in at least G-|N|
/// genomes (mumemto-style relative threshold). Zero remains invalid.
fn resolve_minimum_genomes(raw: i64, total: usize) -> Result<usize> {
    ensure!(total > 0, "at least one genome is required");
    ensure!(raw != 0, "minimum genome count must be non-zero");
    if raw > 0 {
        let n = raw as usize;
        ensure!(
            (1..=total).contains(&n),
            "require 1 <= N <= total genome count ({total})"
        );
        return Ok(n);
    }
    let omit = (-raw) as usize;
    ensure!(
        omit < total,
        "relative -N must satisfy |N| < genome count ({total})"
    );
    Ok(total - omit)
}

// A sequence present in N of G genomes must occur in at least N-(G-S) of
// any S seed genomes. This preserves every qualifying k-mer, including repeats.
// Reusable union states must retain all new genomes, even below the threshold.
fn seed_plan(
    genomes: usize,
    minimum: Option<usize>,
    reusable: bool,
    requested: Option<usize>,
) -> Result<(usize, usize)> {
    // Every strict all-genome match occurs in the first input. A complete
    // reference graph therefore contains every qualifying seed. A bounded
    // prefix further rejects private reference variation without building a
    // graph over the full collection. Exact validation still scans all inputs;
    // subset thresholds and reusable unions need broader candidate coverage.
    let seeds = requested.unwrap_or_else(|| {
        if !reusable && minimum.unwrap_or(genomes) == genomes {
            genomes.min(16)
        } else {
            genomes
        }
    });
    ensure!(
        (1..=genomes).contains(&seeds),
        "seed-genomes must be between 1 and the number of new input genomes"
    );
    if reusable {
        ensure!(
            seeds == genomes,
            "reusable states require union seeds from ALL new input genomes"
        );
        return Ok((seeds, 1));
    }
    let minimum = minimum.unwrap_or(genomes);
    ensure!(
        (1..=genomes).contains(&minimum),
        "minimum genome count is out of range"
    );
    let required=minimum.checked_sub(genomes-seeds).filter(|&n|n>0)
        .context("seed subset is too small: require seed-genomes >= G-N+1 to cover every qualifying support set")?;
    Ok((seeds, required))
}

fn expand_input_list(args: &mut Args) -> Result<()> {
    if let Some(path) = &args.input_list {
        let parent = path.parent().unwrap_or(Path::new("."));
        let text = std::fs::read_to_string(path).context("read genome input list")?;
        args.genomes.extend(
            text.lines()
                .filter(|line| !line.is_empty())
                .map(|line| parent.join(line)),
        );
    }
    ensure!(!args.genomes.is_empty(), "at least one genome is required");
    Ok(())
}

fn resolve_binary_output(path: PathBuf) -> PathBuf {
    let text = path.to_string_lossy();
    let rewritten = text
        .strip_suffix(".mums.zst")
        .or_else(|| text.strip_suffix(".tsv.zst"))
        .or_else(|| text.strip_suffix(".mums"))
        .or_else(|| text.strip_suffix(".tsv"))
        .map(|base| format!("{base}.bumbl"));
    rewritten.map(PathBuf::from).unwrap_or(path)
}

fn main() -> Result<()> {
    let total_started = Instant::now();
    let mut args = Args::parse();
    expand_input_list(&mut args)?;
    ensure!(
        !(args.binary && args.contig_coords),
        "binary bumbl output (-b) is incompatible with --contig-coords"
    );
    if args.binary {
        args.output = resolve_binary_output(args.output);
    }
    let saved = args
        .update
        .as_deref()
        .map(read_state_info)
        .transpose()
        .context("read update state")?;
    // Primary knob is -l (reported floor = discovery seed unless --seed-k/--k).
    // Backward compatible: -k alone still means k = l = that value.
    let min_length = args
        .min_match_len
        .or(args.k)
        .or_else(|| saved.as_ref().map(|s| s.min_length))
        .context("require -l / --min-match-len (or legacy -k) unless --update")?;
    let k = args
        .k
        .or_else(|| saved.as_ref().map(|s| s.k))
        .unwrap_or(min_length);
    ensure!(min_length > 0, "min-match-len (-l) must be positive");
    ensure!(k >= 3, "seed k must be at least 3");
    ensure!(min_length >= k, "min-match-len (-l) must be >= seed k");
    let total_genomes = args.genomes.len() + saved.as_ref().map_or(0, |s| s.genomes);
    let minimum_genomes = match args.minimum_genomes {
        Some(raw) => Some(resolve_minimum_genomes(raw, total_genomes)?),
        None => saved.as_ref().and_then(|s| s.minimum_genomes),
    };
    if let Some(n) = minimum_genomes {
        ensure!(
            (1..=total_genomes).contains(&n),
            "require 1 <= N <= total genome count ({total_genomes})"
        );
    }
    ensure!(
        args.memory.is_finite() && args.memory > 0.0,
        "memory must be finite and positive"
    );
    let threads = args
        .threads
        .unwrap_or_else(|| thread::available_parallelism().map_or(1, usize::from));
    let locate_options = LocateOptions {
        k,
        min_length,
        m: args
            .locate_minimizer_length
            .or_else(|| saved.as_ref().map(|s| s.m))
            .unwrap_or(k.min(21)),
        threads,
        chunk_bases: args.chunk_bases,
        forward_only: args.forward_only || saved.as_ref().is_some_and(|s| s.forward_only),
        mumemto_format: !args.contig_coords,
        binary: args.binary,
    };
    locate_options.validate()?;
    if let Some(info) = &saved {
        ensure!(
            (info.k, info.min_length, info.m, info.forward_only)
                == (
                    k,
                    locate_options.min_length,
                    locate_options.m,
                    locate_options.forward_only
                ),
            "updates must retain the state's k, min-match-len, minimizer length and strand mode"
        );
    }
    let mut seen = std::collections::HashSet::new();
    for genome in &args.genomes {
        ensure!(
            genome.is_file(),
            "genome does not exist: {}",
            genome.display()
        );
        ensure!(
            seen.insert(std::fs::canonicalize(genome)?),
            "the same genome file was supplied more than once"
        );
    }
    if args.update.is_none() {
        ensure!(
            !args
                .state
                .as_ref()
                .is_some_and(|p| p.join("CURRENT").exists()),
            "state already exists; use --update to add genomes"
        );
    }
    let reusable = args.state.is_some() || args.update.is_some();
    let (seed_genomes, seed_presence) = seed_plan(
        args.genomes.len(),
        minimum_genomes,
        reusable,
        args.seed_genomes,
    )?;
    let scratch = if let Some(parent) = &args.temp_dir {
        std::fs::create_dir_all(parent)?;
        tempfile::Builder::new()
            .prefix("mumento2-")
            .tempdir_in(parent)?
    } else {
        tempfile::Builder::new().prefix("mumento2-").tempdir()?
    };
    let scratch_path = if args.keep_temp {
        scratch.keep()
    } else {
        scratch.path().to_path_buf()
    };
    let protos = scratch_path.join("proto-mums.fasta.zst");
    let executable = args.ggcat.unwrap_or_else(local_ggcat);
    let mut build = Command::new(&executable);
    build
        .arg("build")
        .arg("-k")
        .arg(k.to_string())
        .arg("-j")
        .arg(threads.to_string())
        .arg("-m")
        .arg(args.memory.to_string())
        .arg("-t")
        .arg(scratch_path.join("ggcat"))
        .arg("-o")
        .arg(&protos);
    // A reusable state must retain below-threshold k-mers. Future genomes can
    // make them qualify. Updates only build the union graph of the NEW inputs.
    if reusable {
        build.arg("-s").arg("1");
    } else if seed_presence == seed_genomes {
        build.arg("--present-in-all-colors");
    } else {
        build
            .arg("--present-in-at-least-colors")
            .arg(seed_presence.to_string());
    }
    if locate_options.forward_only {
        build.arg("-f");
    }
    if args.keep_temp {
        build.arg("--keep-temp-files");
    }
    if let Some(m) = args.ggcat_minimizer_length {
        build.arg("--minimizer-length").arg(m.to_string());
    }
    if let Some(b) = args.buckets_count_log {
        build.arg("-b").arg(b.to_string());
    }
    if args.input_list.is_some() {
        // Keep the child command short too, including union/state builds that
        // need every input genome rather than a bounded seed prefix.
        let list = scratch_path.join("graph-inputs.list");
        let cwd = std::env::current_dir()?;
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&list)?);
        for genome in &args.genomes[..seed_genomes] {
            let path = cwd.join(genome);
            let text = path.to_str().context("input-list paths must be UTF-8")?;
            ensure!(
                !text.contains(['\n', '\r']),
                "input-list paths cannot contain line breaks"
            );
            writeln!(writer, "{text}")?;
        }
        writer.flush()?;
        build.arg("--input-lists").arg(list);
    } else {
        build.args(&args.genomes[..seed_genomes]);
    }
    eprintln!(
        "Building proto-MUMs with native GGCAT ({} inputs; {} seeds)...",
        seed_genomes,
        if reusable {
            "union"
        } else {
            "genome-presence filtered"
        }
    );
    if seed_genomes < args.genomes.len() {
        eprintln!(
            "Seed subset: {seed_genomes}/{} genomes, presence >= {seed_presence}; full validation still uses every genome",
            args.genomes.len()
        );
    }
    let started = Instant::now();
    let status = build
        .status()
        .with_context(|| format!("run modified GGCAT at {}", executable.display()))?;
    ensure!(status.success(), "native GGCAT build failed ({status})");
    let graph_seconds = started.elapsed().as_secs_f64();
    if let Some(path) = &args.proto_output {
        std::fs::copy(&protos, path).context("save proto-MUM FASTA")?;
    }
    eprintln!("Locating and verifying proto-MUMs with K2Rmini...");
    let stats = if reusable || minimum_genomes.is_some_and(|n| n < total_genomes) {
        locate_shared_mums(
            &protos,
            &args.genomes,
            &args.output,
            Some(&scratch_path),
            locate_options,
            minimum_genomes,
            args.state.as_deref(),
            args.update.as_deref(),
        )?
    } else {
        locate_graph_mums(
            &protos,
            &args.genomes,
            &args.output,
            Some(&scratch_path),
            locate_options,
        )?
    };
    eprintln!(
        "{} MUMs written to {} ({} duplicate outputs merged)",
        stats.verified_mums,
        args.output.display(),
        stats.merged_duplicates
    );
    if args.keep_temp {
        eprintln!("Intermediate files: {}", scratch_path.display());
    }
    let mut phases = vec![("ggcat_build", graph_seconds)];
    phases.extend(stats.timings.phases());
    phases.push(("pipeline_total", total_started.elapsed().as_secs_f64()));
    for (phase, seconds) in &phases {
        eprintln!("[timing]\t{phase}\t{seconds:.6}");
    }
    if let Some(path) = args.timings {
        let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(writer, "phase\twall_seconds")?;
        for (phase, seconds) in phases {
            writeln!(writer, "{phase}\t{seconds:.6}")?;
        }
        writer.flush()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn automatic_graph_seed_plan_preserves_threshold_and_reusable_coverage() {
        for g in 1..=100 {
            let seeds = g.min(16);
            assert_eq!(seed_plan(g, None, false, None).unwrap(), (seeds, seeds));
            assert_eq!(seed_plan(g, Some(g), false, None).unwrap(), (seeds, seeds));
            assert_eq!(seed_plan(g, None, false, Some(g)).unwrap(), (g, g));
            assert_eq!(seed_plan(g, None, true, None).unwrap(), (g, 1));
            for n in 1..g {
                assert_eq!(seed_plan(g, Some(n), false, None).unwrap(), (g, n));
                assert_eq!(seed_plan(g, Some(n), true, None).unwrap(), (g, 1));
            }
        }
    }
    #[test]
    fn subset_presence_bound_covers_every_support_set() {
        for g in 1usize..=10 {
            for n in 1..=g {
                for s in g - n + 1..=g {
                    let (_, q) = seed_plan(g, Some(n), false, Some(s)).unwrap();
                    for support in 0usize..1 << g {
                        if support.count_ones() as usize >= n {
                            assert!((support & ((1 << s) - 1)).count_ones() as usize >= q);
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn unsafe_support_subsets_are_rejected() {
        for g in 2..=10 {
            for n in 1..g {
                for s in 1..=g - n {
                    assert!(seed_plan(g, Some(n), false, Some(s)).is_err());
                }
            }
        }
    }
    #[test]
    fn all_genome_matches_allow_any_nonempty_seed_subset() {
        for s in 1..=100 {
            assert_eq!(seed_plan(100, None, false, Some(s)).unwrap(), (s, s));
        }
    }
    #[test]
    fn reusable_states_cannot_drop_below_threshold_candidates() {
        assert!(seed_plan(10, Some(2), true, Some(9)).is_err());
        assert_eq!(seed_plan(10, Some(2), true, None).unwrap(), (10, 1));
    }
    #[test]
    fn zero_and_oversized_subsets_are_rejected() {
        for s in [0, 11] {
            assert!(seed_plan(10, None, false, Some(s)).is_err());
        }
    }

    #[test]
    fn relative_minimum_genomes_matches_mumemto_style() {
        assert_eq!(resolve_minimum_genomes(3, 5).unwrap(), 3);
        assert_eq!(resolve_minimum_genomes(-1, 5).unwrap(), 4);
        assert_eq!(resolve_minimum_genomes(-2, 5).unwrap(), 3);
        assert!(resolve_minimum_genomes(0, 5).is_err());
        assert!(resolve_minimum_genomes(-5, 5).is_err());
        assert!(resolve_minimum_genomes(6, 5).is_err());
    }

    #[test]
    fn primary_l_or_legacy_k_is_accepted() {
        let with_l = Args::try_parse_from(["mumento2", "-l", "40", "-o", "out", "g.fa"]).unwrap();
        assert_eq!(with_l.min_match_len, Some(40));
        assert_eq!(with_l.k, None);
        let with_k = Args::try_parse_from(["mumento2", "-k", "21", "-o", "out", "g.fa"]).unwrap();
        assert_eq!(with_k.k, Some(21));
        assert_eq!(with_k.min_match_len, None);
        let both = Args::try_parse_from([
            "mumento2", "-l", "40", "-k", "21", "-o", "out", "g.fa",
        ])
        .unwrap();
        assert_eq!((both.min_match_len, both.k), (Some(40), Some(21)));
    }

    #[test]
    fn binary_output_rewrites_mums_suffix() {
        assert_eq!(
            resolve_binary_output(PathBuf::from("out.mums.zst")),
            PathBuf::from("out.bumbl")
        );
        assert_eq!(
            resolve_binary_output(PathBuf::from("out.tsv")),
            PathBuf::from("out.bumbl")
        );
        assert_eq!(
            resolve_binary_output(PathBuf::from("custom.bumbl")),
            PathBuf::from("custom.bumbl")
        );
    }

    #[test]
    fn input_list_preserves_order_spaces_crlf_and_relative_paths() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let list = directory.path().join("inputs.list");
        std::fs::write(
            &list,
            "with spaces.fa\r\n\r\nrelative/second.fa\r\n/absolute/third.fa\n",
        )?;
        let mut args = Args::try_parse_from([
            "mumento2",
            "-k",
            "31",
            "-o",
            "out",
            "--input-list",
            list.to_str().unwrap(),
            "first.fa",
        ])?;
        expand_input_list(&mut args)?;
        assert_eq!(
            args.genomes,
            vec![
                PathBuf::from("first.fa"),
                directory.path().join("with spaces.fa"),
                directory.path().join("relative/second.fa"),
                PathBuf::from("/absolute/third.fa")
            ]
        );
        Ok(())
    }

    #[test]
    fn input_list_accepts_100000_entries_with_a_short_command() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let list = directory.path().join("inputs.list");
        let paths: String = (0..100000).map(|i| format!("genome-{i}.fa\n")).collect();
        std::fs::write(&list, paths)?;
        let mut args = Args::try_parse_from([
            "mumento2",
            "-k",
            "31",
            "-o",
            "out",
            "--input-list",
            list.to_str().unwrap(),
        ])?;
        expand_input_list(&mut args)?;
        assert_eq!(args.genomes.len(), 100000);
        assert_eq!(
            args.genomes[99999],
            directory.path().join("genome-99999.fa")
        );
        std::fs::write(&list, "\n\r\n")?;
        let mut args = Args::try_parse_from([
            "mumento2",
            "-k",
            "31",
            "-o",
            "out",
            "--input-list",
            list.to_str().unwrap(),
        ])?;
        assert!(expand_input_list(&mut args).is_err());
        Ok(())
    }
}
