use std::{
    io::Write,
    path::{Path, PathBuf},
    process::Command,
    thread,
    time::Instant,
};

use anyhow::{Context, Result, ensure};
use clap::Parser;
use k2rmini::locate::{LocateOptions, locate_mums};

#[derive(Parser)]
#[command(
    version,
    about = "Build shared proto-MUMs inside GGCAT, then locate and verify them with K2Rmini"
)]
struct Args {
    /// One FASTA/FASTQ file per genome; all records in a file belong to that genome
    #[arg(required = true)]
    genomes: Vec<PathBuf>,
    /// K-mer length and minimum MUM length
    #[arg(short, long)]
    k: usize,
    /// Zstd-compressed MUM TSV, plus <output>.contigs.tsv.zst
    #[arg(short, long)]
    output: PathBuf,
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
    /// K2Rmini minimizer length, defaults to min(k, 21)
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

fn main() -> Result<()> {
    let total_started = Instant::now();
    let args = Args::parse();
    ensure!(args.k >= 3, "k must be at least 3");
    ensure!(
        args.memory.is_finite() && args.memory > 0.0,
        "memory must be finite and positive"
    );
    let threads = args
        .threads
        .unwrap_or_else(|| thread::available_parallelism().map_or(1, usize::from));
    let locate_options = LocateOptions {
        k: args.k,
        m: args.locate_minimizer_length.unwrap_or(args.k.min(21)),
        threads,
        chunk_bases: args.chunk_bases,
        forward_only: args.forward_only,
    };
    locate_options.validate()?;
    for genome in &args.genomes {
        ensure!(
            genome.is_file(),
            "genome does not exist: {}",
            genome.display()
        );
    }
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
        .arg("--present-in-all-colors")
        .arg("-k")
        .arg(args.k.to_string())
        .arg("-j")
        .arg(threads.to_string())
        .arg("-m")
        .arg(args.memory.to_string())
        .arg("-t")
        .arg(scratch_path.join("ggcat"))
        .arg("-o")
        .arg(&protos);
    if args.forward_only {
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
    build.args(&args.genomes);
    eprintln!("Building proto-MUMs with native GGCAT all-genome presence filtering...");
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
    let stats = locate_mums(
        &protos,
        &args.genomes,
        &args.output,
        Some(&scratch_path),
        locate_options,
    )?;
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
