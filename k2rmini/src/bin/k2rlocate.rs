use clap::Parser;
use k2rmini::locate::{LocateOptions, locate_mums};
use std::{io, path::PathBuf, thread};

#[derive(Parser)]
#[command(about = "Locate GGCAT proto-MUMs using K2Rmini and verify maximal unique matches")]
struct Args {
    /// Original genome files, in output coordinate order
    #[arg(required = true)]
    genomes: Vec<PathBuf>,
    /// GGCAT proto-MUM FASTA
    #[arg(short)]
    patterns: PathBuf,
    /// Zstd-compressed MUM output; also writes <prefix>.lengths
    #[arg(short)]
    output: PathBuf,
    /// Contig-relative 4-col coords (default: mumemto global length/starts/strands)
    #[arg(long)]
    contig_coords: bool,
    /// Seed / dictionary k-mer length
    #[arg(short, default_value_t = 31)]
    k: usize,
    /// Minimum reported MUM length (default: k)
    #[arg(short = 'l', long = "min-match-len")]
    min_match_len: Option<usize>,
    /// Minimizer length, defaults to min(k, 21)
    #[arg(short)]
    m: Option<usize>,
    #[arg(short = 'T', long)]
    threads: Option<usize>,
    #[arg(long, default_value_t = 65_536)]
    chunk_bases: usize,
    #[arg(long)]
    temp_dir: Option<PathBuf>,
    #[arg(long)]
    forward_only: bool,
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    let options = LocateOptions {
        k: args.k,
        min_length: args.min_match_len.unwrap_or(args.k),
        m: args.m.unwrap_or(args.k.min(21)),
        threads: args
            .threads
            .unwrap_or_else(|| thread::available_parallelism().map_or(1, usize::from)),
        chunk_bases: args.chunk_bases,
        forward_only: args.forward_only,
        mumemto_format: !args.contig_coords,
        binary: false,
    };
    let stats = locate_mums(
        &args.patterns,
        &args.genomes,
        &args.output,
        args.temp_dir.as_deref(),
        options,
    )?;
    eprintln!(
        "K2Rmini: {} verified MUMs, {} duplicate outputs merged",
        stats.verified_mums, stats.merged_duplicates
    );
    Ok(())
}
