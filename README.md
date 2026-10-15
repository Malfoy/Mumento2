# Mumento2

**Work in progress.** A Rust implementation inspired by [Mumemto](https://github.com/vikshiv/mumemto), using modified versions of [GGCAT](https://github.com/algbio/ggcat) and [K2Rmini](https://github.com/Malfoy/K2Rmini) to find maximal unique matches shared by multiple genomes.

## Build

Requires a recent stable Rust toolchain, a C/C++ compiler, CMake, pkg-config, and zlib/xz development libraries. Tested on Linux.

```sh
git clone https://github.com/Malfoy/Mumento2.git
cd Mumento2
cargo build --locked --manifest-path ggcat/Cargo.toml -p ggcat_cmdline --profile mumento
cargo build --locked --release
```

## Run

```sh
target/release/mumento2 -k 31 -j 32 -m 1 -o matches.tsv.zst genome1.fa genome2.fa genome3.fa
```

Each input file is one genome; its records are separate contigs. FASTA/FASTQ inputs may be plain, gzip, xz, or zstd compressed. Both strands are searched.

`-k` sets the k-mer size and minimum match length, `-j` sets the worker count, and `-m` sets GGCAT's suggested memory budget in GiB. Output is a zstd-compressed TSV with match lengths, zero-based starts, strands, and contig IDs. Add `--timings timings.tsv` to save phase timings.
