# Mumento2

**Work in progress.** A Rust implementation inspired by [Mumemto](https://github.com/vikshiv/mumemto), using modified [GGCAT](https://github.com/algbio/ggcat) and [K2Rmini](https://github.com/Malfoy/K2Rmini) components to discover and validate maximal unique matches across genome collections.

Mumento2 searches both strands and respects contig boundaries. It reports maximal exact matches that occur once in every participating genome. By default a match must occur in all input genomes; a configurable minimum genome count supports matches shared by a subset, including matches absent from the first genome. Independently valid internal fragments survive when a longer proto-MUM fails validation.

Reusable state preserves union candidates, packed genomes and exact position-validity evidence. Adding genomes revalidates those candidates, allowing matches to shorten or split and newly qualifying matches to appear. State generations are immutable and committed atomically. Collinear filtering and a full-text query index are not currently provided.

The streaming all-genome path uses parallel exact seed scanning and extension, batched scratch I/O, inline alignment caches and progressive pruning of individual proto seed positions. A seed position is removed only when exact evidence proves it absent from an already processed genome. Reusable state retains raw candidates needed for later updates.

Native CPU optimization, optimization level 3, fat LTO and a single code-generation unit are enabled for optimized builds. A validated human LLVM 22 compilation profile is included. Training used 32 chromosome 19 haplotypes and one thread; the profile build helper checks exact source/compiler fingerprints and the profile digest. Performance outside the measured workload requires separate measurements.

Inputs may be plain or compressed FASTA/FASTQ, with each file representing one genome and each record a separate contig. Large collections can use a file list. Output is a compressed TSV containing match lengths, zero-based forward-strand starts, strands and contig IDs, plus contig metadata. Phase timings expose graph construction, indexing, scanning/extension, reduction and output costs.

The test suite contains **231 library tests, 8 driver tests and 205 end-to-end tests**, including independent small-input oracles, internal-fragment recovery, reverse strands, repeats, thresholds and incremental updates. See [test coverage](tests/README.md).

Human measurements use real HPRC **chromosome 19 haplotypes**, following the chromosome scaling experiment in the [Mumemto paper](https://doi.org/10.1186/s13059-025-03644-0). Inputs average approximately 61 Mb; uniqueness is checked within chromosome 19. The 512-input collection represents 276 individuals, preserving the Release 2 prefix and adding public HPRC R3 assemblies. Whole-human-genome performance has not been measured. See the [input provenance and selection](docs/human-chr19-dataset-20261018.md).

The latest measured configuration uses the same frozen **human-PGO executable** at every size, k = 63, minimum MUM length = 63 bp, both strands and strict presence in every input. Rust uses 32 threads, up to eight graph-seed genomes and a process-local **mimalloc 3.2.8** preload. Each size has one measured Rust run; all saved C++ measurements and the completed 512-input Rust run are reused.

| Haplotypes | C++ wall (s) | Rust wall (s) | C++ CPU (s) | Rust CPU (s) | C++ RAM (GiB) | Rust RAM (GiB) |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 103.82 | 9.95 (10.43×) | 103.68 | 48.65 (2.13×) | 2.58 | 2.46 (1.05×) |
| 4 | 125.58 | 10.36 (12.12×) | 372.91 | 55.58 (6.71×) | 8.32 | 2.63 (3.16×) |
| 8 | 147.97 | 10.71 (13.81×) | 1,018.21 | 73.99 (13.76×) | 19.34 | 3.23 (5.98×) |
| 16 | 221.77 | 12.20 (18.17×) | 3,022.01 | 102.83 (29.39×) | 41.48 | 5.11 (8.12×) |
| 32 | 454.91 | 13.02 (34.93×) | 6,393.35 | 180.93 (35.34×) | 44.24 | 9.21 (4.80×) |
| 64 | 868.70 | 18.24 (47.62×) | 13,076.85 | 314.68 (41.56×) | 44.63 | 9.18 (4.86×) |
| 128 | 1,737.04 | 27.90 (62.26×) | 26,275.72 | 579.18 (45.37×) | 45.84 | 10.28 (4.46×) |
| 256 | 3,464.91 | 46.45 (74.59×) | 52,896.79 | 1,126.70 (46.95×) | 46.45 | 10.24 (4.53×) |
| 512 | 6,931.30 | 86.80 (79.86×) | 106,287.83 | 2,335.49 (45.51×) | 46.23 | 9.97 (4.64×) |

Parentheses give C++ divided by Rust: factors above one mean Rust uses less wall time, CPU time or RAM. These are single observations. The fastest observed 512-input configuration takes **86.80 seconds**, versus 6,931.30 seconds for C++: **79.86× faster**. Without the allocator preload, the same PGO executable and seed prefix took **87.53 seconds** and **7.22 GiB**, versus **9.97 GiB** with the preload. The small timing difference comes with a substantial RAM cost.

Every output at all nine sizes matches C++ exactly, including lengths, coordinates and strands; the 512-input output contains **254,892 MUMs**. Sequence, contig-boundary, joint maximality and physical uniqueness audits are inherited only through exact full row-set equality with independently validated C++ output.

For stateless all-genome searches, graph construction uses up to 16 seed genomes by default. Every qualifying MUM occurs in that prefix; validation still scans every input genome. Threshold searches and reusable states retain full candidate coverage. The native Rust SSHash proto-MUM index incorporates PTRHash/SIMD acceleration from REINDEER. Directional stop maps allow exact alignments to cross proto boundaries when every later oriented occurrence is already covered. Copies may disappear; new mappings remain stopping points. Unresolved contexts compare DNA only within these proven bounds. The compact scanner uses a bounded exact alignment cache; collisions cause recomputation, and duplicate evidence is removed before reduction. SSHash uses m = min(23, k−1). Other retained optimizations bound repeat evidence to three distinct oriented mappings per genome, share packed DNA buffers between the scanner and SSHash, use 32-byte repeat-context nodes and exact alignment caches, preserve dictionary-only bucket information at query jumps and skip absent-minimizer windows.

Rust requests 32 native threads. Saved C++ runs used 1, 3, 7 and 15 partition workers at 2, 4, 8 and 16 inputs, respectively; runs from 32 inputs used 16 workers, including merging, constrained by the 48 GiB memory cap. These measurements used a 32-logical-CPU Intel i9-13950HX. CPU is summed user plus system time, including children. Timing includes graph construction, indexing, scanning, reduction and output; downloads, builds, staging and correctness audits are excluded. Single observations are not statistical performance estimates.

See the [complete comparison](docs/human-global-all-sizes-latest-20261018.json), [CSV results](docs/human-global-all-sizes-latest-20261018.csv), [tab-separated email table](docs/human-global-all-sizes-latest-20261018.tsv), [512-input phase timings](docs/human-global-512-latest-20261018.json) and [optimization record](docs/human-optimization-512-20261018.json). All outputs at all nine sizes match the complete independently validated C++ row sets.
