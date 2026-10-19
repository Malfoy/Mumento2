# Mumento2

**Work in progress.** A Rust implementation inspired by [Mumemto](https://github.com/vikshiv/mumemto), using modified [GGCAT](https://github.com/algbio/ggcat) and [K2Rmini](https://github.com/Malfoy/K2Rmini) components to discover and validate maximal unique matches across genome collections.

Mumento2 searches both strands and respects contig boundaries. It reports maximal exact matches that occur once in every participating genome. By default a match must occur in all input genomes; a configurable minimum genome count supports matches shared by a subset, including matches absent from the first genome. Independently valid internal fragments survive when a longer proto-MUM fails validation.

Reusable state preserves union candidates, packed genomes and exact position-validity evidence. Adding genomes revalidates those candidates, allowing matches to shorten or split and newly qualifying matches to appear. State generations are immutable and committed atomically. Collinear filtering and a full-text query index are not currently provided.

The streaming all-genome path uses parallel exact seed scanning and extension, batched scratch I/O, inline alignment caches and progressive pruning of individual proto seed positions. A seed position is removed only when exact evidence proves it absent from an already processed genome. Reusable state retains raw candidates needed for later updates.

Native CPU optimization, optimization level 3, fat LTO and a single code-generation unit are enabled for optimized builds. A validated human LLVM 22 compilation profile is included. Training used 32 chromosome 19 haplotypes and one thread; the profile build helper checks exact source/compiler fingerprints and the profile digest. Performance outside the measured workload requires separate measurements.

Inputs may be plain or compressed FASTA/FASTQ, with each file representing one genome and each record a separate contig. Large collections can use a file list. Output is a compressed TSV containing match lengths, zero-based forward-strand starts, strands and contig IDs, plus contig metadata. Phase timings expose graph construction, indexing, scanning/extension, reduction and output costs.

The test suite contains **237 library tests, 8 driver tests and 205 end-to-end tests**, including independent small-input oracles, internal-fragment recovery, reverse strands, repeats, thresholds and incremental updates. See [test coverage](tests/README.md).

Human measurements use real HPRC **chromosome 19 haplotypes**, following the chromosome scaling experiment in the [Mumemto paper](https://doi.org/10.1186/s13059-025-03644-0). Inputs average approximately 61 Mb; uniqueness is checked within chromosome 19. The 512-input collection represents 276 individuals, preserving the Release 2 prefix and adding public HPRC R3 assemblies. Whole-human-genome performance has not been measured. See the [input provenance and selection](docs/human-chr19-dataset-20261018.md).

The latest measured configuration uses the same frozen **human-PGO executable** at every size, k = 63, minimum MUM length = 63 bp, both strands and strict presence in every input. Rust uses 32 threads, up to 16 graph-seed genomes and the system allocator, with no allocator preload. Each size has one measured Rust run; all saved C++ measurements and the completed 512-input Rust run are reused.

GGCAT uses graph minimizer length 14 and 256 buckets. The SSHash dictionary uses minimizer length 31 at k = 63.

| Haplotypes | C++ wall (s) | Rust wall (s) | C++ CPU (s) | Rust CPU (s) | C++ RAM (GiB) | Rust RAM (GiB) |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 103.82 | 9.92 (10.46×) | 103.68 | 46.37 (2.24×) | 2.58 | 1.69 (1.53×) |
| 4 | 125.58 | 9.59 (13.09×) | 372.91 | 50.92 (7.32×) | 8.32 | 1.98 (4.21×) |
| 8 | 147.97 | 10.15 (14.57×) | 1,018.21 | 66.87 (15.23×) | 19.34 | 2.59 (7.47×) |
| 16 | 221.77 | 11.75 (18.87×) | 3,022.01 | 107.84 (28.02×) | 41.48 | 3.90 (10.63×) |
| 32 | 454.91 | 12.48 (36.45×) | 6,393.35 | 158.93 (40.23×) | 44.24 | 6.62 (6.69×) |
| 64 | 868.70 | 16.36 (53.12×) | 13,076.85 | 254.79 (51.32×) | 44.63 | 6.90 (6.47×) |
| 128 | 1,737.04 | 23.29 (74.58×) | 26,275.72 | 457.99 (57.37×) | 45.84 | 7.13 (6.43×) |
| 256 | 3,464.91 | 37.67 (91.99×) | 52,896.79 | 872.07 (60.66×) | 46.45 | 7.23 (6.43×) |
| 512 | 6,931.30 | 67.21 (103.13×) | 106,287.83 | 1,745.12 (60.91×) | 46.23 | 7.33 (6.30×) |

Parentheses give C++ divided by Rust: factors above one mean Rust uses less wall time, CPU time or RAM. These are single observations. The fastest observed 512-input configuration takes **67.21 seconds**, versus 6,931.30 seconds for C++: **103.13× faster**.

Every output at all nine sizes matches C++ exactly, including lengths, coordinates and strands; the 512-input output contains **254,892 MUMs**. Sequence, contig-boundary, joint maximality and physical uniqueness audits are inherited only through exact full row-set equality with independently validated C++ output.

For stateless all-genome searches, graph construction uses up to 16 seed genomes by default. Every qualifying MUM occurs in that prefix; validation still scans every input genome. Threshold searches and reusable states retain full candidate coverage. The native Rust SSHash proto-MUM index incorporates PTRHash/SIMD acceleration from REINDEER. Directional stop maps allow exact alignments to cross proto boundaries when every later oriented occurrence is already covered. Copies may disappear; new mappings remain stopping points. Unresolved contexts compare DNA only within these proven bounds. The compact scanner uses a bounded exact alignment cache; collisions cause recomputation, and duplicate evidence is removed before reduction. SSHash uses m = 31 at k = 63 and m = min(23, k−1) otherwise. SIMD minimizer batches use 16 KiB to improve per-worker cache locality; exact point-query differential tests cover batch boundaries, jumps, resets and absent-window skips. Already verified alignments resolve repeated-seed context leaves through exact occurrence keys and complete context-depth coverage; every failed proof retains the original DNA-checked tree lookup. Only occurrences represented by actual context leaves are eligible. A bounded bitmap rank directory finds their original records without binary searches, while checking the full position key and group membership. Sparse or large coordinate spaces retain the original search. Other retained optimizations bound repeat evidence to three distinct oriented mappings per genome, share packed DNA buffers between the scanner and SSHash, use 32-byte repeat-context nodes and exact alignment caches, preserve dictionary-only bucket information at query jumps and skip absent-minimizer windows.

Rust requests 32 native threads. Saved C++ runs used 1, 3, 7 and 15 partition workers at 2, 4, 8 and 16 inputs, respectively; runs from 32 inputs used 16 workers, including merging, constrained by the 48 GiB memory cap. These measurements used a 32-logical-CPU Intel i9-13950HX. CPU is summed user plus system time, including children. Timing includes graph construction, indexing, scanning, reduction and output; downloads, builds, staging and correctness audits are excluded. Single observations are not statistical performance estimates.

See the [complete comparison](docs/human-global-all-sizes-latest-20261018.json), [CSV results](docs/human-global-all-sizes-latest-20261018.csv), [tab-separated email table](docs/human-global-all-sizes-latest-20261018.tsv), [512-input phase timings](docs/human-global-512-latest-20261018.json) and [optimization record](docs/human-optimization-512-20261018.json). All outputs at all nine sizes match the complete independently validated C++ row sets.
