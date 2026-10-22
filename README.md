# Mumento2

**Work in progress.** A Rust implementation inspired by [Mumemto](https://github.com/vikshiv/mumemto), using modified [GGCAT](https://github.com/algbio/ggcat) and [K2Rmini](https://github.com/Malfoy/K2Rmini) components to discover and validate maximal unique matches across genome collections.

Mumento2 searches both strands and respects contig boundaries. It reports maximal exact matches that occur once in every participating genome. By default a match must occur in all input genomes; a configurable minimum genome count supports matches shared by a subset, including matches absent from the first genome. Independently valid internal fragments survive when a longer proto-MUM fails validation.

Reusable state preserves union candidates, packed genomes and exact position-validity evidence. Adding genomes revalidates those candidates, allowing matches to shorten or split and newly qualifying matches to appear. State generations are immutable and committed atomically. Collinear filtering and a full-text query index are not currently provided.

The streaming all-genome path uses parallel exact seed scanning and extension, batched scratch I/O, inline alignment caches and progressive pruning of individual proto seed positions. A seed position is removed only when exact evidence proves it absent from an already processed genome. Reusable state retains raw candidates needed for later updates.

Native CPU optimization, optimization level 3, fat LTO and a single code-generation unit are enabled for optimized builds. A validated human LLVM 22 compilation profile is included. Training used 32 chromosome 19 haplotypes and one thread; the profile build helper checks exact source/compiler fingerprints and the profile digest. Performance outside the measured workload requires separate measurements.

Inputs may be plain or compressed FASTA/FASTQ, with each file representing one genome and each record a separate contig. Large collections can use a file list. Primary length knob is `-l` / `--min-match-len` (sets discovery seed `k = l`). Optional advanced `--seed-k` / `-k` may set a smaller seed than `-l` (more RAM/CPU; never recovers MUMs shorter than `k`). Legacy `-k` alone still means `k = l`. By default output is mumemto-shaped (zstd) `length`, global contig-concat starts, and strands with no header, plus a mumemto-format `.lengths` file. `-b` / `--binary` writes mumemto-compatible `.bumbl` (global coords only; incompatible with `--contig-coords`). `--contig-coords` selects the previous contig-relative 4-col TSV (still headerless) and also writes `.contigs.tsv.zst`. Partial matches use `-N` (negative values are relative to the genome count). Updateable runs use `--state`/`--mergeable` and `--update` (mumento2-only metadata; see [docs/HANDOFF.md](docs/HANDOFF.md)). Phase timings expose graph construction, indexing, scanning/extension, reduction and output costs.

The test suite contains **242 library tests, 8 driver tests and 206 end-to-end tests**, including independent small-input oracles, internal-fragment recovery, reverse strands, repeats, thresholds and incremental updates. See [test coverage](tests/README.md).

Human measurements use real HPRC **chromosome 19 haplotypes**, following the chromosome scaling experiment in the [Mumemto paper](https://doi.org/10.1186/s13059-025-03644-0). Inputs average approximately 61 Mb; uniqueness is checked within chromosome 19. The 512-input collection represents 276 individuals, preserving the Release 2 prefix and adding public HPRC R3 assemblies. Whole-human-genome performance has not been measured. See the [input provenance and selection](docs/human-chr19-dataset-20261018.md).

The benchmark table records the frozen v101 executable; the subsequent small-k SSHash correction has not been benchmarked.

The latest measured configuration uses the same frozen **human-PGO executable** at every size, k = 63, minimum MUM length = 63 bp, both strands and strict presence in every input. Rust uses 32 threads, up to 16 graph-seed genomes and the system allocator, with no allocator preload. Each size has one measured Rust run; all saved C++ measurements and the completed 512-input Rust run are reused.

GGCAT uses graph minimizer length 14 and 256 buckets. The SSHash dictionary uses minimizer length 31 at k = 63.

| Haplotypes | C++ wall (s) | Rust wall (s) | C++ CPU (s) | Rust CPU (s) | C++ RAM (GiB) | Rust RAM (GiB) |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 103.82 | 10.24 (10.13×) | 103.68 | 47.55 (2.18×) | 2.58 | 1.88 (1.37×) |
| 4 | 125.58 | 9.99 (12.57×) | 372.91 | 53.25 (7.00×) | 8.32 | 2.17 (3.83×) |
| 8 | 147.97 | 10.57 (14.00×) | 1,018.21 | 69.66 (14.62×) | 19.34 | 2.82 (6.86×) |
| 16 | 221.77 | 12.10 (18.33×) | 3,022.01 | 111.90 (27.01×) | 41.48 | 4.15 (9.99×) |
| 32 | 454.91 | 12.94 (35.14×) | 6,393.35 | 161.27 (39.64×) | 44.24 | 6.77 (6.53×) |
| 64 | 868.70 | 16.26 (53.44×) | 13,076.85 | 252.63 (51.76×) | 44.63 | 7.17 (6.22×) |
| 128 | 1,737.04 | 22.69 (76.57×) | 26,275.72 | 441.01 (59.58×) | 45.84 | 7.30 (6.28×) |
| 256 | 3,464.91 | 36.53 (94.86×) | 52,896.79 | 836.25 (63.25×) | 46.45 | 7.46 (6.23×) |
| 512 | 6,931.30 | 66.34 (104.48×) | 106,287.83 | 1,700.70 (62.50×) | 46.23 | 7.59 (6.09×) |

Parentheses give C++ divided by Rust: factors above one mean Rust uses less wall time, CPU time or RAM. These are single observations. The fastest observed 512-input configuration takes **66.34 seconds**, versus 6,931.30 seconds for C++: **104.48× faster**.

Every output at all nine sizes matches C++ exactly, including lengths, coordinates and strands; the 512-input output contains **254,892 MUMs**. Sequence, contig-boundary, joint maximality and physical uniqueness audits are inherited only through exact full row-set equality with independently validated C++ output.

Verified singleton contexts can cover additional repeated-seed boundaries inside an already exact reference alignment. Each covered boundary retains its original full-key singleton leaf proof, requires at most k+64 context bases, excludes conflicting physical mappings, and keeps 64 verified bases beyond every skipped seed. Ordinary boundary checks apply near alignment ends.

A bounded reference-coordinate seed directory can reuse the next proto seed when its entire k-mer is already proved by a verified DNA alignment. Conflicting or unsupported mappings retain SSHash lookup, and every occurrence boundary still passes the original context restriction. The optional directory is capped at 512 MiB.

For stateless all-genome searches, graph construction uses up to 16 seed genomes by default. Every qualifying MUM occurs in that prefix; validation still scans every input genome. Threshold searches and reusable states retain full candidate coverage. The native Rust SSHash proto-MUM index incorporates PTRHash/SIMD acceleration from REINDEER. Directional stop maps allow exact alignments to cross proto boundaries when every later oriented occurrence is already covered. Copies may disappear; new mappings remain stopping points. Unresolved contexts compare DNA only within these proven bounds. The compact scanner uses a bounded exact alignment cache; collisions cause recomputation, and duplicate evidence is removed before reduction. SSHash uses m = 31 at k = 63 and m = min(23, max(1, floor(k/2))) otherwise. Small-k, even-k, forward-only, threshold and incremental searches retain SSHash for 3 ≤ k ≤ 64; overlapping proto positions are represented by packed occurrence lists, with sorting used to construct a distinct-seed SPSS. SIMD minimizer batches use 16 KiB to improve per-worker cache locality; exact point-query differential tests cover batch boundaries, jumps, resets and absent-window skips. Already verified alignments resolve repeated-seed context leaves through exact occurrence keys and complete context-depth coverage; every failed proof retains the original DNA-checked tree lookup. Only occurrences represented by actual context leaves are eligible. A bounded bitmap rank directory finds their original records without binary searches, while checking the full position key and group membership. Sparse or large coordinate spaces retain the original search. Other retained optimizations bound repeat evidence to three distinct oriented mappings per genome, share packed DNA buffers between the scanner and SSHash, use 32-byte repeat-context nodes and exact alignment caches, preserve dictionary-only bucket information at query jumps and skip absent-minimizer windows.

Rust requests 32 native threads. Saved C++ runs used 1, 3, 7 and 15 partition workers at 2, 4, 8 and 16 inputs, respectively; runs from 32 inputs used 16 workers, including merging, constrained by the 48 GiB memory cap. These measurements used a 32-logical-CPU Intel i9-13950HX. CPU is summed user plus system time, including children. Timing includes graph construction, indexing, scanning, reduction and output; downloads, builds, staging and correctness audits are excluded. Single observations are not statistical performance estimates.

See the [complete comparison](docs/human-global-all-sizes-latest-20261018.json), [CSV results](docs/human-global-all-sizes-latest-20261018.csv), [tab-separated email table](docs/human-global-all-sizes-latest-20261018.tsv), [512-input phase timings](docs/human-global-512-latest-20261018.json) and [optimization record](docs/human-optimization-512-20261018.json). All outputs at all nine sizes match the complete independently validated C++ row sets.
