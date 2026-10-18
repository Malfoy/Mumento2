# Mumento2

**Work in progress.** A Rust implementation inspired by [Mumemto](https://github.com/vikshiv/mumemto), using modified [GGCAT](https://github.com/algbio/ggcat) and [K2Rmini](https://github.com/Malfoy/K2Rmini) components to discover and validate maximal unique matches across genome collections.

Mumento2 searches both strands and respects contig boundaries. It reports maximal exact matches that occur once in every participating genome. By default a match must occur in all input genomes; a configurable minimum genome count supports matches shared by a subset, including matches absent from the first genome. Independently valid internal fragments survive when a longer proto-MUM fails validation.

Reusable state preserves union candidates, packed genomes and exact position-validity evidence. Adding genomes revalidates those candidates, allowing matches to shorten or split and newly qualifying matches to appear. State generations are immutable and committed atomically. Collinear filtering and a full-text query index are not currently provided.

The streaming all-genome path uses parallel exact seed scanning and extension, batched scratch I/O, inline alignment caches and progressive pruning of individual proto seed positions. A seed position is removed only when exact evidence proves it absent from an already processed genome. Reusable state retains raw candidates needed for later updates.

Native CPU optimization, optimization level 3, fat LTO and a single code-generation unit are enabled for optimized builds. Profile-guided compilation support and a validated LLVM 22 training profile are included. Profiles are tied to source/compiler versions and workload; improvements on other datasets are not established.

Inputs may be plain or compressed FASTA/FASTQ, with each file representing one genome and each record a separate contig. Large collections can use a file list. Output is a compressed TSV containing match lengths, zero-based forward-strand starts, strands and contig IDs, plus contig metadata. Phase timings expose graph construction, indexing, scanning/extension, reduction and output costs.

The test suite contains **230 library tests, 8 driver tests and 205 end-to-end tests**, including independent small-input oracles, internal-fragment recovery, reverse strands, repeats, thresholds and incremental updates. See [test coverage](tests/README.md).

Human measurements use real HPRC **chromosome 19 haplotypes**, following the chromosome scaling experiment in the [Mumemto paper](https://doi.org/10.1186/s13059-025-03644-0). Inputs average approximately 61 Mb; uniqueness is checked within chromosome 19. The 512-input collection represents 276 individuals, preserving the Release 2 prefix and adding public HPRC R3 assemblies. Whole-human-genome performance has not been measured. See the [input provenance and selection](docs/human-chr19-dataset-20261018.md).

The published **v24** build was measured on every collection size below, with both strands, strict presence in every input, k = 63 and minimum MUM length = 63 bp. The same frozen executable is used at all nine sizes. Each size has one measured Rust run; the completed v24 512-input run and all saved C++ measurements are reused.

| Haplotypes | C++ wall (s) | Rust wall (s) | Wall gain | C++ CPU (s) | Rust CPU (s) | CPU gain | C++ RAM (GiB) | Rust RAM (GiB) | RAM gain |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 103.82 | 11.05 | 9.39× | 103.68 | 55.32 | 1.87× | 2.58 | 1.68 | 1.54× |
| 4 | 125.58 | 11.41 | 11.01× | 372.91 | 63.50 | 5.87× | 8.32 | 1.91 | 4.37× |
| 8 | 147.97 | 12.11 | 12.22× | 1,018.21 | 84.90 | 11.99× | 19.34 | 2.53 | 7.65× |
| 16 | 221.77 | 14.13 | 15.70× | 3,022.01 | 138.57 | 21.81× | 41.48 | 3.82 | 10.87× |
| 32 | 454.91 | 15.39 | 29.57× | 6,393.35 | 227.42 | 28.11× | 44.24 | 6.51 | 6.80× |
| 64 | 868.70 | 21.43 | 40.54× | 13,076.85 | 392.97 | 33.28× | 44.63 | 7.02 | 6.36× |
| 128 | 1,737.04 | 32.72 | 53.09× | 26,275.72 | 722.16 | 36.38× | 45.84 | 7.15 | 6.41× |
| 256 | 3,464.91 | 55.86 | 62.03× | 52,896.79 | 1,380.54 | 38.32× | 46.45 | 7.35 | 6.32× |
| 512 | 6,931.30 | 105.12 | 65.94× | 106,287.83 | 2,845.34 | 37.36× | 46.23 | 7.37 | 6.27× |

Gains are C++ divided by Rust: values above one mean Rust uses less wall time, CPU time or RAM. These are single observations, not statistical performance estimates.

All **254,892 MUMs** match C++ exactly, including lengths, coordinates and strands. The comparison inherits sequence, contig-boundary, joint maximality and physical uniqueness audits only through exact full row-set equality with independently validated C++ output.

For stateless all-genome searches, graph construction uses up to 16 seed genomes by default. Every qualifying MUM occurs in that prefix; validation still scans every input genome. Threshold searches and reusable states retain full candidate coverage. The native Rust SSHash proto-MUM index incorporates PTRHash/SIMD acceleration from REINDEER. Exact aligned intervals can skip across proto boundaries only until the next repeat entry or change in its physical occurrence mappings. Separate forward and reverse stop maps coalesce unchanged complete reference occurrence sets. Unresolved contexts compare DNA only to the next occurrence-set boundary. Other retained optimizations bound repeat evidence to three distinct oriented mappings per genome, share packed DNA buffers between the scanner and SSHash, use 32-byte repeat-context nodes and exact alignment caches, preserve dictionary-only bucket information at query jumps and skip absent-minimizer windows.

Rust requests 32 native threads. Saved C++ runs used 1, 3, 7 and 15 partition workers at 2, 4, 8 and 16 inputs, respectively; runs from 32 inputs used 16 workers, including merging, constrained by the 48 GiB memory cap. These measurements used a 32-logical-CPU Intel i9-13950HX. CPU is summed user plus system time, including children. Timing includes graph construction, indexing, scanning, reduction and output; downloads, builds, staging and correctness audits are excluded. Single observations are not statistical performance estimates.

See the [complete comparison](docs/human-global-all-sizes-latest-20261018.json), [CSV results](docs/human-global-all-sizes-latest-20261018.csv), [512-input phase timings](docs/human-global-512-latest-20261018.json) and [optimization record](docs/human-optimization-512-20261018.json). All outputs at all nine sizes match the complete independently validated C++ row sets.
