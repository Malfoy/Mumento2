# Mumento2

**Work in progress.** A Rust implementation inspired by [Mumemto](https://github.com/vikshiv/mumemto), using modified [GGCAT](https://github.com/algbio/ggcat) and [K2Rmini](https://github.com/Malfoy/K2Rmini) components to discover and validate maximal unique matches across genome collections.

Mumento2 searches both strands and respects contig boundaries. It reports maximal exact matches that occur once in every participating genome. By default a match must occur in all input genomes; a configurable minimum genome count supports matches shared by a subset, including matches absent from the first genome. Independently valid internal fragments survive when a longer proto-MUM fails validation.

Reusable state preserves union candidates, packed genomes and exact position-validity evidence. Adding genomes revalidates those candidates, allowing matches to shorten or split and newly qualifying matches to appear. State generations are immutable and committed atomically. Collinear filtering and a full-text query index are not currently provided.

The streaming all-genome path uses parallel exact seed scanning and extension, batched scratch I/O, inline alignment caches and progressive pruning of individual proto seed positions. A seed position is removed only when exact evidence proves it absent from an already processed genome. Reusable state retains raw candidates needed for later updates.

Native CPU optimization, optimization level 3, fat LTO and a single code-generation unit are enabled for optimized builds. Profile-guided compilation support and a validated LLVM 22 training profile are included. Profiles are tied to source/compiler versions and workload; improvements on other datasets are not established.

Inputs may be plain or compressed FASTA/FASTQ, with each file representing one genome and each record a separate contig. Large collections can use a file list. Output is a compressed TSV containing match lengths, zero-based forward-strand starts, strands and contig IDs, plus contig metadata. Phase timings expose graph construction, indexing, scanning/extension, reduction and output costs.

The test suite contains **230 library tests, 8 driver tests and 205 end-to-end tests**, including independent small-input oracles, internal-fragment recovery, reverse strands, repeats, thresholds and incremental updates. See [test coverage](tests/README.md).

In single measured runs on a 32-logical-CPU Intel i9-13950HX, the retained algorithm completed 18,800 bacterial assemblies in **26.66 seconds**, versus **34.44 seconds** before progressive pruning. On a fixed 10,000-assembly collection, the profile-guided executable took **16.08 seconds**, versus **16.96 seconds** for the ordinary optimized executable. Both produced the same fully validated results. These are workload-specific measurements, not statistical performance estimates; threshold/update mode and human datasets require separate measurements.

Human measurements use real HPRC **chromosome 19 haplotypes**, following the chromosome scaling experiment in the [Mumemto paper](https://doi.org/10.1186/s13059-025-03644-0). Inputs average approximately 61 Mb; uniqueness is checked within chromosome 19. The 512-input collection represents 276 individuals, preserving the Release 2 prefix and adding public HPRC R3 assemblies. Whole-human-genome performance has not been measured.

The best measured Rust build is **v24**, on all **512 haplotypes**, with both strands, strict presence in every input, k = 63 and minimum MUM length = 63 bp. Each distinct Rust candidate is measured once; the saved C++ measurement is reused after checking input manifests and executable hashes.

| Tool | Wall seconds | CPU seconds | Peak process-tree RSS (GiB) |
|---|---:|---:|---:|
| C++ Mumemto 1.4.2 | 6,931.30 | 106,287.83 | 46.23 |
| Rust Mumento2 v24 | 105.12 | 2,845.34 | 7.37 |
| C++ ÷ Rust | **65.94× faster** | **37.36× less CPU time** | **6.27× less RAM** |

All **254,892 MUMs** match C++ exactly, including lengths, coordinates and strands. The comparison inherits sequence, contig-boundary, joint maximality and physical uniqueness audits only through exact full row-set equality with independently validated C++ output.

For stateless all-genome searches, graph construction uses up to 16 seed genomes by default. Every qualifying MUM occurs in that prefix; validation still scans every input genome. Threshold searches and reusable states retain full candidate coverage. The native Rust SSHash proto-MUM index incorporates PTRHash/SIMD acceleration from REINDEER. Exact aligned intervals can skip across proto boundaries only until the next repeat entry or change in its physical occurrence mappings. Separate forward and reverse stop maps coalesce unchanged complete reference occurrence sets. Unresolved contexts compare DNA only to the next occurrence-set boundary. Other retained optimizations bound repeat evidence to three distinct oriented mappings per genome, share packed DNA buffers between the scanner and SSHash, use 32-byte repeat-context nodes and exact alignment caches, preserve dictionary-only bucket information at query jumps and skip absent-minimizer windows.

Rust uses 32 native threads. Saved C++ runs use 16 concurrent partition workers, including merging, constrained by the 48 GiB memory cap. These measurements used a 32-logical-CPU Intel i9-13950HX. CPU is summed user plus system time, including children. Timing includes graph construction, indexing, scanning, reduction and output; downloads, builds, staging and correctness audits are excluded. Single observations are not statistical performance estimates.

See the [latest 512-input comparison](docs/human-global-512-latest-20261018.json) and [optimization record](docs/human-optimization-512-20261018.json). The [archived scaling comparison](docs/human-global-all-sizes-latest-20261018.json) records the same earlier frozen v13 build at 2–512 inputs; those smaller sizes have not been remeasured with v24.
