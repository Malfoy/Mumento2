# Mumento2

**Work in progress.** A Rust implementation inspired by [Mumemto](https://github.com/vikshiv/mumemto), using modified [GGCAT](https://github.com/algbio/ggcat) and [K2Rmini](https://github.com/Malfoy/K2Rmini) components to discover and validate maximal unique matches across genome collections.

Mumento2 searches both strands and respects contig boundaries. It reports maximal exact matches that occur once in every participating genome. By default a match must occur in all input genomes; a configurable minimum genome count supports matches shared by a subset, including matches absent from the first genome. Independently valid internal fragments survive when a longer proto-MUM fails validation.

Reusable state preserves union candidates, packed genomes and exact position-validity evidence. Adding genomes revalidates those candidates, allowing matches to shorten or split and newly qualifying matches to appear. State generations are immutable and committed atomically. Collinear filtering and a full-text query index are not currently provided.

The streaming all-genome path uses parallel exact seed scanning and extension, batched scratch I/O, inline alignment caches and progressive pruning of individual proto seed positions. A seed position is removed only when exact evidence proves it absent from an already processed genome. Reusable state retains raw candidates needed for later updates.

Native CPU optimization, optimization level 3, fat LTO and a single code-generation unit are enabled for optimized builds. Profile-guided compilation support and a validated LLVM 22 training profile are included. Profiles are tied to source/compiler versions and workload; improvements on other datasets are not established.

Inputs may be plain or compressed FASTA/FASTQ, with each file representing one genome and each record a separate contig. Large collections can use a file list. Output is a compressed TSV containing match lengths, zero-based forward-strand starts, strands and contig IDs, plus contig metadata. Phase timings expose graph construction, indexing, scanning/extension, reduction and output costs.

The test suite contains **213 library tests, 7 driver tests and 193 end-to-end tests**, including independent small-input oracles, internal-fragment recovery, reverse strands, repeats, thresholds and incremental updates. See [test coverage](tests/README.md).

In single measured runs on a 32-logical-CPU Intel i9-13950HX, the retained algorithm completed 18,800 bacterial assemblies in **26.66 seconds**, versus **34.44 seconds** before progressive pruning. On a fixed 10,000-assembly collection, the profile-guided executable took **16.08 seconds**, versus **16.96 seconds** for the ordinary optimized executable. Both produced the same fully validated results. These are workload-specific measurements, not statistical performance estimates; threshold/update mode and human datasets require separate measurements.

Human measurements use real HPRC Release 2 **chromosome 19 haplotypes**, one per distinct individual, following the chromosome scaling experiment in the [Mumemto paper](https://doi.org/10.1186/s13059-025-03644-0). These are chromosome inputs, approximately 61 Mb each; uniqueness is checked within those inputs. Nested prefixes were tested once each, with both strands, strict presence in every input, Rust k = 63 and minimum MUM length = 63 bp.

| Haplotypes | C++ seconds | C++ peak RSS (GiB) | Validated C++ MUMs | Rust |
|---:|---:|---:|---:|---|
| 2 | 103.82 | 2.58 | 57,456 | 600 s timeout |
| 4 | 125.58 | 8.32 | 91,348 | 600 s timeout |
| 8 | 147.97 | 19.34 | 120,060 | 600 s timeout |
| 16 | 221.77 | 41.48 | 140,229 | 600 s timeout |
| 32 | Stopped at memory cap | 48.06 at abort | — | Not attempted |

The CPU budget was 32. C++ Mumemto 1.4.2 has a single-threaded core; its partitioned workflow used 1, 3, 7, 15 and 31 workers respectively, with merging included in timing. The 32-input partitioned attempt exceeded the 48 GiB process-tree RSS budget, ending the series. Rust requested 32 threads but remained in serial reference verification. Downloads, staging and validation were excluded from timing. This is a newer panel and tool version than the paper's exact experiment.

All completed C++ outputs passed sequence, coordinate, contig-boundary, joint maximality and independent physical uniqueness checks. Rust produced no completed human output, so cross-tool equality is unestablished. Reference self-alignment over repetitive proto-MUM occurrences is the main unresolved bottleneck: at four inputs this phase alone consumed at least 586 seconds. The initial two-input k = 31 trial also timed out in Rust while C++ completed in 104.10 seconds. These observations do not establish whole-human performance. Input provenance, binary hashes, resource limits and separate experimental configurations are recorded in the [human benchmark results](docs/human-benchmark-20261018.json).
