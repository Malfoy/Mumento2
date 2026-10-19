# Correctness test design

The suite contains **237 library tests**, **8 driver tests**, and **205 end-to-end tests** (193 discovery methods plus 12 pipeline methods). There are **183 new, individually named scenarios in each layer**: 127 discovery cases, 28 incremental workflows, and 28 rejection/atomicity cases. The earlier 21 library and 6 pipeline tests are retained. Three additional pipeline tests exercise seed-subset coverage, threshold matches absent from genome zero, internal fragments, and rejected unsafe subsets. Driver tests exhaustively check the subset support bound for collections of up to ten genomes and reject unsafe settings. Two additional library tests compare every emitted fragment with seed skipping enabled and disabled across k-mer sizes, strands, chunks, repeats, ambiguity breaks, and index additions. Additional tests cover batched context encoding, truncated scratch records, vectorizable DNA validation across all byte values, file-list path handling, 100,000-entry lists, and an end-to-end compressed-input list with spaces and duplicate rejection. Five additional library tests check bitmap boundaries, inline-cache equivalence, batch scheduling, pruned seed islands, and restoration of raw candidates when extending a union index. All-genome discovery corpus cases also run with pruning disabled and forced batch sizes one and two against the independent oracle. A scenario is a separate Rust `#[test]` and Python `unittest` method, not a repetition counted inside a loop.

An exhaustive physical-interval oracle checks the bounded three-occurrence reducer across duplicate mappings, nested intervals, opposite strands and multiple shards. Separate wide-coordinate checks guard mapping arithmetic near u64 limits.

Additional compact-index tests compare SIMD batches with exact point queries across the 128-KiB batch boundary, mutations, ambiguity breaks, reverse strands and alignment jumps. Packed streaming tests cover k=31/33/61/63. The vendored SSHash suite adds scalar ASCII comparisons for packed minimizer values over odd/even k, m=1–31, unaligned starts and short/long views; it has 131 passing tests and two pre-existing ignored tests. Bitmap span tests compare both directions against a per-position boundary oracle. End-to-end periodic and overlapping-repeat cases check that occurrence-boundary skipping preserves smaller valid MUMs.

All fixtures are small synthetic sequences. The runner does not use bacterial benchmark data, collect performance comparisons, or run benchmark repetitions.

## Run

Requires the project build dependencies and Python **3.14+** (`compression.zstd`) on Linux (advisory-lock tests use `fcntl`); no Python packages need installing.

```sh
bash tests/run.sh
```

Run one named test when investigating a failure:

```sh
cargo test --locked --manifest-path k2rmini/Cargo.toml --lib update_new_only_region_second_batch
PYTHONPATH=tests python3 -m unittest test_corpus.CorpusTests.test_update_new_only_region_second_batch -v
```

The end-to-end tests invoke the actual Mumento2 executable and modified GGCAT, then decompress and inspect their files. `MUMENTO2_BIN` and `GGCAT_BIN` can select other builds; debug builds are used by default. No graph or locator is mocked.

## Expected results and assertions

[generate_cases.py](generate_cases.py) defines the reviewed scenarios, materializes their complete genomes and expected matches in [corpus/cases.json](corpus/cases.json), and generates the Rust test registrations. `--check` verifies that committed artifacts match the definitions. Regenerate after deliberately changing a scenario:

```sh
python3 tests/generate_cases.py
python3 tests/generate_cases.py --check
```

The Python oracle enumerates all unambiguous substrings and builds physical occurrence tables, without using graph paths, seeds, hashes, or alignment ranges. It checks support, single-copy status in every genome containing the match, and common left/right extension across the full support set. Palindromes have one physical occurrence and both possible extension orientations. A second, independent Rust oracle searches each candidate in the original sequences and checks the committed expected results during library tests. For the >128-KiB reader fixtures, both layers verify that genomes have identical contig multisets up to strand, then count exact occurrences of complete contigs: every proper unique substring extends inside its corresponding contig, so only complete unique contigs can be maximal. This gives a size-linear oracle without enumerating billions of substrings.

Every successful case checks complete match-set equality, not just the number of rows. Assertions also cover lengths, contig-relative coordinates, DNA reconstruction on both strands, genome order, missing-value consistency, owner orientation, duplicate suppression, deterministic coordinate order, and contig names/lengths. CLI discovery tests inspect the proto FASTA and phase-timing schema without asserting performance. All-genome library cases additionally compare the generalized API with the original streaming API.

Every update step is checked against the oracle. Old source files are deleted before adding genomes. After the final step, a separate full rebuild must produce the same matches. Tests check source-read counts, inherited parameters, genome counts, threshold changes, source preservation when branching, and seed-index growth when nothing novel is added.

Rejected operations must fail for the expected reason, preserve the existing output, and leave the active checkpoint pointer unchanged. Checkpoints are created by a successful real run before being damaged or modified for a failure test.

## Coverage

| Area | Adversarial cases |
| --- | --- |
| Threshold support | N=1, exact N, more/fewer than N, all genomes, missing initial/final genomes, nonadjacent support, disjoint and overlapping support sets, nested support changing at either or both ends |
| Uniqueness | Whole-string repeats, repeats in an extra genome, prefix/interior/suffix decoys, overlapping/tandem repeats, duplicated contigs, reverse-complement copies, palindromes, no unique k-mer anywhere in a valid full MUM |
| Exact maximality | Boundary/interior SNPs, multiple SNPs, insertion, deletion, inversion, rearranged contigs, truncated/overlapping proto ranges, valid internal fragments, ambiguous DNA boundaries |
| Coordinates and strands | Reversed first/later genomes, alternating strands, forward-only searches, both orientations at one physical position, leading ambiguity, contig offsets, empty records, repeated names |
| Index boundaries | k=3/4/15/29/31/32/33/63/64/65/67/129, match lengths k−1/k/k+1, minimizer length 1 and k, chunks 1/6/7/8/29/64/65, one worker and excess workers |
| FASTA/FASTQ I/O | Plain/gzip/xz/zstd and mixed codecs, concatenated compression members/frames with a short first read, wrapping at 1/7/63/64/65, CRLF and wrapped CRLF, missing final newline, Unicode/descriptive headers, quality strings containing `@`/`+`, CR/LF crossing SIMD blocks, headers and sequences exceeding the reader buffer, long quality-skip paths |
| Empty/short inputs | First/last/all documents empty, ambiguity-only documents, many short contigs, no contig reaching k, empty records preserving contig IDs |
| Incremental discovery | Preserve/shorten/split/remove old matches, later repeats, old below-threshold regions becoming eligible, matches present only in additions, discoveries delayed to a second batch, batched/sequential additions, changed N, inherited options, branching, no novel k-mers |
| Rejection and atomicity | Invalid N/minimizers/workers/chunks, absent/repeated inputs, incompatible update parameters, active locks, existing destinations, invalid/path-traversing CURRENT, absent or damaged manifests/seeds/packed genomes/context records, unsupported versions, absent output parent |

These checks establish correctness on the enumerated synthetic cases; they do not establish large-dataset memory use or speed. Low-N scaling and production-dataset performance remain separate benchmark work.

The first-level alignment cache is checked against interval boundaries, both strands, deliberate slot collisions and wide coordinates. Packed jump tests revisit identical minimizers to catch unsafe singleton shortcuts. SIMD-batch tests compare every window omitted by absent-minimizer skipping with independent exact dictionary point queries.

Cross-proto alignment skips are compared with disabled skips for complete context-set equality. Automatic graph seed-prefix coverage is checked against a full-input graph, including a seventeenth document that splits matches. An exhaustive substring oracle guards MUMs crossing nested repeat ends. Monotonic compact-scan alignment caches are compared with both the full-history cache and basewise extension on both strands.

A packed-DNA compatibility test restores bytes produced by the previous packed-seq implementation and checks every base plus the serialized encoding at empty, partial-byte, word-boundary and long lengths.

Directional repeat-stop tests check unchanged two-copy mappings, canonical strand changes, changed partners and relative strands, three copies and overlapping proto groups. Forward and reverse spans are checked against separate per-position bit oracles, and complete context-set equality includes queries beginning inside repeat runs.

Compact repeat-context leaves preserve full 64-bit locations and absent-location markers across wide-coordinate and depth boundaries.

Indexed-seed coordinate tests cover both canonical orientations and query directions, multiple and empty contigs, mixed case and ambiguity boundaries. Full run-boundary checks remain in refinement and final uniqueness validation.

An independent complete affine occurrence-set oracle checks every directional repeat-stop bit across copy counts, nested truncations, mutations, reverse copies and contig/ambiguity boundaries. Conflicting physical groups remain conservative. Separate multi-copy scans compare every emitted context with disabled skipping, including unresolved repeated queries and differing repeat tails.

Alignment-witness tests compare shortcut results directly with the original context tree over every represented occurrence, both directions, truncated queries, ambiguous runs, nested repeats and coordinates above 32 bits. Full context sets are also compared with an exhaustive scan with witness hints disabled.

Leaf occurrence directory tests cover exact record IDs, absent keys, word boundaries, random sparse keys, conflicting records and wide-coordinate fallback. Existing alignment-witness oracles exercise the directory in both orientations.

Context singleton skip tests exhaustively check the 64-base proof reserve across word boundaries and both scan directions, and reject conflicting full mappings, over-depth leaves, two-copy leaves, missing leaves and wrong original keys. Independent all-seed oracles disable both original and new shortcuts while checking diverging repeat copies, ambiguity, both strands, k7,31,33,63 and alignments near context boundaries.
