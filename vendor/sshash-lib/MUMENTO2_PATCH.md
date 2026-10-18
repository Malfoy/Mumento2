Based on COMBINE-lab `sshash-lib` 0.7.1 (BSD-3-Clause), with the local
`../REINDEER/vendor/sshash-rs` acceleration fork imported. The exact imported
source hashes are recorded in `benchmarks/human/optimization-sshash-20261018/reindeer-source.json`.

The REINDEER fork supplies PTRHash for minimizer and skew MPHFs and matching
SIMD minimizers for construction and queries. Mumento2 adds an entry point for
precomputed minimizers so the existing streaming verifier can preserve its
exact DNA checks and alignment skips. SIMD batches cover at most 128 KiB of
windows per worker, avoiding a full-document minimizer array. Packed document views can now supply both SIMD minimizers and streaming
windows directly in the ACTG two-bit alphabet. Initial windows are read in
word-sized pieces; adjacent windows retain the rolling k-mer state. Scalar
ASCII minimizer checks, unaligned packed-view tests and debug point-query
equality checks cover this adapter. The locator now uses the same re-exported packed-seq and seq-hash interfaces
as SSHash, so compact queries borrow the existing packed genome buffer.
No packed-buffer copy or ASCII reconstruction is needed at that interface.
Legacy serialized packed-genome bytes remain compatible and are tested. Dictionary file
formats follow the fork; older standalone SSHash dictionaries must be rebuilt.

The needletail dependency uses `default-features = false`. Mumento2 supplies
already parsed proto-MUM sequences and uses its existing helicase input layer;
this avoids bringing a second, incompatible liblzma native library into the
dependency graph. Build phase wall times are printed for profiling. See
`UPSTREAM.md` for the imported fork's pinned upstream revision and changes.

Packed jumps preserve immutable dictionary bucket evidence while discarding
position-dependent extension and singleton shortcuts. Mumento2 uses the
minimizer-absence flag to omit only windows whose minimizer is provably absent;
this does not change the dictionary or its exact DNA verification.
