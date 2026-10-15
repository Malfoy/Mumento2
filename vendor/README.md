# Shared dependencies

`helicase/` is Helicase 0.1.1 from
[`583e73568fc8c9cd17424170a875fef5afc9c4e1`](https://github.com/imartayan/helicase/tree/583e73568fc8c9cd17424170a875fef5afc9c4e1).
Both GGCAT and K2Rmini use this copy for FASTA/FASTQ parsing. Its gzip, xz,
and zstd decoders detect compression from input bytes. Synchronous reader
inputs no longer require `Send`, allowing borrowed archive entry readers.
The upstream MIT license is included.

`parallel-processor/` is version 0.2.7, from upstream revision
`b14611d64a34e57d52756276b5af69d551d1bd4a`. GGCAT patches crates.io to use
this copy. Bucket writers and readers use zstd, including formerly plain
lock-free buckets. Bucket headers use zstd skippable frames; checkpoint
indexes and payloads are compressed frames. The upstream MIT and Apache-2.0
licenses are included. These disk formats are incompatible with upstream
GGCAT intermediates.
