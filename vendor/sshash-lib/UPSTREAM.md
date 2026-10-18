SSHash-rs Cargo dependency derived from COMBINE-lab/sshash-rs revision
2480b1cbabce907d3dcaa336714a75d11b0d0ecb (0.7.1).

Local changes replace the internal PHast MPHFs with PTRHash, use SIMD
minimizers consistently in the builder and lookup paths, accept REINDEER
k-mer sizes including even values and 64, and expose SPSS construction
and batched querying through the original Dictionary. This is the SSHash
library dependency, not a second REINDEER application source.

README.md is the upstream reference documentation (with its directory tree
marked as text for rustdoc). Its PHast dependency description refers to the
upstream version; this vendored Cargo manifest resolves PTRHash instead.
