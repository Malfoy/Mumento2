# Mumento2 handoff notes (2026-10)

Developer notes for the k≥20 / partial-MUM / updateable-merge / mumemto-I/O work.
Algorithm recap: [Mumento2 codebase explanation](486981e8-3202-4c9a-9d5b-b98f3e51df8c).

## Architecture (unchanged)

1. **GGCAT** builds colored de Bruijn proto-MUMs (shared k-mer unitigs).
2. **K2Rmini locate** remaps every proto k-mer into genomes, extends, and reduces to maximal unique matches.
3. **SSHash** is the compact proto dictionary on the all-genome graph path (`locate_graph_mums`).

Not provided: collinear blocks, full-text query index, classic mumemto merge interop.

## Changes in this handoff

### 0. Primary `-l`, optional advanced seed `-k`, and `.bumbl`

| Lever | Role |
|-------|------|
| **`-l` / `--min-match-len`** | **Primary knob.** Reported MUM floor; default discovery seed is `k = l`. |
| **`--seed-k` / `-k`** | Advanced override of seed/dictionary k. Default: `k = l`. Must satisfy **`l >= k`**. Legacy `-k` alone still means `k = l`. |

Completeness is set by the seed: only MUMs of length ≥ **k** can be discovered. Prefer `-l` alone (`k = l`). Setting `k < l` discovers denser seeds then drops short survivors — more RAM/CPU, never recovers MUMs shorter than k. Larger k (with `l = k`) → fewer shared k-mers → usually less RAM and faster locate.

**`-b` / `--binary`**: write mumemto-compatible `.bumbl` (LE: uint16 flags with length32 + optional partial, n_seqs, n_mums, uint32 lengths, int64 starts with -1 absent, MSB-first strands; no collinear blocks). Always global coords; **errors if combined with `--contig-coords`**. Still writes `{prefix}.lengths`. If `-o` ends in `.mums` / `.mums.zst` / `.tsv` / `.tsv.zst`, the path is rewritten to `.bumbl`.

Updates inherit both `k` and `min_length` from state (`STATE_VERSION = 3`).

### 1. SSHash for min length ~20 (merged with upstream `7670d40`)

Upstream commit **Keep small-k and incremental proto indexes on SSHash** lands the fuller fix:

- Minimizer `m = 31` at k=63; else clamp to `[1, min(k-1, 23)]` of floor(k/2), with **k≤31 → max(floor(2k/3), floor(k/2))** (k=21 → m=14). Override: `--sshash-m` / `-M`.
- New [`postings.rs`](../k2rmini/src/locate/stream/postings.rs): packed occurrence lists behind a distinct-seed SSHash SPSS for overlapping/union protos.
- Small-k, even-k, forward-only, threshold and incremental paths retain SSHash for **3 ≤ k ≤ 64**.

Our handoff keeps even-k on the **graph compact** path (`from_graph` no longer rejects even k) and layers mumemto I/O / relative `-N` / state docs on top. Prefer upstream’s SSHash architecture when in doubt; do not reintroduce dense hash-posting dictionary fallbacks.

### 2. Partial MUMs (`-N`)

Already implemented via `locate_shared_mums` / `reduce_shared`.

- Added **relative** `-N` like mumemto `-k`: `-N -1` ⇒ present in all but one genome.
- Zero remains invalid (`None` / omit flag = strict all-genome).
- Absent genomes use `.` in contig-coords mode; empty fields in mumemto global mode.
- Nested shorter matches with higher support remain maximal.

### 3. Updateable merge (mumento2 ↔ mumento2 only)

Use the existing state path; no classic mumemto `.thresh` / MUMs-of-MUMs merge.

```bash
# Initial updateable run
mumento2 -k 63 --state /path/state -o /path/out.mums.zst genomes_a/*.fa

# Add assemblies (union-set semantics; compared to fresh rebuild in tests)
mumento2 --update /path/state -o /path/out2.mums.zst genomes_b/*.fa

# Optional fork (preserve source state)
mumento2 --update /path/state --state /path/fork -o ... genomes_b/*.fa
```

`--mergeable` is an alias for `--state`.

#### On-disk schema (`STATE_VERSION = 3`)

| Path | Role |
|------|------|
| `CURRENT` | Points at active `generation-*` |
| `generation-*/state.json` | `k`, `min_length`, `m`, `forward_only`, `minimum_genomes` (`null` = all incl. future), `genomes`, `patterns` |
| `generation-*/protos.fa.zst` | Union candidate proto sequences |
| `generation-*/genomes/{id}.zst` | Packed DNA (`M2REF`), contig meta, ambiguity runs, validity fragments |
| `generation-*/contexts/{anchor}/{query}` | Exact position-validity evidence reused across updates |

TSV-only merge of two finished MUM files is **incorrect**: maximality and uniqueness need packed genomes + evidence. Always update via `--state`/`--update` (new FASTAs required for added assemblies).

No backwards compatibility with version-1/2 state dirs or classic mumemto merge metadata.

### 4. Mumemto-compatible I/O

| Mode | Flag | MUM rows | Sidecars |
|------|------|----------|----------|
| **Default** | (none) | 3-col global: `length\\tstarts\\tstrands` (contig-concat offsets), **no header**, zstd | `{prefix}.lengths` |
| Contig-relative | `--contig-coords` | 4-col: `length\\tstarts\\tstrands\\tcontigs`, **no header**, zstd | `{prefix}.lengths` + `{output}.contigs.tsv.zst` |
| Binary | `-b` / `--binary` | mumemto `.bumbl` (global starts; `-1` absent) | `{prefix}.lengths` |

`.lengths` format (plain text, mumemto-compatible):

```
path * total_bases
path contig_name length
...
```

Prefix is derived by stripping `.zst` / `.tsv` / `.mums` / `.bumbl` from `-o`.

## How to test

```bash
conda activate rust
export CARGO_TARGET_DIR=$PWD/target   # keep builds off tiny rootfs /tmp caches
cargo test --locked --bin mumento2 -j 4
cargo test --locked --manifest-path k2rmini/Cargo.toml --lib -j 4 -- --test-threads=4
cargo build --locked -p mumento2 -j 4
# GGCAT: prefer an existing optimized binary, e.g. ggcat/target/mumento/ggcat

export MUMENTO2_BIN=$PWD/target/debug/mumento2
export GGCAT_BIN=$PWD/ggcat/target/mumento/ggcat
/usr/bin/python3.14 -m unittest discover -s tests -p 'test_*.py' -v
/usr/bin/python3.14 tests/pipeline.py -v
```

Verified in this handoff: driver + K2Rmini lib tests (incl. `-l` / bumbl round-trip), pipeline tests for `-l` filtering and `-b` `.bumbl` (+ `.lengths`).

Chr19 k/l efficiency sweep (4 haplotypes, SLURM array, `/usr/bin/time -v`): see [`benchmarks/kl_sweep/README.md`](../benchmarks/kl_sweep/README.md).

Scratch for manual runs: [`artifacts/handoff-tests/`](../artifacts/handoff-tests/) (gitignored under `/artifacts/`).

## Known limits

- k ≤ 63; large collections at short k still grow the proto graph (more memory than k=63).
- Shared / update path is not SSHash-backed.
- Whole-genome (vs chr19) performance not re-benchmarked in this handoff.
- No collinear filtering / query index.

## Small-k performance campaign (`small-k-perf`)

Chr19 AGC 16-hap campaign notes live under [`benchmarks/small_k_perf/`](../benchmarks/small_k_perf/). Locked cheap defaults: adaptive SSHash `m` (`--sshash-m` / `-M` override), compact prune env-off, sparse `ReferenceSeeds` for wide coords, k-adaptive seed-genome cap. Heavier scan/build experiments were dropped.
