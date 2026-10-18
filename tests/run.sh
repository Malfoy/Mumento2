#!/usr/bin/env bash
# Small synthetic correctness tests only; no benchmark data or timed runs.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 tests/generate_cases.py --check
cargo test --locked --bin mumento2 -j "${MUMENTO2_TEST_BUILD_JOBS:-4}"
cargo test --locked --manifest-path k2rmini/Cargo.toml --lib -j "${MUMENTO2_TEST_BUILD_JOBS:-4}" -- --test-threads=4
cargo build --locked -j "${MUMENTO2_TEST_BUILD_JOBS:-4}"
cargo build --locked --manifest-path ggcat/Cargo.toml -p ggcat_cmdline -j "${MUMENTO2_TEST_BUILD_JOBS:-4}"
python3 -m unittest discover -s tests -p 'test_*.py' -v
python3 tests/pipeline.py -v
