"""Compare native GGCAT selection with an independent per-genome counter.

Run after `cargo build -p ggcat_cmdline`:
    python3 tests/all_genomes.py
Set GGCAT_BIN to exercise a different build.
"""

import collections
import compression.zstd as zstd
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("GGCAT_BIN", ROOT / "target/debug/ggcat")).resolve()


def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def kmers(records, k, forward=False):
    counts = collections.Counter()
    for seq in records:
        seq = seq.upper()
        for i in range(len(seq) - k + 1):
            word = seq[i:i + k]
            if set(word) <= set("ACGT"):
                counts[word if forward else min(word, rc(word))] += 1
    return counts


def fasta(path, records):
    path.write_text("".join(f">c{i}\n{s}\n" for i, s in enumerate(records)))


def read_fasta(path):
    records = []
    for line in zstd.decompress(path.read_bytes()).decode().splitlines():
        if line.startswith(">"):
            records.append("")
        elif line:
            records[-1] += line
    return records


class AllGenomePresenceTests(unittest.TestCase):
    def check_graph(self, genomes, k, *, forward=False, buckets=2, threads=2, grouped=False,
                    native=True, extra=()):
        with tempfile.TemporaryDirectory(prefix="ggcat-all-genomes-") as directory:
            directory = Path(directory)
            paths = []
            mapping = []
            for genome, records in enumerate(genomes):
                groups = [[s] for s in records] if grouped else [records]
                for contig, sequences in enumerate(groups):
                    path = directory / f"g{genome}-{contig}.fasta"
                    fasta(path, sequences)
                    paths.append(path)
                    mapping.append(f"g{genome}\t{path}\n")
            output = directory / "unitigs.fasta"
            cmd = [str(BIN), "build", "-k", str(k), "--minimizer-length", str(min(3, k - 2)),
                   "-b", str(buckets), "-j", str(threads), "-m", "0.1", "-t", str(directory / "tmp"),
                   "-o", str(output)]
            cmd += ["--present-in-all-colors"] if native else ["-s", "2"]
            if forward:
                cmd.append("-f")
            if grouped:
                mapping_path = directory / "mapping.tsv"
                mapping_path.write_text("".join(mapping))
                cmd += ["-d", str(mapping_path)]
            else:
                cmd += list(map(str, paths))
            cmd += list(extra)
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout[-4000:] + result.stderr[-4000:])
            counts = [kmers(records, k, forward) for records in genomes]
            all_words = set().union(*counts)
            expected = {w for w in all_words if all(c[w] >= 1 for c in counts)} if native else {
                w for w in all_words if sum(c[w] for c in counts) >= 2
            }
            observed = set(kmers(read_fasta(output), k, forward))
            self.assertEqual(observed - expected, set(), f"unexpected k-mers, k={k}")
            self.assertEqual(expected - observed, set(), f"missing k-mers, k={k}")

    def test_repeats_missing_colors_and_contig_boundaries(self):
        self.check_graph([
            ["ACGTTGCATGTCGCATGATGCATGAGAGCT", "ACGTTGC", "AAAAAAA", "AAAAAAA"],
            ["ACGTTGCATGTCGCATGATGCATGAGAGCT", "GGGGGGG"],
            ["ACGTTGCATGTCGCATGATGCATGAGAGCT", "AAAAAAA"],
        ], 7)

    def test_grouped_files_are_one_genome(self):
        self.check_graph([
            ["ACGTTGCATGTCGCATGATGCATGAGAGCT", "ACGTTGC"],
            ["ACGTTGCATGTCGCATGATGCATGAGAGCT", "CCCCCCC"],
        ], 7, grouped=True)

    def test_shared_kmers_remain_eligible_at_any_copy_count(self):
        sequence = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        self.check_graph([[sequence, sequence], [sequence] * 3, [rc(sequence)]], 7)

    def test_reverse_complements_and_ambiguity(self):
        sequence = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
        for forward in (False, True):
            with self.subTest(forward=forward):
                self.check_graph([[sequence + "NNN" + sequence[:8].lower()], [rc(sequence)]],
                                 7, forward=forward)

    def test_ascii_reader_normalizes_all_ambiguity_codes(self):
        sequence = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
        separated = ''.join(sequence + base for base in "RYSWKMBDHVNXryswkmbdhvnx") + sequence.lower()
        self.check_graph([[separated], [rc(separated)]], 7)

    def test_empty_genome(self):
        self.check_graph([["ACGTTGCATGTCGCATGATGCATGAGAGCT"], [""]], 7)

    def test_one_genome(self):
        self.check_graph([["ACGTTGCATGTCGCATGATGCATGAGAGCT", "ACGTTGC"]], 7)

    def test_randomized_against_bruteforce(self):
        rng = random.Random(163)
        reference = "".join(rng.choices("ACGT", k=500))
        genomes = [
            [reference, reference[45:83], "A" * 35],
            [reference[:210] + "NNN" + reference[210:], reference[160:190]],
            [rc(reference), "".join(rng.choices("ACGT", k=75))],
        ]
        for k in (3, 4, 5, 6, 7, 8, 15, 31, 33, 63, 64, 67):
            with self.subTest(k=k):
                self.check_graph(genomes, k, buckets=3, threads=4)

    def test_short_input_parallel_shutdown(self):
        # Workers can start after the sole shared address has already closed.
        # They must still participate in final bucketing barriers.
        rng = random.Random(163)
        sequence = "".join(rng.choices("ACGT", k=500))
        for attempt in range(20):
            with self.subTest(attempt=attempt):
                self.check_graph([[sequence], [sequence], [rc(sequence)]],
                                 5, buckets=3, threads=4)

    def test_existing_abundance_mode(self):
        self.check_graph([["ACGTTGCATGTCGCATGATGCATGAGAGCT"], ["ACGTTGC"]], 7, native=False)


if __name__ == "__main__":
    unittest.main()
