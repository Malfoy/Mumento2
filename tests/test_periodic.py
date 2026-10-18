"""Independent full-pipeline oracles for the compact periodic-seed path."""
import compression.zstd as zstd
import os
import random
from pathlib import Path
import subprocess
import tempfile
import unittest

from oracle import expected_matches, rc, verify_output

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get('MUMENTO2_BIN', ROOT / 'target/debug/mumento2')).resolve()
GGCAT = Path(os.environ.get('GGCAT_BIN', ROOT / 'ggcat/target/debug/ggcat')).resolve()


class PeriodicTests(unittest.TestCase):
    def check_genomes(self, genomes, k):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            inputs = []
            for g, contigs in enumerate(genomes):
                path = directory / f'g{g}.fa'
                path.write_text(''.join(f'>c{i}\n{s}\n' for i, s in enumerate(contigs)))
                inputs.append(path)
            output = directory / 'out.zst'
            process = subprocess.run([str(BIN), '--ggcat', str(GGCAT), '-k', str(k),
                '-j', '3', '-m', '0.1', '-o', str(output), *map(str, inputs)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90)
            self.assertEqual(process.returncode, 0, process.stderr.decode())
            verify_output(self, zstd.decompress(output.read_bytes()).decode(), genomes,
                          k, None, False, expected_matches(genomes, k))

    def test_repeated_proto_fragments_keep_inner_mums_at_occurrence_boundaries(self):
        rng = random.Random(20261018512)
        for cut in (7, 31, 61, 63, 64, 65, 127, 129):
            core = ''.join(rng.choices('ACGT', k=290))
            reference = ['AA' + core + 'CC', 'TG' + core[cut:-7] + 'GA',
                         rc('CT' + core[7:-cut] + 'AG')]
            genomes = [reference, ['GT' + core + 'AC'], [rc('CA' + core + 'TG')]]
            with self.subTest(cut=cut): self.check_genomes(genomes, 7)

    def test_repeated_nonperiodic_contexts_and_physical_palindromes(self):
        rng = random.Random(20261018)
        for case in range(48):
            core = ''.join(rng.choices('ACGT', k=15 + case % 13))
            k = [3,5,7,9][case % 4]
            genome = ['A' + core + 'C', 'G' + core + 'T', 'A' + core + 'G']
            if case % 3 == 0:
                genomes = [genome, [rc(genome[0]), genome[1]], [genome[0]]]
            elif case % 3 == 1:
                genomes = [genome, ['T' + core + 'C', genome[1]], [core]]
            else:
                palindrome = core + rc(core)
                genomes = [['C' + palindrome + 'G'], ['C' + palindrome], [palindrome]]
            with self.subTest(case=case): self.check_genomes(genomes, k)

    def check_period(self, motif):
        p = len(motif)
        k = 2 * p + 1
        core = motif * 5 + motif[:p // 2]
        variants = [
            [[core], [core], [rc(core)]],
            [[core], [core + motif], [core]],
            [['G' + core + 'T'], ['C' + core + 'A'], [rc('A' + core + 'C')]],
            [[core, core[:-1]], [rc(core)], [core]],
            [[core], [core, core], [rc(core)]],
            [['AC' + core + 'TG'], [rc('AC' + core + 'TG')], ['AC' + core + 'TG']],
            [[core + 'N' + motif * 2], [core], [core.lower()]],
            [[core], [core[p:] + 'N' + core[:p + k]], [rc(core)]],
            [['A' + core + 'CTG', 'G' + core + 'TCA'], ['A' + core + 'CTG'], [rc('A' + core + 'CTG')]],
        ]
        for ordinal, genomes in enumerate(variants):
            with self.subTest(period=p, variant=ordinal), tempfile.TemporaryDirectory() as directory:
                directory = Path(directory)
                inputs = []
                for g, contigs in enumerate(genomes):
                    path = directory / f'g{g}.fa'
                    path.write_text(''.join(f'>c{i}\n{s}\n' for i, s in enumerate(contigs)))
                    inputs.append(path)
                output = directory / 'out.zst'
                process = subprocess.run([str(BIN), '--ggcat', str(GGCAT), '-k', str(k),
                    '-j', '3', '-m', '0.1', '-o', str(output), *map(str, inputs)],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90)
                self.assertEqual(process.returncode, 0, process.stderr.decode())
                expected = expected_matches(genomes, k)
                verify_output(self, zstd.decompress(output.read_bytes()).decode(), genomes,
                              k, None, False, expected)

    def test_period_1_homopolymers(self): self.check_period('A')
    def test_period_2_reverse_complement_palindromes(self): self.check_period('AT')
    def test_period_3_phase_and_strand(self): self.check_period('ACG')
    def test_period_4_phase_and_strand(self): self.check_period('AACT')
    def test_period_5_phase_and_strand(self): self.check_period('ATTCC')
    def test_period_6_phase_and_strand(self): self.check_period('AACCGT')
    def test_period_7_phase_and_strand(self): self.check_period('AACTGGC')
    def test_period_8_phase_and_strand(self): self.check_period('AAAACGTT')


if __name__ == '__main__': unittest.main()
