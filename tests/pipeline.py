"""End-to-end tests of the native GGCAT → K2Rmini Rust driver."""

import gzip
import lzma
import compression.zstd as zstd
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("MUMENTO2_BIN", ROOT / "target/debug/mumento2")).resolve()
GGCAT = Path(os.environ.get("GGCAT_BIN", ROOT / "ggcat/target/debug/ggcat")).resolve()


def rc(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class PipelineTests(unittest.TestCase):
    def test_repeated_prefix_skip_preserves_unique_mums_crossing_its_end(self):
        from oracle import expected_matches, verify_output
        rng = random.Random(170035)
        dna = lambda n: ''.join(rng.choices('ACGT', k=n))
        repeated = dna(120)
        lefts = [dna(30) for _ in range(3)]
        tails = [dna(40) for _ in range(4)]
        reference = [lefts[i] + repeated + tails[i] for i in range(3)]
        # This copy starts inside the repeated prefix and alone matches the
        # query beyond its end. A skip must leave a crossing seed to recover it.
        reference.append(repeated[40:] + tails[3])
        # Its right suffix also occurs in other documents/contigs: the first
        # crossing seed alone cannot establish this MUM's left uniqueness.
        reference.extend([repeated[45:] + tails[3], repeated[50:] + tails[3]])
        query = [lefts[0] + repeated + tails[3]]
        for reverse_reference in [False, True]:
            for reverse_query in [False, True]:
                genomes = [[rc(s) if reverse_reference else s for s in reference],
                           [rc(s) if reverse_query else s for s in query]]
                for k in [7, 31, 63]:
                    with self.subTest(k=k, rr=reverse_reference, rq=reverse_query), tempfile.TemporaryDirectory(prefix='mumento2-repeat-prefix-') as scratch:
                        scratch = Path(scratch)
                        inputs = self.write_genomes(scratch, genomes, 'g')
                        self.run_tool(scratch, 'out.zst', inputs, '-k', str(k), '--seed-genomes', '1')
                        text = zstd.decompress((scratch/'out.zst').read_bytes()).decode()
                        expected = expected_matches(genomes, k)
                        verify_output(self, text, genomes, k, None, False, expected)
                        self.assertTrue(any(len(m['sequence']) >= 120 and m['hits'][0][0] == 3 for m in expected))

    def test_automatic_bounded_seed_graph_validates_every_document(self):
        rng = random.Random(170034)
        core = ''.join(rng.choices('ACGT', k=1000))
        # The last document lies beyond the automatic seed prefix. Its contig
        # breaks must still restrict long candidates into their actual MUMs.
        genomes = [[core] for _ in range(16)] + [[core[:400], core[370:730], core[700:]]]
        with tempfile.TemporaryDirectory(prefix='mumento2-auto-seeds-') as scratch:
            scratch = Path(scratch)
            inputs = self.write_genomes(scratch, genomes, 'g')
            automatic, run = self.run_tool(scratch, 'auto.zst', inputs, '-k', '31')
            full, _ = self.run_tool(scratch, 'full.zst', inputs, '-k', '31', '--seed-genomes', '17')
            self.assertEqual(set(automatic), set(full))
            self.assertEqual({int(r.split('\t')[0]) for r in automatic}, {300, 360, 400})
            self.assertIn('Seed subset: 16/17', run.stderr)

    def test_input_list_runs_compressed_genomes_with_spaces_and_rejects_duplicates(self):
        rng = random.Random(180031)
        core = ''.join(rng.choices('ACGT', k=240))
        with tempfile.TemporaryDirectory(prefix='mumento2-list-') as directory:
            directory = Path(directory)
            inputs = self.write_genomes(directory, [[core], [rc(core)], [core]], 'genome with spaces')
            listing = directory / 'genomes.list'
            listing.write_bytes(('\r\n'.join(p.name for p in inputs) + '\r\n\r\n').encode())
            rows, _ = self.run_tool(directory, 'out.zst', [], '-k', '31', '--input-list', str(listing),
                                    '--keep-temp', '--temp-dir', str(directory), '--seed-genomes', '3')
            self.assertEqual(rows, ['240\t0,0,0\t+,-,+\t0,0,0'])
            graph_list = next(directory.glob('mumento2-*/graph-inputs.list'))
            self.assertEqual(graph_list.read_text().splitlines(), list(map(str, inputs)))
            listing.write_text(inputs[0].name + '\n' + inputs[0].name + '\n')
            run = subprocess.run([str(BIN), '-k', '31', '-o', str(directory/'rejected.zst'),
                                  '--input-list', str(listing)], capture_output=True, text=True, timeout=15)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn('more than once', run.stderr)
            self.assertFalse((directory/'rejected.zst').exists())

    def test_seed_subset_preserves_all_genome_matches_and_internal_fragments(self):
        rng = random.Random(170031)
        core = ''.join(rng.choices('ACGT', k=1000))
        genomes = [[core], [rc(core)], [core[:400], core[370:730], core[700:]]]
        with tempfile.TemporaryDirectory(prefix='mumento2-seed-subset-') as scratch:
            scratch = Path(scratch)
            inputs = self.write_genomes(scratch, genomes, 'g')
            full, _ = self.run_tool(scratch, 'full.zst', inputs, '-k', '31',
                                    '--seed-genomes', str(len(inputs)))
            self.assertEqual({int(r.split('\t')[0]) for r in full}, {300, 360, 400})
            automatic, _ = self.run_tool(scratch, 'automatic.zst', inputs, '-k', '31')
            self.assertEqual(set(automatic), set(full))
            for count in (1, 2):
                subset, _ = self.run_tool(scratch, f'subset-{count}.zst', inputs,
                                          '-k', '31', '--seed-genomes', str(count))
                self.assertEqual(set(subset), set(full))

    def test_seed_subset_preserves_threshold_matches_absent_from_first_genome(self):
        rng = random.Random(170032)
        a, b, c = [''.join(rng.choices('ACGT', k=n)) for n in (240, 210, 180)]
        genomes = [[a, c, c], [rc(a), b], [a, b], [rc(b)]]
        with tempfile.TemporaryDirectory(prefix='mumento2-seed-threshold-') as scratch:
            scratch = Path(scratch)
            inputs = self.write_genomes(scratch, genomes, 'g')
            for minimum, count in ((2, 3), (3, 2)):
                full, _ = self.run_tool(scratch, f'full-{minimum}.zst', inputs,
                                        '-k', '31', '-N', str(minimum))
                subset, _ = self.run_tool(scratch, f'subset-{minimum}.zst', inputs,
                                          '-k', '31', '-N', str(minimum), '--seed-genomes', str(count))
                self.assertEqual(set(subset), set(full))
                self.assertTrue(any(r.startswith('210\t.,') for r in subset))

    def test_seed_subset_rejects_incomplete_threshold_coverage_and_state(self):
        rng = random.Random(170033)
        core = ''.join(rng.choices('ACGT', k=100))
        with tempfile.TemporaryDirectory(prefix='mumento2-seed-reject-') as scratch:
            scratch = Path(scratch)
            inputs = self.write_genomes(scratch, [[core]]*4, 'g')
            for extra, message in [(['-N', '2'], 'too small'),
                                   (['--state', str(scratch/'state')], 'reusable')]:
                result = subprocess.run([str(BIN), '--ggcat', str(GGCAT), '-k', '31',
                                         '--seed-genomes', '2', '-o', str(scratch/'out.zst'),
                                         *extra, *map(str, inputs)], capture_output=True, text=True, timeout=15)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr.lower())
                self.assertFalse((scratch/'out.zst').exists())

    def run_tool(self, directory, output_name, inputs, *extra):
        output = directory / output_name
        run = subprocess.run([str(BIN), '--ggcat', str(GGCAT), '-j', '3', '-m', '0.1',
                              '-o', str(output), *extra, *map(str, inputs)],
                             capture_output=True, text=True, timeout=90)
        self.assertEqual(run.returncode, 0, run.stdout[-4000:] + run.stderr[-4000:])
        rows = zstd.decompress(output.read_bytes()).decode().splitlines()[2:]
        return rows, run

    @staticmethod
    def write_genomes(directory, genomes, prefix):
        paths = []
        for g, records in enumerate(genomes):
            path = directory / f'{prefix}-{g}.fa.zst'
            text = ''.join(f'>c{c}\n{s}\n' for c, s in enumerate(records))
            path.write_bytes(zstd.compress(text.encode()))
            paths.append(path)
        return paths

    def test_thresholds_include_missing_reference_nested_support_and_reject_repeats(self):
        rng = random.Random(5091)
        a, b, c, d, nested = [''.join(rng.choices('ACGT', k=n)) for n in (240, 210, 180, 150, 200)]
        genomes = [[a, b, c, c, nested], [a, b, c, d, nested],
                   [rc(a), c, d, rc(nested[:120])], [d]]
        with tempfile.TemporaryDirectory(prefix='mumento2-threshold-') as scratch:
            scratch = Path(scratch)
            inputs = self.write_genomes(scratch, genomes, 'g')
            expected = {
                '240\t0,0,0,.\t+,+,-,.\t0,0,0,.',
                '210\t0,0,.,.\t+,+,.,.\t1,1,.,.',
                '150\t.,0,0,0\t.,+,+,+\t.,3,2,0',
                '200\t0,0,.,.\t+,+,.,.\t4,4,.,.',
                '120\t0,0,0,.\t+,+,-,.\t4,4,3,.',
            }
            for k, minimum in ((21, 2), (21, 3), (21, 4), (20, 2)):
                rows, _ = self.run_tool(scratch, f'k{k}-N{minimum}.zst', inputs, '-k', str(k), '-N', str(minimum))
                wanted = {r for r in expected if sum(p != '.' for p in r.split('\t')[1].split(',')) >= minimum}
                self.assertEqual(set(rows), wanted)

    def test_incremental_threshold_matches_fresh_build_and_new_only_matches(self):
        rng = random.Random(2751)
        core, private, below, novel = [''.join(rng.choices('ACGT', k=n)) for n in (240, 160, 180, 190)]
        old = [[core, private], [rc(core), below]]
        new = [[core[:80], core[70:160], core[150:], private, novel], [rc(novel), below]]
        with tempfile.TemporaryDirectory(prefix='mumento2-update-') as scratch:
            scratch = Path(scratch)
            old_inputs = self.write_genomes(scratch, old, 'old')
            state = scratch / 'state'
            rows, _ = self.run_tool(scratch, 'initial.zst', old_inputs, '-k', '21', '-N', '2', '--state', str(state))
            self.assertEqual(rows, ['240\t0,0\t+,-\t0,0'])
            for path in old_inputs:
                path.unlink()
            new_inputs = self.write_genomes(scratch, new, 'new')
            updated, run = self.run_tool(scratch, 'updated.zst', new_inputs, '--update', str(state))
            self.assertIn('2 new source reads', run.stderr)
            self.assertIn('2 cached genomes', run.stderr)
            inputs = self.write_genomes(scratch, old + new, 'fresh')
            fresh, _ = self.run_tool(scratch, 'fresh.zst', inputs, '-k', '21', '-N', '2')
            self.assertEqual(set(updated), set(fresh))
            self.assertTrue({160, 180, 190, 240} <= {int(r.split('\t')[0]) for r in updated})
            self.assertIn('190\t.,.,0,0\t.,.,+,-\t.,.,4,0', updated)

    def test_incremental_all_genomes_splits_old_mum_and_inherits_options(self):
        rng = random.Random(1536)
        core = ''.join(rng.choices('ACGT', k=240))
        old = [[core], [core]]
        new = [[core[:80], core[70:160], core[150:]]]
        with tempfile.TemporaryDirectory(prefix='mumento2-update-all-') as scratch:
            scratch = Path(scratch)
            inputs = self.write_genomes(scratch, old, 'old')
            state = scratch / 'state'
            self.run_tool(scratch, 'initial.zst', inputs, '-k', '21', '--forward-only', '--state', str(state))
            for path in inputs:
                path.unlink()
            inputs = self.write_genomes(scratch, new, 'new')
            rows, _ = self.run_tool(scratch, 'updated.zst', inputs, '--update', str(state))
            self.assertEqual(set(rows), {'80\t0,0,0\t+,+,+\t0,0,0',
                                         '90\t70,70,0\t+,+,+\t0,0,1',
                                         '90\t150,150,0\t+,+,+\t0,0,2'})

    def test_long_proto_recovers_internal_fragments_at_contig_breaks(self):
        rng = random.Random(9182)
        core = ''.join(rng.choices('ACGT', k=1000))
        # All graph 31-mers survive, but no complete 1,000-base occurrence
        # exists in genome 1. Its middle fragment contains neither proto end.
        pieces = [core[:400], core[370:730], core[700:]]
        with tempfile.TemporaryDirectory(prefix='mumento2-fragments-') as scratch:
            scratch = Path(scratch)
            inputs = []
            for i, records in enumerate([[core], [rc(s) for s in pieces]]):
                path = scratch / f'g{i}.fa.zst'
                path.write_bytes(zstd.compress(''.join(f'>c{c}\n{s}\n' for c, s in enumerate(records)).encode()))
                inputs.append(str(path))
            output = scratch / 'mums.tsv.zst'
            protos = scratch / 'protos.fa.zst'
            run = subprocess.run([str(BIN), '--ggcat', str(GGCAT), '-k', '31', '-j', '4',
                                  '-m', '0.1', '--proto-output', str(protos), '-o', str(output), *inputs],
                                 capture_output=True, text=True, timeout=90)
            self.assertEqual(run.returncode, 0, run.stdout[-4000:] + run.stderr[-4000:])
            proto_text = zstd.decompress(protos.read_bytes()).decode()
            self.assertEqual([len(s) for s in proto_text.splitlines() if not s.startswith('>')], [1000])
            rows = zstd.decompress(output.read_bytes()).decode().splitlines()[2:]
            self.assertEqual(rows, ['400\t0,0\t+,-\t0,0', '360\t370,0\t+,-\t0,1', '300\t700,0\t+,-\t0,2'])

    def test_mum_with_no_unique_kmer_is_retained(self):
        rng = random.Random(20261015)
        block = ''.join(rng.choices('ACGT', k=80))
        k = 21
        # Every k-mer has its own extra contig, so no k-mer is single-copy.
        # The complete main sequence still occurs exactly once per genome.
        decoys = [block[i:i + k] for i in range(len(block) - k + 1)]
        with tempfile.TemporaryDirectory(prefix='mumento2-repeated-kmers-') as scratch:
            scratch = Path(scratch)
            inputs = []
            for genome in range(2):
                records = [block if genome == 0 else rc(block), *decoys]
                text = ''.join(f'>c{i}\n{s}\n' for i, s in enumerate(records)).encode()
                path = scratch / f'g{genome}.fasta.zst'
                path.write_bytes(zstd.compress(text))
                inputs.append(str(path))
            output = scratch / 'mums.tsv.zst'
            run = subprocess.run([str(BIN), '--ggcat', str(GGCAT), '-k', str(k),
                                  '-j', '4', '-m', '0.1', '-o', str(output), *inputs],
                                 capture_output=True, text=True, timeout=90)
            self.assertEqual(run.returncode, 0, run.stdout[-4000:] + run.stderr[-4000:])
            rows = zstd.decompress(output.read_bytes()).decode().splitlines()[2:]
            self.assertEqual(rows, ['80\t0,0\t+,-\t0,0'])

    def test_original_workflow_with_rearranged_and_reversed_contigs(self):
        rng = random.Random(425)
        blocks = ["".join(rng.choices("ACGT", k=length)) for length in (200, 230, 260)]
        genomes = [blocks, [rc(blocks[2]), blocks[0], blocks[1]], [blocks[1], rc(blocks[0]), blocks[2]]]
        with tempfile.TemporaryDirectory(prefix="mumento2-test-") as scratch:
            scratch = Path(scratch)
            inputs = []
            for g, records in enumerate(genomes):
                text = "".join(f">g{g}c{c}\n{record}\n" for c, record in enumerate(records)).encode()
                suffix, data = [("gz", gzip.compress(text)), ("xz", lzma.compress(text)),
                                ("zst", zstd.compress(text))][g]
                path = scratch / f"genome-{g}.fasta.{suffix}"
                path.write_bytes(data)
                inputs.append(str(path))
            results = []
            for threads, chunk, forward in [(1, 17, False), (4, 71, False), (4, 19, True)]:
                output = scratch / f"result-{threads}-{forward}.tsv.zst"
                command = [str(BIN), "--ggcat", str(GGCAT), "-k", "21", "-j", str(threads),
                           "--chunk-bases", str(chunk), "-m", "0.1", "-o", str(output), *inputs]
                if forward:
                    command.append("--forward-only")
                    command += ["--keep-temp", "--temp-dir", str(scratch / "retained")]
                run = subprocess.run(command, capture_output=True, text=True, timeout=90)
                self.assertEqual(run.returncode, 0, run.stdout[-4000:] + run.stderr[-4000:])
                text = zstd.decompress(output.read_bytes()).decode()
                lines = [line for line in text.splitlines() if not line.startswith("#")]
                self.assertEqual(lines[0], "length\tstarts\tstrands\tcontigs")
                rows = lines[1:]
                self.assertEqual(len(rows), 1 if forward else 3, text)
                observed = set()
                for row in rows:
                    length, starts, strands, contigs = row.split("\t")
                    length = int(length)
                    starts = list(map(int, starts.split(",")))
                    contigs = list(map(int, contigs.split(",")))
                    strands = strands.split(",")
                    self.assertEqual(starts, [0, 0, 0])
                    sequences = [genomes[g][c][s:s + length] for g, (c, s) in enumerate(zip(contigs, starts))]
                    sequences = [rc(s) if strand == "-" else s for s, strand in zip(sequences, strands)]
                    self.assertEqual(sequences, [sequences[0]] * 3)
                    observed.add(sequences[0])
                self.assertEqual(observed, {blocks[1]} if forward else set(blocks))
                mapping = Path(str(output) + ".contigs.tsv.zst")
                self.assertIn("genome\tinput\tcontig", zstd.decompress(mapping.read_bytes()).decode())
                if not forward:
                    results.append(sorted(rows))
                else:
                    # Includes native bucket containers, colormaps, indirect
                    # sequence storage, checkpoints, and proto-MUM FASTA.
                    for path in (scratch / "retained").rglob("*"):
                        if path.is_file():
                            data = path.read_bytes()
                            self.assertIn(data[:4], [b"\x28\xb5\x2f\xfd", b"\x50\x2a\x4d\x18"], str(path))
                            zstd.decompress(data)
            self.assertEqual(*results)


if __name__ == "__main__":
    unittest.main()
