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
