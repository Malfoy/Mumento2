"""Independent, exhaustive multi-MUM oracle for small correctness fixtures.

No graph, seed, hash, alignment-range or production locator code is used.
Occurrence tables enumerate actual substrings and count physical positions.
"""
from collections import defaultdict

DNA = frozenset('ACGT')
TRANS = str.maketrans('ACGTacgt', 'TGCAtgca')


def rc(sequence):
    return sequence.translate(TRANS)[::-1]


def dna_runs(sequence):
    sequence = sequence.upper()
    begin = 0
    for end, base in enumerate(sequence + '!'):
        if base not in DNA:
            if begin < end:
                yield begin, sequence[begin:end]
            begin = end + 1


def expected_matches(genomes, k, minimum=None, forward_only=False):
    """Return canonical DNA and all unique physical (contig,start) positions.

    All genomes containing a string must be single-copy. Maximality requires
    that its entire support set cannot extend left or right together; a longer
    string losing a genome does not suppress the shorter string. Palindromes
    have one physical position with two possible extension orientations.
    """
    minimum = len(genomes) if minimum is None else minimum
    occurrences = defaultdict(lambda: defaultdict(set))
    normalized = [[s.upper() for s in genome] for genome in genomes]
    for g, genome in enumerate(normalized):
        for c, sequence in enumerate(genome):
            for origin, run in dna_runs(sequence):
                for start in range(len(run)):
                    for end in range(start + k, len(run) + 1):
                        pattern = run[start:end]
                        if not forward_only:
                            pattern = min(pattern, rc(pattern))
                        occurrences[pattern][g].add((c, origin + start))
    results = []
    for pattern, documents in occurrences.items():
        if len(documents) < minimum or any(len(p) != 1 for p in documents.values()):
            continue
        left_sets, right_sets = [], []
        reverse_pattern = rc(pattern)
        for g, positions in documents.items():
            c, start = next(iter(positions))
            sequence = normalized[g][c]
            physical = sequence[start:start + len(pattern)]
            orientations = [physical != pattern]
            if not forward_only and pattern == reverse_pattern:
                orientations = [False, True]
            left, right = set(), set()
            for reverse in orientations:
                li = start + len(pattern) if reverse else start - 1
                ri = start - 1 if reverse else start + len(pattern)
                if 0 <= li < len(sequence) and sequence[li] in DNA:
                    left.add(rc(sequence[li]) if reverse else sequence[li])
                if 0 <= ri < len(sequence) and sequence[ri] in DNA:
                    right.add(rc(sequence[ri]) if reverse else sequence[ri])
            left_sets.append(left)
            right_sets.append(right)
        if set.intersection(*left_sets) or set.intersection(*right_sets):
            continue
        hits = [list(next(iter(documents[g]))) if g in documents else None
                for g in range(len(genomes))]
        results.append({'sequence': pattern, 'hits': hits})
    return sorted(results, key=lambda m: (m['sequence'], repr(m['hits'])))


def match_keys(matches):
    return {(m['sequence'], tuple(None if p is None else tuple(p) for p in m['hits']))
            for m in matches}


def expected_identical_contigs(genomes, k, forward_only=False):
    """Linear-size oracle for buffer-refill fixtures with identical contig sets.

    Equal contig multisets (allowing RC) imply every unique proper substring
    extends in every genome along its containing contig. Only complete unique
    contigs can be MUMs. Verify that premise, then count their exact occurrences.
    This avoids quadratic substring enumeration on 128-KiB parser fixtures.
    """
    normalized = [[s.upper() for s in g] for g in genomes]
    canonical = lambda s: s if forward_only else min(s, rc(s))
    reference = sorted(canonical(s) for s in normalized[0])
    assert all(sorted(canonical(s) for s in g) == reference for g in normalized)
    assert all(set(s) <= DNA for g in normalized for s in g)
    matches = []
    for p in sorted(set(reference)):
        if len(p) < k:
            continue
        hits = []
        for genome in normalized:
            found = set()
            for c, s in enumerate(genome):
                for orientation in {p} if forward_only else {p, rc(p)}:
                    start = s.find(orientation)
                    while start >= 0:
                        found.add((c, start))
                        start = s.find(orientation, start + 1)
            hits.append(list(next(iter(found))) if len(found) == 1 else None)
        if all(p is not None for p in hits):
            matches.append({'sequence': p, 'hits': hits})
    return matches


def verify_output(test, text, genomes, k, minimum, forward_only, expected):
    """Check structure, coordinates, DNA, strands, uniqueness and completeness.

    Expects contig-relative 4-col rows with no header (mumento2 --contig-coords).
    """
    lines = [line for line in text.splitlines() if line and not line.startswith('#')]
    observed = []
    order = []
    for row in lines:
        fields = row.split('\t')
        test.assertEqual(len(fields), 4, row)
        length = int(fields[0])
        test.assertGreaterEqual(length, k, row)
        starts, strands, contigs = (s.split(',') for s in fields[1:])
        test.assertEqual([len(starts), len(strands), len(contigs)], [len(genomes)] * 3, row)
        hits, sequences = [], []
        for g, (p, strand, c) in enumerate(zip(starts, strands, contigs)):
            if p == '.':
                test.assertEqual((strand, c), ('.', '.'), row)
                hits.append(None)
                continue
            start, contig = int(p), int(c)
            test.assertGreaterEqual(start, 0, row)
            test.assertTrue(0 <= contig < len(genomes[g]), row)
            sequence = genomes[g][contig].upper()
            test.assertLessEqual(start + length, len(sequence), row)
            actual = sequence[start:start + length]
            test.assertTrue(set(actual) <= DNA, row)
            test.assertIn(strand, ('+', '-'), row)
            if forward_only:
                test.assertEqual(strand, '+', row)
            if not sequences:
                test.assertEqual(strand, '+', row)
            sequences.append(rc(actual) if strand == '-' else actual)
            hits.append([contig, start])
        test.assertGreaterEqual(len(sequences), minimum or len(genomes), row)
        test.assertTrue(all(s == sequences[0] for s in sequences), row)
        pattern = sequences[0] if forward_only else min(sequences[0], rc(sequences[0]))
        observed.append({'sequence': pattern, 'hits': hits})
        owner = next(i for i, p in enumerate(hits) if p is not None)
        contig, start = hits[owner]
        absolute = sum(len(s) for s in genomes[owner][:contig]) + start
        order.append((owner, absolute, length))
    test.assertEqual(order, sorted(order), 'output must have deterministic coordinate order')
    keys = match_keys(observed)
    test.assertEqual(len(keys), len(observed), 'duplicate MUM rows')
    test.assertEqual(keys, match_keys(expected))
    return keys
