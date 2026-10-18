#!/usr/bin/env python3
"""Reproducible, reviewed scenario catalogue; never runs production tools.

Regenerate with python3 tests/generate_cases.py, or verify with --check.
Every catalogue entry becomes one separately named Rust and Python test.
"""
import argparse
import json
from pathlib import Path
import random
from oracle import expected_matches, expected_identical_contigs, rc

ROOT = Path(__file__).resolve().parents[1]


def build_cases():
    cases = []
    rng = random.Random(701031)
    def dna(n):
        return ''.join(rng.choices('ACGT', k=n))
    a, b, c, d = [dna(n) for n in (80, 65, 57, 43)]
    p = dna(12)
    palindrome = p + rc(p)
    def add(name, genomes, *, k=7, minimum=None, forward=False, description=None, **kwargs):
        case = dict(name=name, kind='discovery', description=description or name.replace('_', ' '),
                    genomes=genomes, k=k, minimum=minimum, forward_only=forward,
                    threads=3, chunk_bases=17)
        case.update(kwargs)
        case['expected'] = (expected_identical_contigs(genomes, k, forward) if case.get('oracle_mode')=='identical_contigs' else expected_matches(genomes, k, minimum, forward))
        cases.append(case)
    # Presence/support geometry, including support changes on both sides.
    add('identical_two', [[a], [a]])
    add('identical_five', [[a]] * 5)
    add('unrelated_documents', [[a], [b], [c]])
    add('one_document_unique_runs', [[a, b]], minimum=1)
    add('private_runs_n_one', [[a], [b], [c]], minimum=1)
    add('pair_absent_first', [[b], [a], [a]], minimum=2)
    add('pair_absent_last', [[a], [a], [b]], minimum=2)
    add('pair_absent_first_two', [[b], [c], [a], [a]], minimum=2)
    add('support_nonadjacent', [[a], [b], [a], [c], [a]], minimum=3)
    add('exact_three_of_five', [[b], [a], [a], [c], [a]], minimum=3)
    add('above_three_of_five', [[a], [a], [a], [b], [a]], minimum=3)
    add('below_threshold_two_of_five', [[a], [a], [b], [c], [d]], minimum=3)
    add('nested_right_support', [[a], [a], [a[:51]], [b]], minimum=2)
    add('nested_left_support', [[a], [a], [a[29:]], [b]], minimum=2)
    add('nested_both_support', [[a], [a[10:70]], [a[25:55]], [b]], minimum=2)
    add('disjoint_support_sets', [[a], [a], [b], [b]], minimum=2)
    add('incomparable_support_sets', [[a], [a, b], [b], [c]], minimum=2)
    add('common_and_private_contigs', [[a, b], [a, c], [a, d]], minimum=2)
    add('alternating_six_support', [[a], [b], [a], [c], [a], [d]], minimum=3)
    add('late_pair_empty_first', [[], [b], [a], [a]], minimum=2)
    # Strand orientation, physical uniqueness of palindromes and RC repeats.
    add('one_reverse_document', [[a], [rc(a)], [a]])
    add('first_document_reverse', [[rc(a)], [a], [a]])
    add('alternating_reverse_five', [[a], [rc(a)], [a], [rc(a)], [a]], minimum=3)
    add('forward_only_opposite_strands', [[a], [rc(a)]], forward=True)
    add('forward_only_mixed_support', [[a], [rc(a)], [a]], minimum=2, forward=True)
    add('forward_only_rc_is_distinct', [[a, rc(a)], [a]], forward=True)
    add('both_strands_same_document_repeat', [[a, rc(a)], [a]])
    add('full_palindrome_single_physical_hit', [[palindrome], [palindrome]], k=5)
    add('palindrome_seed_asymmetric_flanks', [['AACGTACGGA'], ['TCCGTACGTT']], k=4)
    add('reverse_hit_with_contig_offsets', [[b, 'NN' + a], [c, 'NNNN' + rc(a)]])
    add('reverse_lowercase_input', [[a.lower()], [rc(a).lower()]])
    # Exact extension must stop at actual mismatches, ends and rearrangements.
    mutate = lambda s, i: s[:i] + ('A' if s[i] != 'A' else 'C') + s[i+1:]
    add('snp_first_base', [[a], [mutate(a, 0)]])
    add('snp_last_base', [[a], [mutate(a, len(a)-1)]])
    add('snp_middle_base', [[a], [mutate(a, 40)]])
    add('two_internal_snps', [[a], [mutate(mutate(a, 25), 54)]])
    add('single_base_insertion', [[a], [a[:40] + 'T' + a[40:]]])
    add('single_base_deletion', [[a], [a[:40] + a[41:]]])
    add('internal_inversion', [[a], [a[:20] + rc(a[20:60]) + a[60:]]])
    add('translocated_contigs', [[a, b, c], [c, a, b]])
    add('overlapping_contig_fragments', [[a], [a[:32], a[26:60], a[54:]]], patterns=[a])
    add('adjacent_contigs_never_bridge', [[a], [a[:40], a[40:]]], patterns=[a])
    add('ambiguity_internal_n', [[a], [a[:40] + 'N' + a[41:]]])
    add('ambiguity_iupac_internal', [[a], [a[:20] + 'RYWSKMBDHV' + a[30:]]])
    add('ambiguity_lowercase_n', [[a], [a[:40].lower() + 'nn' + a[40:].lower()]])
    add('leading_ambiguity_offsets', [['NNNN' + a], ['RYN' + a]])
    add('trailing_ambiguity', [[a + 'NNRY'], [a + 'NN']])
    add('mixed_case_dna', [[a.lower()], [''.join(x.lower() if i % 2 else x for i,x in enumerate(a))]])
    # Whole-string repeats, overlapping repeats and nonunique k-mer seeds.
    add('full_repeat_one_genome', [[a, a], [a]])
    add('repeat_in_additional_support_genome', [[a], [a], [a, a], [b]], minimum=2)
    add('repeated_prefix_decoy', [[a, a[:31]], [a]])
    add('repeated_suffix_decoy', [[a, a[-31:]], [a]])
    add('repeated_interior_decoy', [[a, a[21:52]], [a]])
    add('same_seed_different_full_context', [['ACGATCGAANACGTTAGGC'], ['ACGTTAGGCNACGATCGAA']], k=3)
    add('every_kmer_has_decoy_contig', [[a, *[a[i:i+7] for i in range(len(a)-6)]], [rc(a), *[a[i:i+7] for i in range(len(a)-6)]]], seed_mode='kmers')
    add('overlapping_homopolymer', [['AAAAAAAAAAA'], ['AAAAAAAAA']], k=3)
    add('tandem_repeat_unique_flanks', [['GACACACACT'], ['GACACACACT']], k=3)
    add('duplicated_contig_names_independent_positions', [[a, a], [a]], header='same')
    add('palindrome_two_physical_positions', [[palindrome, palindrome], [palindrome]], k=4)
    add('reverse_complement_copy_is_repeat', [[a + 'NN' + rc(a)], [a]])
    add('both_flanks_disambiguate_repeated_seed', [['TACGATNNTACGACNNGACGAT'], ['TACGATNNTACGAGNNCACGAT']], k=3, seed_mode='kmers')
    add('duplicate_proto_seeds_deduplicate_output', [[a], [a], [rc(a)]], seed_mode='duplicates')
    add('reverse_proto_seeds_keep_physical_coordinates', [[a], [rc(a)]], seed_mode='reverse')
    # Packed word widths, exact minimum-length edges and SIMD chunk seams.
    for k in (3, 4, 15, 29, 31, 32, 33, 63, 64, 65, 67, 129):
        s = dna(k + 19)
        add(f'k_{k}_packed_word_boundary', [[s], [rc(s)]], k=k,
            description=f'Canonical k={k} seeds spanning packed words and reverse coordinates')
    s = dna(31)
    add('match_length_k_minus_one', [[s[:30]], [s[:30]]], k=31)
    add('match_length_exactly_k', [[s], [s]], k=31)
    add('match_length_k_plus_one', [[s + 'A'], [s + 'A']], k=31)
    for chunk in (1, 6, 7, 8, 29, 64, 65):
        s = dna(91)
        add(f'chunk_{chunk}_exact_extension', [['NN' + s], [rc(s) + 'NN']], chunk_bases=chunk,
            description=f'Seed/extension crossing chunk boundary {chunk}, including reverse hit')
    add('many_short_contigs_before_valid_hit', [['AC', 'T', 'NN', a], ['T', 'A', a]])
    add('empty_first_document', [[], [a], [a]], minimum=2)
    add('empty_last_document_all_required', [[a], [a], []])
    add('all_documents_empty', [[], []])
    add('ambiguity_only_documents', [['NNRYWSKMBDHV'], ['NNNN']])
    add('all_contigs_shorter_than_k', [['ACG', 'TACG'], ['ACGT', 'CG']], k=7)
    add('empty_record_preserves_contig_id', [['', a], [a, '']])
    add('minimizer_one', [[a], [rc(a)]], m=1)
    add('minimizer_equal_k', [[a], [rc(a)]], m=7)
    add('one_worker_many_support_sets', [[a,b], [a], [b], [a[:40]]], minimum=2, threads=1)
    add('more_workers_than_documents', [[a], [rc(a)]], threads=8)
    # Actual parsing/compression variants, not variations of random seeds alone.
    for fmt in ('gzip', 'xz', 'zstd', 'mixed', 'gzip_members', 'xz_members', 'zstd_frames'):
        add(f'fasta_{fmt}', [[a[:11], b], [b, a[:11]]], compression=fmt)
    for fmt in ('plain', 'gzip', 'xz', 'zstd'):
        add(f'fastq_{fmt}', [[a,b], [rc(a),b]], fastx='fastq', compression=fmt)
    for wrap in (1, 7, 63, 64, 65):
        add(f'fasta_wrapped_width_{wrap}', [[a], [rc(a)]], wrap=wrap)
    add('fasta_crlf', [[a,b], [rc(a),b]], crlf=True)
    add('fasta_wrapped_crlf', [[a,b], [rc(a),b]], crlf=True, wrap=7)
    add('fasta_crlf_packed_block_boundary', [[a,b], [rc(a),b]], crlf=True, wrap=59)
    add('fasta_gzip_crlf', [[a,b], [rc(a),b]], crlf=True, compression='gzip', wrap=7)
    add('fastq_crlf', [[a,b], [rc(a),b]], crlf=True, fastx='fastq')
    add('fastq_crlf_dna_newline_crosses_simd_block', [[a[:58],b], [rc(a[:58]),b]], crlf=True, fastx='fastq')
    add('fasta_no_final_newline', [[a,b], [a,b]], final_newline=False)
    add('fasta_header_descriptions', [[a,b], [rc(a),b]], header='description')
    add('fasta_unicode_header', [[a,b], [a,rc(b)]], header='unicode')
    add('fastq_quality_contains_header_markers', [[a,b], [rc(a),b]], fastx='fastq', quality='@+')
    # Real reader refills and the >=1024-byte FASTQ quality-skip path.
    long_record = dna(131137)
    add('fasta_reader_buffer_refill', [[long_record], [rc(long_record)]], k=31, oracle_mode='identical_contigs')
    add('fasta_wrapped_reader_buffer_refill', [[long_record], [rc(long_record)]], k=31, wrap=71, compression='mixed', oracle_mode='identical_contigs')
    add('fastq_long_quality_skip_refreshes_lexer', [[long_record[:4099],b], [rc(long_record[:4099]),b]], k=31, fastx='fastq', quality='@+', oracle_mode='identical_contigs')
    add('fastq_gzip_quality_skip_buffer_refill', [[long_record,b], [rc(long_record),b]], k=31, fastx='fastq', compression='gzip', oracle_mode='identical_contigs')
    add('fasta_header_larger_than_reader_buffer', [[a,b],[rc(a),b]], header_padding=131100)
    add('fastq_header_larger_than_reader_buffer', [[a,b],[rc(a),b]], fastx='fastq', header_padding=131100)
    # Incremental operations: every step is checked against the oracle; final
    # state is also compared to a separate full rebuild with real new files.
    def update(name, old, additions, *, minimum=None, k=7, forward=False, **extra):
        case = dict(name=name, kind='update', description=name.replace('_', ' '), genomes=old,
                    k=k, minimum=minimum, forward_only=forward, threads=3, chunk_bases=11, **extra)
        case['expected'] = expected_matches(old, k, minimum, forward)
        all_genomes, steps, threshold = list(old), [], minimum
        for addition in additions:
            step = addition if isinstance(addition, dict) else {'genomes': addition}
            step = dict(step)
            if step.get('minimum') is not None:
                threshold = step['minimum']
            all_genomes.extend(step['genomes'])
            step['expected'] = expected_matches(all_genomes, k, threshold, forward)
            steps.append(step)
        case['updates'] = steps
        cases.append(case)
    update('update_all_preserves_full_match', [[a],[a]], [[[a]]])
    update('update_all_new_reverse_genome', [[a],[rc(a)]], [[[rc(a)]]])
    update('update_all_prefix_survives', [[a],[a]], [[[a[:43]]]])
    update('update_all_suffix_survives', [[a],[a]], [[[a[30:]]]])
    update('update_all_interior_survives', [[a],[a]], [[[a[20:61]]]])
    update('update_all_splits_at_new_contigs', [[a],[a]], [[[a[:32],a[26:60],a[54:]]]])
    update('update_all_splits_at_new_ambiguity', [[a],[a]], [[[a[:40]+'N'+a[41:]]]])
    update('update_all_splits_at_new_snp', [[a],[a]], [[[mutate(a,40)]]])
    update('update_all_repeated_new_genome_invalidates', [[a],[a]], [[[a,a]]])
    update('update_all_empty_new_genome_removes_matches', [[a],[a]], [[[]]])
    update('update_fixed_threshold_keeps_old_only_match', [[a],[a]], [[[b]]], minimum=2)
    update('update_below_threshold_old_region_qualifies', [[a,b],[a]], [[[b]]], minimum=2)
    update('update_new_only_shared_region', [[a],[a]], [[[b],[rc(b)]]], minimum=2)
    update('update_new_only_region_second_batch', [[a],[a]], [[[b]], [[rc(b)]]], minimum=2)
    update('update_new_contig_order', [[a,b],[b,a]], [[[rc(b),rc(a)]]])
    update('update_n_one_private_new_regions', [[a],[b]], [[[c]]], minimum=1)
    update('update_threshold_increases', [[a,b],[a,b]], [{'genomes':[[a]], 'minimum':3}], minimum=2)
    update('update_threshold_decreases', [[a,b],[a],[c]], [{'genomes':[[b]], 'minimum':2}])
    update('update_sequential_three_additions', [[a],[a]], [[[a[:65]]], [[a[10:]]], [[a[20:60]]]])
    update('update_batch_multiple_documents', [[a],[rc(a)]], [[[a[:60]],[a[20:]]]])
    update('update_fork_preserves_source', [[a],[a]], [{'genomes':[[b]], 'fork':True}], minimum=2)
    update('update_forward_only_mode_inherited', [[a],[a]], [[[a,b]]], forward=True)
    update('update_locator_minimizer_inherited', [[a],[a]], [[[rc(a)]]], m=1)
    update('update_compressed_fastq_inputs', [[a],[rc(a)]], [[[a,b]], [[rc(b)]]], minimum=2, fastx='fastq', compression='mixed')
    update('update_empty_old_genome_late_pair', [[],[a]], [{'genomes':[[a]], 'minimum':2}])
    update('update_long_k_and_partial_range', [[a],[a]], [[[a[5:75]]]], k=63)
    update('update_rearrangement_with_no_novel_kmers', [[a,b],[a,b]], [[[b,a]]], minimum=2)
    update('update_duplicate_seed_growth_is_bounded', [[a],[a]], [[[a]], [[a]], [[rc(a)]]], minimum=2, stable_patterns=True)
    # Error paths must fail without changing an existing generation/output.
    errors = ('minimum_zero', 'minimum_exceeds_documents', 'threads_zero', 'chunk_zero',
              'minimizer_zero', 'minimizer_exceeds_k', 'minimizer_exceeds_word',
              'missing_input', 'duplicate_input', 'missing_state', 'existing_initial_state',
              'update_changed_k', 'update_changed_minimizer', 'update_changed_strand',
              'update_duplicate_document', 'unsupported_state_version', 'missing_current',
              'invalid_current_path', 'missing_state_manifest', 'invalid_state_json',
              'zero_state_genomes', 'invalid_proto_dna', 'invalid_packed_magic',
              'truncated_context_record', 'locked_state', 'existing_fork_destination',
              'missing_output_parent', 'missing_proto_archive')
    reasons = {
        'minimum_zero':'require 1 <= N', 'minimum_exceeds_documents':'require 1 <= N',
        'threads_zero':'threads and chunk-bases must be positive',
        'chunk_zero':'threads and chunk-bases must be positive',
        'minimizer_zero':'require k >= m', 'minimizer_exceeds_k':'require k >= m',
        'minimizer_exceeds_word':'require k >= m', 'duplicate_input':'same genome file',
        'existing_initial_state':'state already exists', 'update_changed_k':'updates must retain',
        'update_changed_minimizer':'updates must retain', 'update_changed_strand':'updates must retain',
        'update_duplicate_document':'same genome file', 'unsupported_state_version':'unsupported or invalid',
        'invalid_current_path':'invalid state generation name', 'invalid_state_json':'key must be a string',
        'zero_state_genomes':'unsupported or invalid', 'invalid_proto_dna':'proto-MUMs must contain only',
        'invalid_packed_magic':'invalid packed genome state', 'truncated_context_record':'invalid context record length',
        'locked_state':'locked', 'existing_fork_destination':'destination state already exists',
    }
    for error in errors:
        cases.append(dict(name='error_'+error, kind='error', error=error,
                          description=f'Reject {error.replace("_", " ")} and preserve committed data',
                          error_contains=reasons.get(error,'No such file'),
                          cli_error_contains='genome does not exist' if error=='missing_input' else reasons.get(error,'No such file'),
                          genomes=[[a],[a]], k=7, minimum=None, forward_only=False, threads=2, chunk_bases=11))
    assert len({c['name'] for c in cases}) == len(cases)
    assert sum(c['kind']=='discovery' for c in cases) >= 100
    return cases


def artifacts(cases):
    json_text = json.dumps(cases, indent=2, ensure_ascii=False) + '\n'
    rust = '// Generated by tests/generate_cases.py; each entry is a separate test.\n'
    rust += '\n'.join(f'corpus_test!({c["name"]});' for c in cases) + '\n'
    return {ROOT/'tests/corpus/cases.json':json_text,
            ROOT/'k2rmini/src/locate/stream/shared/corpus_cases.rs':rust}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args=parser.parse_args()
    cases=build_cases()
    for path, text in artifacts(cases).items():
        if args.check:
            if not path.exists() or path.read_text()!=text:
                raise SystemExit(f'Out of date: {path}')
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
    counts={kind:sum(c['kind']==kind for c in cases) for kind in ('discovery','update','error')}
    print(f'{len(cases)} named scenarios per layer: {counts}')

if __name__=='__main__': main()
