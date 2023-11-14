from exonize.exonize import *

test_object = Exonize('',
                      '',
                      'test',
                      min_exon_length=20,
                      cds_overlapping_threshold=0.8,
                      masking_perc_threshold=0.8,
                      self_hit_threshold=0.5)


def test_full_matches():
    # tuple format: (fragment_id, gene_id, query_start, query_end, target_start, target_end, type)
    matches_a = [(1, 'gene_test', 5, 10, 25, 30, 'INS_CDS'),
                 (2, 'gene_test', 20, 35, 0, 15, 'TRUNC'),
                 (3, 'gene_test', 4, 10, 23, 29, 'TRUNC')]
    test_object.cds_overlapping_threshold = 0.8
    assert test_object.assign_pair_ids(matches_a) == [(1, 'MATCH', 1), (1, 'RECIPROCAL', 2), (1, 'OVERLAPPING', 3)]
    matches_b = [(1, 'gene_test', 5, 10, 25, 30, 'FULL'),
                 (2, 'gene_test', 25, 30, 5, 10, 'FULL'),
                 (3, 'gene_test', 4, 10, 23, 29, 'TRUNC')]
    assert test_object.assign_pair_ids(matches_b) == [(1, 'MATCH', 1), (1, 'RECIPROCAL', 2), (1, 'OVERLAPPING', 3)]
    matches_c = [(1, 'gene_test', 5, 10, 25, 30, 'FULL'),
                 (2, 'gene_test', 25, 30, 5, 10, 'FULL'),
                 (3, 'gene_test', 4, 10, 50, 56, 'TRUNC')]
    assert test_object.assign_pair_ids(matches_c) == [(1, 'MATCH', 1), (1, 'RECIPROCAL', 2), (1, 'NON-RECIPROCAL', 3)]
