from exonize.blast_searcher import BLASTsearcher

blast_engine = BLASTsearcher(
    data_container=None,
    masking_percentage_threshold=0.8,
    sleep_max_seconds=40,
    self_hit_threshold=0.5,
    min_exon_length=20,
    cds_overlapping_threshold=0.8,
    evalue_threshold=1e-5,
    debug_mode=False,
)


def test_get_overlap_percentage():
    pass


def test_get_shorter_longer_intervals():
    pass


def test_compute_identity():
    pass


def test_reformat_tblastx_frame_strand():
    assert blast_engine.reformat_tblastx_frame_strand(frame=1) == '+', 0
    assert blast_engine.reformat_tblastx_frame_strand(frame=-1) == '-', 0


def test_reverse_sequence_bool():
    assert blast_engine.reverse_sequence_bool(strand="+") is False
    assert blast_engine.reverse_sequence_bool(strand="-") is True


def test_get_first_overlapping_intervals():
    pass


def test_resolve_overlaps_coords_list():
    pass


def test_get_candidate_cds_coordinates():
    pass


def test_align_cds():
    pass


def check_for_masking():
    pass


def test_fetch_dna_sequence():
    pass

