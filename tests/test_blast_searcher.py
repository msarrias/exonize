from unittest.mock import Mock
from exonize.blast_searcher import BLASTsearcher
import portion as P
import pytest


blast_engine = BLASTsearcher(
    data_container=Mock(),
    sleep_max_seconds=40,
    self_hit_threshold=0.5,
    min_exon_length=20,
    cds_overlapping_threshold=0.8,
    evalue_threshold=1e-5,
    debug_mode=False,
)


blast_engine.data_container.gene_hierarchy_dictionary = dict(
    gene_1=dict(
        coordinates=P.open(1, 10),
        chrom='1',
        strand='+',
        mRNAs=dict(
            transcript_1=dict(
                coordinate=P.open(0, 127),
                strand='+',
                structure=[
                    dict(
                        id='CDS1_t1',
                        coordinate=P.open(0, 127),
                        frame=0,
                        type='CDS'
                    ),
                    dict(
                        id='CDS2_t1',
                        coordinate=P.open(4545, 4682),
                        frame=2,
                        type='CDS'
                    ),
                    dict(
                        id='CDS3_t1',
                        coordinate=P.open(6460, 6589),
                        frame=0,
                        type='CDS'
                    ),
                    dict(
                        id='CDS4_t1',
                        coordinate=P.open(7311, 7442),
                        frame=0,
                        type='CDS'
                    )
                ]
            ),
            transcript_2=dict(
                strand='+',
                structure=[
                    dict(
                        id='CDS1_t2',
                        coordinate=P.open(0, 127),
                        frame=0,
                        type='CDS'
                    ),
                    dict(
                        id='CDS2_t2',
                        coordinate=P.open(6460, 6589),
                        frame=2,
                        type='CDS'
                    ),
                    dict(
                        id='CDS3_t2',
                        coordinate=P.open(7311, 7442),
                        frame=2,
                        type='CDS'
                    )
                ]
            )
        )
    )
)


def test_get_overlap_percentage():
    # test no overlap
    assert blast_engine.get_overlap_percentage(
        intv_i=P.open(0, 1),
        intv_j=P.open(10, 100)
    ) == 0
    # test full overlap
    assert blast_engine.get_overlap_percentage(
        intv_i=P.open(10, 100),
        intv_j=P.open(10, 100)
    ) == 1
    # test partial overlap
    assert blast_engine.get_overlap_percentage(
        intv_i=P.open(10, 100),
        intv_j=P.open(15, 85)
    ) == (85 - 15) / (85 - 15)

    assert blast_engine.get_overlap_percentage(
        intv_i=P.open(15, 85),
        intv_j=P.open(10, 100)
    ) == (85 - 15) / (100 - 10)
    pass


def test_get_single_candidate_cds_coordinate():
    # interval i is contained in interval j
    assert blast_engine.get_single_candidate_cds_coordinate(
        intv_i=P.open(10, 100),
        intv_j=P.open(15, 85)
    ) == (
        P.open(15, 85)
    )
    # interval i overlaps interval j on the left
    assert blast_engine.get_single_candidate_cds_coordinate(
        intv_i=P.open(10, 50),
        intv_j=P.open(15, 55)
    ) == (
        P.open(10, 50)
    )
    # interval i overlaps interval j on the right
    assert blast_engine.get_single_candidate_cds_coordinate(
        intv_i=P.open(10, 55),
        intv_j=P.open(45, 105)
    ) == (
        P.open(10, 55)
    )


def test_compute_identity():
    assert blast_engine.compute_identity(
        sequence_i="ACGT",
        sequence_j="ACGT"
    ) == 1.0
    assert blast_engine.compute_identity(
        sequence_i="ACGT",
        sequence_j="AC-T"
    ) == 0.75
    assert blast_engine.compute_identity(
        sequence_i="ACGT",
        sequence_j="ATAT"
    ) == 0.5
    assert blast_engine.compute_identity(
        sequence_i="ACGT",
        sequence_j="TGCA"
    ) == 0.0
    with pytest.raises(ValueError):
        blast_engine.compute_identity(
            sequence_i="ACGT",
            sequence_j="ACG"
        )


def test_reformat_tblastx_frame_strand():
    assert blast_engine.reformat_tblastx_frame_strand(frame=1) == (0, '+')
    assert blast_engine.reformat_tblastx_frame_strand(frame=-1) == (0, '-')


def test_reverse_sequence_bool():
    assert blast_engine.reverse_sequence_bool(strand="+") is False
    assert blast_engine.reverse_sequence_bool(strand="-") is True


def test_get_first_overlapping_intervals():
    test_a = [
        P.open(0, 100),
        P.open(180, 300),
    ]
    res_a = (None, None)
    assert blast_engine.get_first_overlapping_intervals(
        sorted_intervals=test_a
    ) == res_a

    test_b = [
        P.open(180, 300),
        P.open(200, 300),
        P.open(250, 900),
    ]
    res_b = (
        P.open(180, 300),
        P.open(200, 300)
    )
    assert blast_engine.get_first_overlapping_intervals(
        sorted_intervals=test_b
    ) == res_b


def test_resolve_overlaps_between_coordinates():
    blast_engine.cds_overlapping_threshold = 0.8
    test = [
        P.open(0, 100),
        P.open(180, 300),
        P.open(200, 300),
        P.open(600, 900),
        P.open(700, 2000)
    ]
    res_a = [
        P.open(0, 100),
        P.open(200, 300),
        P.open(600, 900),
        P.open(700, 2000)
    ]
    assert blast_engine.resolve_overlaps_between_coordinates(
        sorted_cds_coordinates=test
    ) == res_a

    blast_engine.cds_overlapping_threshold = 0.3
    res_b = [
        P.open(0, 100),
        P.open(200, 300),
        P.open(600, 900),
        P.open(700, 2000)
    ]
    assert blast_engine.resolve_overlaps_between_coordinates(
        sorted_cds_coordinates=test
    ) == res_b

    blast_engine.cds_overlapping_threshold = 0.001
    res_c = [
        P.open(0, 100),
        P.open(200, 300),
        P.open(600, 900)
    ]
    assert blast_engine.resolve_overlaps_between_coordinates(
        sorted_cds_coordinates=test
    ) == res_c


def test_get_candidate_cds_coordinates():
    res_a_i = {
        P.open(0, 127): '0',
        P.open(4545, 4682): '2',
        P.open(6460, 6589): '0_2',
        P.open(7311, 7442): '0_2'
    }

    res_a_ii = [
        P.open(0, 127),
        P.open(4545, 4682),
        P.open(6460, 6589),
        P.open(7311, 7442)
    ]
    blast_res_a = blast_engine.get_candidate_cds_coordinates('gene_1')
    assert blast_res_a['cds_frame_dict'] == res_a_i
    assert blast_res_a['candidates_cds_coordinates'] == res_a_ii


def test_fetch_dna_sequence():
    blast_engine.data_container.genome_dictionary = {
        "chr1": "ATGC" * 100  # example sequence
    }

    sequence = blast_engine.fetch_dna_sequence(
        chromosome="chr1",
        annotation_start=0,
        annotation_end=8,
        trim_start=0,
        trim_end=4,
        strand="+"
    )
    assert sequence == "ATGC"
    sequence = blast_engine.fetch_dna_sequence(
        chromosome="chr1",
        annotation_start=0,
        annotation_end=8,
        trim_start=0,
        trim_end=4,
        strand="-"
    )
    assert sequence == "GCAT"  # Reverse complement of "ATGC"
    sequence = blast_engine.fetch_dna_sequence(
        chromosome="chr1",
        annotation_start=0,
        annotation_end=12,
        trim_start=2,
        trim_end=10,
        strand="+"
    )
    assert sequence == "GCATGCAT"  # Trimmed sequence
