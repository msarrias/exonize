from unittest.mock import Mock
from exonize.blast_searcher import BLASTsearcher
from exonize.classifier_handler import ClassifierHandler
import portion as P

blast_engine = BLASTsearcher(
    data_container=Mock(),
    masking_percentage_threshold=0.8,
    sleep_max_seconds=40,
    self_hit_threshold=0.5,
    min_exon_length=20,
    cds_overlapping_threshold=0.8,
    evalue_threshold=1e-5,
    debug_mode=False,
)
classifier_handler = ClassifierHandler(
    blast_engine=blast_engine,
    cds_overlapping_threshold=0.8,
)


def test_recover_query_cds():
    transcript_dictionary = {
        "structure": [
            {"id": "CDS1", "coordinate": P.open(74, 123), "frame": "0", "type": "CDS"},
            {"id": "CDS2", "coordinate": P.open(50, 150), "frame": "0", "type": "CDS"},
            {"id": "CDS3", "coordinate": P.open(200, 300), "frame": "0", "type": "CDS"},
            {"id": "CDS4", "coordinate": P.open(76, 123), "frame": "0", "type": "CDS"}
        ]
    }
    cds_coordinate = P.open(76, 123)
    cds_annot = classifier_handler.recover_query_cds(
        transcript_dictionary=transcript_dictionary,
        query_coordinate=cds_coordinate
    )
    assert cds_annot == ('CDS4', P.open(76, 123), '0')


def test_find_overlapping_annotation():
    classifier_handler.cds_overlapping_threshold = 0.6

    # Example transcript dictionary and CDS coordinate
    transcript_dictionary = {
        "structure": [
            {"id": "CDS1", "coordinate": P.open(0, 100), "frame": "0", "type": "CDS"},
            {"id": "CDS2", "coordinate": P.open(48, 150), "frame": "0", "type": "CDS"},
            {"id": "CDS2", "coordinate": P.open(50, 150), "frame": "0", "type": "CDS"},
            {"id": "CDS3", "coordinate": P.open(200, 300), "frame": "0", "type": "CDS"}
        ]
    }
    cds_coordinate = P.open(50, 125)
    best_matching_cds = classifier_handler.find_overlapping_annotation(
        transcript_dictionary=transcript_dictionary,
        cds_coordinate=cds_coordinate
    )
    # Expected result
    expected_result = ("CDS2", P.open(50, 150), 0)
    assert best_matching_cds == expected_result


def test_indentify_insertion_target():
    pass


def test_indentify_truncation_target():
    pass


def test_identify_obligate_pair():
    pass


def test_identify_neither_pair():
    pass


def test_insert_full_length_duplication_tuple():
    pass


def test_get_interval_dictionary():
    pass
