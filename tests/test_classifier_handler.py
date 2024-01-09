from unittest.mock import Mock
from exonize.classifier_handler import ClassifierHandler

classifier_handler = ClassifierHandler(
    blast_engine=Mock(),
    cds_overlapping_threshold=0.8,
)


def test_find_overlapping_annotations():
    pass


def test_identify_query():
    pass


def test_target_out_of_mrna():
    pass


def test_indetify_full_target():
    pass


def test_filter_structure_by_interval_and_type():
    pass


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
