from unittest.mock import Mock
from exonize.counter_handler import CounterHandler

counter_handler = CounterHandler(
    blast_engine=Mock(),
    cds_overlapping_threshold=0.8,
)


def test_get_average_overlap_percentage():
    pass


def test_get_shorter_intv_overlapping_percentage():
    pass


def test_get_candidate_reference_dictionary():
    pass


def test_overlap_condition():
    pass


def test_get_overlapping_clusters():
    pass


def test_build_reference_dictionary():
    pass


def test_assign_event_ids():
    pass
