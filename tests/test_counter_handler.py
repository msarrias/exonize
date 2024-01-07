from exonize.blast_searcher import BLASTsearcher
from exonize.counter_handler import CounterHandler

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

counter_handler = CounterHandler(
    blast_engine=blast_engine,
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