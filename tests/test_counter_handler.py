from exonize.blast_searcher import BLASTsearcher
from exonize.sqlite_handler import SqliteHandler
from exonize.data_preprocessor import DataPreprocessor
from exonize.counter_handler import CounterHandler

database_interface = SqliteHandler(
    results_database_path='',
    timeout_database=30,
)

data_container = DataPreprocessor(
            logger_obj=None,
            database_interface=database_interface,
            working_directory='',
            gff_file_path='',
            specie_identifier='test',
            genome_file_path='',
            genome_pickled_file_path='',
            debug_mode=False,
            evalue_threshold=1e-5,
)

blast_engine = BLASTsearcher(
    data_container=data_container,
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
