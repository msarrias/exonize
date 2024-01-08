from exonize.blast_searcher import BLASTsearcher
from exonize.classifier_handler import ClassifierHandler
from exonize.sqlite_handler import SqliteHandler
from exonize.data_preprocessor import DataPreprocessor

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

classifier_handler = ClassifierHandler(
    blast_engine=blast_engine,
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
