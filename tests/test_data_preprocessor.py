from exonize.data_preprocessor import DataPreprocessor
from exonize.sqlite_handler import SqliteHandler

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
            hard_masking=False,
            evalue_threshold=1e-5,
)


def construct_mrna_sequences():
    pass


def test_check_for_overhangs():
    pass


def test_construct_peptide_sequences():
    pass
