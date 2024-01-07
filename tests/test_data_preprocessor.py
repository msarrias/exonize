from exonize.data_preprocessor import DataPreprocessor

data_container = DataPreprocessor(
            logger_obj=None,
            database_interface=None,
            working_directory=None,
            gff_file_path=None,
            specie_identifier='test',
            genome_file_path='',
            genome_pickled_file_path=None,
            debug_mode=False,
            evalue_threshold=1e-5,
)


def construct_mrna_sequences():
    pass


def test_check_for_overhangs():
    pass


def test_construct_peptide_sequences():
    pass

