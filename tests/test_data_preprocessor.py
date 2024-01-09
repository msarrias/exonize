from exonize.data_preprocessor import DataPreprocessor
from unittest.mock import Mock
import portion as P

data_container = DataPreprocessor(
            logger_obj=Mock(),
            database_interface=Mock(),
            working_directory='',
            gff_file_path='',
            specie_identifier='test',
            genome_file_path='',
            genome_pickled_file_path='',
            debug_mode=False,
            hard_masking=False,
            evalue_threshold=1e-5,
)

data_container.genome_dictionary = {
    '1': 'TAAAATCTAGACAGAAGCATTCTCAGAAACTTCTTTGTGCTGTATGTCCTCAATTAACAG',
    '2': 'AGTTGAACCTTTGTTTCGATACAGCATTTTGGAAACATTCCTTTAGTAGAATCTGCAAGT'
}

data_container.gene_hierarchy_dictionary = dict(
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


def test_construct_mrna_sequence():
    pass


def test_check_for_overhangs():
    pass


def test_construct_peptide_sequences():
    pass
