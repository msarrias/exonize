from exonize.data_preprocessor import DataPreprocessor
from unittest.mock import Mock
import pytest
from Bio.Seq import Seq
import portion as P
from pathlib import Path

data_container = DataPreprocessor(
    logger_obj=Mock(),
    database_interface=Mock(),
    working_directory=Path(''),
    gff_file_path=Path(''),
    output_prefix='test',
    genome_file_path=Path(''),
    self_hit_threshold=0.5,
    cds_overlapping_threshold=0.8,
    query_overlapping_threshold=0.9,
    min_exon_length=30,
    debug_mode=False,
    evalue_threshold=1e-5,
)

data_container.genome_dictionary = {
    "chr1": "ATGC" * 100  # Simulate a genome sequence for testing
}


def test_construct_mrna_sequence():
    cds_coordinates_list = [
        {"coordinate": P.open(0, 4)},
        {"coordinate": P.open(4, 8)}
    ]
    expected_sequence = "ATGCATGC"
    assert data_container.construct_mrna_sequence(
        chromosome="chr1",
        gene_strand="+",
        cds_coordinates_list=cds_coordinates_list
    ) == expected_sequence
    # Test case for negative strand
    cds_coordinates_list = [
        {"coordinate": P.open(0, 4)},
        {"coordinate": P.open(4, 8)}
    ]
    expected_sequence = Seq(
        data_container.genome_dictionary["chr1"][4:8] +
        data_container.genome_dictionary["chr1"][0:4]
    ).reverse_complement()
    assert data_container.construct_mrna_sequence(
        chromosome="chr1",
        gene_strand="-",
        cds_coordinates_list=cds_coordinates_list
    ) == expected_sequence


def test_trim_sequence_to_codon_length():
    sequence = "ATGCATGCAT"  # Length 10, 1 overhang
    expected_trimmed_sequence = "ATGCATGCA"  # trim 1 base
    assert data_container.trim_sequence_to_codon_length(
        sequence=sequence,
        is_final_cds=True,
        gene_id='gene_1',
        transcript_id='t_1'
    ) == expected_trimmed_sequence
    with pytest.raises(ValueError):
        data_container.trim_sequence_to_codon_length(
            sequence=sequence,
            is_final_cds=False,
            gene_id='gene_1',
            transcript_id='t_1'
        )

    sequence = "ATGCATGCA"  # Length 9, 0 overhang
    expected_trimmed_sequence = "ATGCATGCA"  # no trimming
    assert data_container.trim_sequence_to_codon_length(
        sequence=sequence,
        is_final_cds=True,
        gene_id='gene_1',
        transcript_id='t_1') == expected_trimmed_sequence


def test_construct_peptide_sequences():
    mrna_sequence = "ATGCATGCAT"  # Example mRNA sequence
    cds_coordinates_list = [
        {
            "coordinate": P.open(0, 3),
            "frame": 0,
            "id": "CDS1"
        },
        {
            "coordinate": P.open(3, 6),
            "frame": 0,
            "id": "CDS2"
        }
    ]
    expected_peptide_sequence = "MH"  # Expected translation of ATGCATGCAT
    peptide_sequence, cds_list_tuples = data_container.construct_peptide_sequences(
        gene_id="gene1",
        transcript_id="transcript1",
        mrna_sequence=mrna_sequence,
        cds_coordinates_list=cds_coordinates_list
    )
    assert peptide_sequence == expected_peptide_sequence
