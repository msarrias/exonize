from exonize.data_preprocessor import DataPreprocessor
from unittest.mock import Mock
import portion as P
from pathlib import Path
from Bio.Seq import Seq
import pytest

data_container = DataPreprocessor(
    gene_annot_feature='gene',
    cds_annot_feature='CDS',
    transcript_annot_feature='mRNA',
    sequence_base=1,
    frame_base=0,
    min_exon_length=20,
    logger_obj=Mock(),
    database_interface=Mock(),
    working_directory=Path(''),
    gff_file_path=Path(''),
    output_prefix='test',
    genome_file_path=Path(''),
    debug_mode=False,
    global_search=False,
    local_search=False,
    csv=False
)

data_container.genome_dictionary = {
        "chr1": "ATGC" * 100  # Simulate a genome sequence for testing
    }


def test_resolve_overlaps_between_coordinates():
    cds_overlapping_threshold = 0.8
    test = [
        P.open(0, 100),
        P.open(180, 300),
        P.open(200, 300),
        P.open(600, 900),
        P.open(700, 2000)
    ]
    res_a = [
        P.open(0, 100),
        P.open(200, 300),
        P.open(600, 900),
        P.open(700, 2000)
    ]
    clusters = data_container.get_overlapping_clusters(
        target_coordinates_set=set((coordinate, None) for coordinate in test),
        threshold=cds_overlapping_threshold
    )

    assert set(data_container.flatten_clusters_representative_exons(
            cluster_list=clusters)) == set(res_a)

    cds_overlapping_threshold = 0.3
    res_b = [
        P.open(0, 100),
        P.open(200, 300),
        P.open(600, 900),
        P.open(700, 2000)
    ]
    clusters = data_container.get_overlapping_clusters(
        target_coordinates_set=set((coordinate, None) for coordinate in test),
        threshold=cds_overlapping_threshold
    )

    assert set(data_container.flatten_clusters_representative_exons(
            cluster_list=clusters)) == set(res_b)

    cds_overlapping_threshold = 0.001
    res_c = [
        P.open(0, 100),
        P.open(200, 300),
        P.open(600, 900)
    ]
    clusters = data_container.get_overlapping_clusters(
        target_coordinates_set=set((coordinate, None) for coordinate in test),
        threshold=cds_overlapping_threshold
    )

    assert set(data_container.flatten_clusters_representative_exons(
            cluster_list=clusters)) == set(res_c)


def test_get_overlapping_clusters():
    # Case 1: Overlapping intervals
    target_coordinates_set_1 = {
        (P.open(0, 50), 0.9),
        (P.open(40, 100), 0.8),
        (P.open(200, 300), 0.7)
    }
    expected_clusters_1 = [
        [(P.open(0, 50), 0.9),
         (P.open(40, 100), 0.8)],
        [(P.open(200, 300), 0.7)]
    ]

    # Case 2: Non-overlapping intervals
    target_coordinates_set_2 = {
        (P.open(0, 50), 0.9),
        (P.open(60, 110), 0.8),
        (P.open(120, 170), 0.7)
    }
    expected_clusters_2 = [
        [(P.open(0, 50), 0.9)],
        [(P.open(60, 110), 0.8)],
        [(P.open(120, 170), 0.7)]
    ]

    # Test and assert
    assert data_container.get_overlapping_clusters(
        target_coordinates_set=target_coordinates_set_1,
        threshold=0
    ) == expected_clusters_1
    assert data_container.get_overlapping_clusters(
        target_coordinates_set=target_coordinates_set_2,
        threshold=0
    ) == expected_clusters_2


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
