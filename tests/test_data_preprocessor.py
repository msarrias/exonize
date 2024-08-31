from exonize.data_preprocessor import DataPreprocessor
from unittest.mock import Mock
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
    draw_event_multigraphs=False,
    csv=False,
)


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
