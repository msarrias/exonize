from unittest.mock import Mock
from exonize.counter_handler import CounterHandler
from exonize.blast_searcher import BLASTsearcher
import portion as P

blast_engine = BLASTsearcher(
    data_container=Mock(),
    sleep_max_seconds=40,
    self_hit_threshold=0.5,
    min_exon_length=20,
    cds_overlapping_threshold=0.8,
    evalue_threshold=1e-5,
    debug_mode=False,

)
counter_handler = CounterHandler(
    blast_engine=blast_engine,
    database_interface=Mock(),
    cds_overlapping_threshold=0.8,
    draw_event_multigraphs=False,
)


def test_get_overlapping_clusters():
    # Case 1: Overlapping intervals
    target_coordinates_set_1 = {
        (P.open(0, 50), 0.9),
        (P.open(40, 100), 0.8),
        (P.open(200, 300), 0.7)
    }
    expected_clusters_1 = [
        [(P.open(0, 50), 0.9), (P.open(40, 100), 0.8)],
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
    assert counter_handler.get_overlapping_clusters(
        target_coordinates_set=target_coordinates_set_1,
        threshold=0
    ) == expected_clusters_1
    assert counter_handler.get_overlapping_clusters(
        target_coordinates_set=target_coordinates_set_2,
        threshold=0.
    ) == expected_clusters_2


def test_build_reference_dictionary():
    query_coordinates = {
        P.open(0, 100),
        P.open(200, 250)
    }
    target_coordinates = {
        (P.open(0, 50), 0.9),
        (P.open(40, 90), 0.8),
        (P.open(200, 250), 0.7),
        (P.open(210, 250), 0.7),
        (P.open(220, 270), 0.6),
        (P.open(400, 450), 0.5)
    }
    cds_candidates_dictionary = {
        'candidates_cds_coordinates': query_coordinates
    }
    overlapping_targets = counter_handler.get_overlapping_clusters(
        target_coordinates_set=target_coordinates,
        threshold=counter_handler.cds_overlapping_threshold
    )

    expected_output = {
        P.open(200, 250): {
            'reference_coordinate': P.open(200, 250),
            'reference_type': 'full'
        },
        P.open(210, 250): {
            'reference_coordinate': P.open(210, 250),
            'reference_type': 'insertion'
        },
        P.open(0, 50): {
            'reference_coordinate': P.open(0, 50),
            'reference_type': 'insertion'
        },
        P.open(40, 90): {
            'reference_coordinate': P.open(40, 90),
            'reference_type': 'insertion'
        },
        P.open(220, 270): {
            'reference_coordinate': P.open(220, 270),
            'reference_type': 'truncation'
        },
        P.open(400, 450): {
            'reference_coordinate': P.open(400, 450),
            'reference_type': 'deactivated'
        }  # Assuming no CDS overlap
    }

    assert counter_handler.build_reference_dictionary(
        cds_candidates_dictionary=cds_candidates_dictionary,
        clusters_list=overlapping_targets
    ) == expected_output


def test_assign_event_ids():
    pass
