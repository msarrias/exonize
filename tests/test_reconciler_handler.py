from unittest.mock import Mock
from exonize.reconciler_handler import ReconcilerHandler
from exonize.data_preprocessor import DataPreprocessor
from exonize.blast_searcher import BLASTsearcher
from pathlib import Path
import portion as P

data_container = DataPreprocessor(
    gene_annot_feature='gene',
    cds_annot_feature='CDS',
    transcript_annot_feature='mRNA',
    min_exon_length=20,
    logger_obj=Mock(),
    database_interface=Mock(),
    working_directory=Path(''),
    gff_file_path=Path(''),
    output_prefix='test',
    genome_file_path=Path(''),
    debug_mode=False,
    csv=False
)

blast_engine = BLASTsearcher(
    data_container=data_container,
    sleep_max_seconds=40,
    self_hit_threshold=0.5,
    min_exon_length=20,
    exon_clustering_overlap_threshold=0.8,
    debug_mode=False,
)

counter_handler = ReconcilerHandler(
    blast_engine=blast_engine,
    targets_clustering_overlap_threshold=0.8,
    query_coverage_threshold=0.8,
)


def test_build_reference_dictionary():
    query_coordinates = {
        P.open(0, 100),
        P.open(200, 250)
    }
    target_coordinates = {
        (P.open(0, 50), 0.9),
        (P.open(0, 48), 0.93),
        (P.open(40, 90), 0.8),
        (P.open(200, 250), 0.7),
        (P.open(210, 250), 0.7),
        (P.open(220, 250), 0.7),
        (P.open(220, 270), 0.6),
        (P.open(215, 270), 0.7),
        (P.open(219, 270), 0.8),
        (P.open(400, 450), 0.4),
        (P.open(402, 450), 0.5),
        (P.open(420, 450), 0.5)
    }
    cds_candidates_dictionary = {
        'candidates_cds_coordinates': query_coordinates
    }
    overlapping_targets = counter_handler.data_container.get_overlapping_clusters(
        target_coordinates_set=target_coordinates,
        threshold=counter_handler.targets_clustering_overlap_threshold
    )

    expected_output = {
        P.open(200, 250): {
            'reference': P.open(200, 250),
            'mode': 'FULL'
        },
        P.open(210, 250): {
            'reference': P.open(200, 250),
            'mode': 'FULL'
        },
        P.open(220, 250): {
            'reference': P.open(220, 250),
            'mode': 'INSERTION'
        },
        P.open(0, 50): {
            'reference': P.open(0, 50),
            'mode': 'INSERTION'
        },
        P.open(0, 48): {
            'reference': P.open(0, 50),
            'mode': 'INSERTION'
        },
        P.open(40, 90): {
            'reference': P.open(40, 90),
            'mode': 'INSERTION'
        },
        P.open(220, 270): {
            'reference': P.open(220, 270),
            'mode': 'TRUNCATION'
        },
        P.open(215, 270): {
            'reference': P.open(220, 270),
            'mode': 'TRUNCATION'
        },
        P.open(219, 270): {
            'reference': P.open(220, 270),
            'mode': 'TRUNCATION'
        },
        P.open(400, 450): {
            'reference': P.open(400, 450),
            'mode': 'CANDIDATE'
        },
        P.open(402, 450): {
            'reference': P.open(400, 450),
            'mode': 'CANDIDATE'
        },
        P.open(420, 450): {
            'reference': P.open(420, 450),
            'mode': 'CANDIDATE'
        }
        # Assuming no CDS overlap
    }
    assert counter_handler.get_matches_reference_mode_dictionary(
        cds_candidates_set=set(cds_candidates_dictionary['candidates_cds_coordinates']),
        clusters_list=overlapping_targets,
        gene_cds_set=set()
    ) == expected_output
