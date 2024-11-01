from unittest.mock import Mock
from exonize.exonize import Exonize
from pathlib import Path
import portion as P


exonize_obj = Exonize(
    gff_file_path=Path('mock_gff.gff3'),
    genome_file_path=Path('mock_genome.fa'),
    gene_annot_feature='gene',
    cds_annot_feature='CDS',
    transcript_annot_feature='mRNA',
    sequence_base=1,
    frame_base=0,
    evalue_threshold=0.01,
    self_hit_threshold=0.5,
    query_coverage_threshold=0.8,
    min_exon_length=20,
    exon_clustering_overlap_threshold=0.8,
    targets_clustering_overlap_threshold=0.9,
    output_prefix="mock_specie",
    csv=False,
    enable_debug=False,
    soft_force=False,
    hard_force=False,
    sleep_max_seconds=0,
    cpus_number=1,
    timeout_database=60,
    output_directory_path=Path("."),
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
    overlapping_targets = exonize_obj.data_container.get_overlapping_clusters(
        target_coordinates_set=target_coordinates,
        threshold=exonize_obj.environment.targets_clustering_overlap_threshold
    )

    expected_output = {
        P.open(200, 250): {
            'reference': P.open(200, 250),
            'mode': exonize_obj.environment.full
        },
        P.open(210, 250): {
            'reference': P.open(200, 250),
            'mode': exonize_obj.environment.full
        },
        P.open(220, 250): {
            'reference': P.open(220, 250),
            'mode': exonize_obj.environment.partial_insertion
        },
        P.open(0, 50): {
            'reference': P.open(0, 50),
            'mode': exonize_obj.environment.partial_insertion
        },
        P.open(0, 48): {
            'reference': P.open(0, 50),
            'mode': exonize_obj.environment.partial_insertion
        },
        P.open(40, 90): {
            'reference': P.open(40, 90),
            'mode': exonize_obj.environment.partial_insertion
        },
        P.open(220, 270): {
            'reference': P.open(220, 270),
            'mode': exonize_obj.environment.inter_boundary
        },
        P.open(215, 270): {
            'reference': P.open(220, 270),
            'mode': exonize_obj.environment.inter_boundary
        },
        P.open(219, 270): {
            'reference': P.open(220, 270),
            'mode': exonize_obj.environment.inter_boundary
        },
        P.open(400, 450): {
            'reference': P.open(400, 450),
            'mode': exonize_obj.environment.intronic
        },
        P.open(402, 450): {
            'reference': P.open(400, 450),
            'mode': exonize_obj.environment.intronic
        },
        P.open(420, 450): {
            'reference': P.open(420, 450),
            'mode': exonize_obj.environment.intronic
        }
        # Assuming no CDS overlap
    }
    assert exonize_obj.event_reconciler.get_matches_reference_mode_dictionary(
        cds_candidates_set=set(cds_candidates_dictionary['candidates_cds_coordinates']),
        clusters_list=overlapping_targets,
        gene_cds_set=set()
    ) == expected_output
