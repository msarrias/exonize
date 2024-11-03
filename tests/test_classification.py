import os
import portion as P
import sqlite3
from exonize.exonize_handler import Exonize
import shutil
from pathlib import Path


def create_exonize_test1():
    results_db_path = Path("mock_results.db")
    if results_db_path.exists():
        os.remove("mock_results.db")
    exonize_obj = Exonize(
        gff_file_path=Path('mock_gff.gff3'),
        genome_file_path=Path('mock_genome.fa'),
        gene_annot_feature='gene',
        cds_annot_feature='CDS',
        transcript_annot_feature='mRNA',
        evalue_threshold=0.01,
        self_hit_threshold=0.5,
        query_coverage_threshold=0.9,
        exon_clustering_overlap_threshold=0.91,
        targets_clustering_overlap_threshold=0.9,
        output_prefix="mock_specie",
        output_directory_path=Path("."),
        )
    shutil.rmtree("mock_specie_exonize", ignore_errors=True)
    exonize_obj.environment.results_database_path = results_db_path
    exonize_obj.data_container.initialize_database()
    gene_hierarchy_dictionary = {
        'gene_0': {
            'coordinate': P.open(0, 2500),
            'chrom': 'X',
            'strand': '+',
            'mRNAs': {
                'transcript_g0_1': {
                    'coordinate': P.open(1, 3000),
                    'strand': '+',
                    'structure': [
                        {'id': 'cds1_g0_t1', 'coordinate': P.open(1, 200), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron1_g0_t1', 'coordinate': P.open(201, 399), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds2_g0_t1', 'coordinate': P.open(400, 500), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron2_g0_t1', 'coordinate': P.open(501, 799), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds3_g0_t1', 'coordinate': P.open(800, 1500), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron3_g0_t1', 'coordinate': P.open(1501, 1799), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds4_g0_t1', 'coordinate': P.open(1800, 2000), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron4_g0_t1', 'coordinate': P.open(2001, 2399), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds5_g0_t1', 'coordinate': P.open(2400, 2500), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron5_g0_t1', 'coordinate': P.open(2501, 2599), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds6_g0_t1', 'coordinate': P.open(2600, 3000), 'frame': 0, 'type': 'CDS'},

                    ]
                },
                'transcript_g0_2': {
                    'coordinate': P.open(500, 1700),
                    'strand': '+',
                    'structure': [
                        {'id': 'cds1_g0_t2', 'coordinate': P.open(600, 700), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron1_g0_t2', 'coordinate': P.open(701, 849), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds2_g0_t2', 'coordinate': P.open(850, 950), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron2_g0_t2', 'coordinate': P.open(951, 999), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds3_g0_t2', 'coordinate': P.open(1000, 1500), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron3_g0_t2', 'coordinate': P.open(1501, 1599), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds4_g0_t2', 'coordinate': P.open(1600, 1700), 'frame': 0, 'type': 'CDS'},

                    ]
                }
            }
        },
        'gene_1': {
            'coordinate': P.open(0, 3000),
            'chrom': 'Y',
            'strand': '+',
            'mRNAs': {
                'transcript_g1_1': {
                    'coordinate': P.open(1, 2400),
                    'strand': '+',
                    'structure': [
                        {'id': 'cds1_g1_t1', 'coordinate': P.open(1, 200), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron1_g1_t1', 'coordinate': P.open(201, 249), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds2_g1_t1', 'coordinate': P.open(250, 350), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron2_g1_t1', 'coordinate': P.open(351, 399), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds3_g1_t1', 'coordinate': P.open(400, 500), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron3_g1_t1', 'coordinate': P.open(501, 549), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds4_g1_t1', 'coordinate': P.open(550, 600), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron4_g1_t1', 'coordinate': P.open(601, 679), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds5_g1_t1', 'coordinate': P.open(680, 780), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron5_g1_t1', 'coordinate': P.open(781, 839), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds6_g1_t1', 'coordinate': P.open(840, 900), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron6_g1_t1', 'coordinate': P.open(901, 949), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds7_g1_t1', 'coordinate': P.open(950, 1000), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron7_g1_t1', 'coordinate': P.open(1001, 1079), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds8_g1_t1', 'coordinate': P.open(1080, 1120), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron8_g1_t1', 'coordinate': P.open(1121, 1199), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds9_g1_t1', 'coordinate': P.open(1200, 1400), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron9_g1_t1', 'coordinate': P.open(1401, 1749), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds10_g1_t1', 'coordinate': P.open(1750, 1900), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron10_g1_t1', 'coordinate': P.open(1901, 2299), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds11_g1_t1', 'coordinate': P.open(2300, 2400), 'frame': 0, 'type': 'CDS'},  #
                    ]
                },
                'transcript_g1_2': {
                    'coordinate': P.open(1, 3000),
                    'strand': '+',
                    'structure': [
                        {'id': 'cds1_g1_t2', 'coordinate': P.open(1, 200), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron0_g1_t2', 'coordinate': P.open(201, 249), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds2_g1_t2', 'coordinate': P.open(250, 350), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron1_g1_t2', 'coordinate': P.open(351, 399), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds3_g1_t2', 'coordinate': P.open(400, 500), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron2_g1_t2', 'coordinate': P.open(501, 549), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds4_g1_t2', 'coordinate': P.open(550, 600), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron3_g1_t2', 'coordinate': P.open(601, 839), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds5_g1_t2', 'coordinate': P.open(840, 900), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron4_g1_t2', 'coordinate': P.open(901, 949), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds6_g1_t2', 'coordinate': P.open(950, 1000), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron5_g1_t2', 'coordinate': P.open(1001, 1199), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds7_g1_t2', 'coordinate': P.open(1200, 1400), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron6_g1_t2', 'coordinate': P.open(1401, 1419), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds8_g1_t2', 'coordinate': P.open(1420, 1460), 'frame': 0, 'type': 'CDS'},  #
                        {'id': 'intron7_g1_t2', 'coordinate': P.open(1461, 1519), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds9_g1_t2', 'coordinate': P.open(1520, 1560), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron8_g1_t2', 'coordinate': P.open(1561, 1749), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds10_g1_t2', 'coordinate': P.open(1750, 2100), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron9_g1_t2', 'coordinate': P.open(2101, 2299), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds11_g1_t2', 'coordinate': P.open(2300, 2400), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron10_g1_t2', 'coordinate': P.open(2401, 2799), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds12_g1_t2', 'coordinate': P.open(2800, 3000), 'frame': 0, 'type': 'CDS'}
                    ]
                },
                'transcript_g1_3': {
                    'coordinate': P.open(1, 3000),
                    'strand': '+',
                    'structure': [
                        {'id': 'cds1_g1_t3', 'coordinate': P.open(1, 200), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron1_g1_t3', 'coordinate': P.open(201, 399), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds2_g1_t3', 'coordinate': P.open(400, 500), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron2_g1_t3', 'coordinate': P.open(501, 1519), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds3_g1_t3', 'coordinate': P.open(1520, 1700), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron3_g1_t3', 'coordinate': P.open(1701, 1799), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds4_g1_t3', 'coordinate': P.open(1800, 2100), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron4_g1_t3', 'coordinate': P.open(2101, 2299), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds5_g1_t3', 'coordinate': P.open(2300, 2400), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron5_g1_t3', 'coordinate': P.open(2401, 2539), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds6_g1_t3', 'coordinate': P.open(2540, 2600), 'frame': 0, 'type': 'CDS'},
                        {'id': 'intron6_g1_t3', 'coordinate': P.open(2601, 2799), 'frame': 0, 'type': 'intron'},
                        {'id': 'cds7_g1_t3', 'coordinate': P.open(2800, 3000), 'frame': 0, 'type': 'CDS'}
                    ]
                }
            }
        }
    }
    mock_gene0 = ['gene_0', 'X', '+', len(gene_hierarchy_dictionary['gene_0']['mRNAs']), 0, 2500]
    fragments_gene0 = [
        ('gene_0', 1, 200, 0, 0, '+', 0, '+', 0, 0, 1e-5, 200, 1, 200, 1300, 1500, '-', '-', '-', 0, 0),
        ('gene_0', 400, 500, 0, 0, '+', 0, '+', 0, 0, 1e-5, 200, 1, 100, 850, 950, '-', '-', '-', 0, 0),
        ('gene_0', 850, 950, 0, 0, '+', 0, '+', 0, 0, 1e-5, 200, 1, 100, 400, 500, '-', '-', '-', 0, 0),
        ('gene_0', 600, 700, 0, 0, '+', 0, '+', 0, 0, 1e-5, 200, 1, 100, 100, 200, '-', '-', '-', 0, 0),
        ('gene_0', 600, 700, 0, 0, '+', 0, '+', 0, 0, 1e-5, 200, 1, 100, 1400, 1500, '-', '-', '-', 0, 0),
        ('gene_0', 600, 700, 0, 0, '+', 0, '+', 0, 0, 1e-5, 200, 1, 100, 1750, 1850, '-', '-', '-', 0, 0),
        ('gene_0', 600, 700, 0, 0, '+', 0, '+', 0, 0, 1e-5, 200, 1, 100, 2200, 2300, '-', '-', '-', 0, 0)
    ]
    mock_gene1 = ['gene_1', 'Y', '+', len(gene_hierarchy_dictionary['gene_1']['mRNAs']), 0, 3000]
    fragments_gene1 = [
        ('gene_1', 1, 200, 0, 0, '+', 0, '+', 0, 0, 1e-5, 200, 1, 200, 1200, 1400, '-', '-', '-', 0, 0),
        ('gene_1', 1200, 1400, 0, 0, '+', 0, '+', 0, 0, 1e-5, 200, 1, 200, 1, 200, '-', '-', '-', 0, 0),
        ('gene_1', 400, 500, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 1, 100, 2300, 2400, '-', '-', '-', 0, 0),
        ('gene_1', 2300, 2400, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 1, 100, 400, 500, '-', '-', '-', 0, 0),
        ('gene_1', 840, 900, 0, 0, '+', 0, '+', 0, 0, 1e-5, 60, 1, 60, 2540, 2600, '-', '-', '-', 0, 0),
        ('gene_1', 2540, 2600, 0, 0, '+', 0, '+', 0, 0, 1e-5, 60, 1, 60, 840, 900, '-', '-', '-', 0, 0),
        ('gene_1', 550, 600, 0, 0, '+', 0, '+', 0, 0, 1e-5, 50, 1, 50, 950, 1000, '-', '-', '-', 0, 0),
        ('gene_1', 950, 1000, 0, 0, '+', 0, '+', 0, 0, 1e-5, 50, 1, 50, 550, 600, '-', '-', '-', 0, 0),
        ('gene_1', 250, 350, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 1, 100, 680, 780, '-', '-', '-', 0, 0),
        ('gene_1', 680, 780, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 1, 100, 250, 350, '-', '-', '-', 0, 0),
        ('gene_1', 1080, 1120, 0, 0, '+', 0, '+', 0, 0, 1e-5, 40, 1, 40, 1420, 1460, '-', '-', '-', 0, 0),
        ('gene_1', 1420, 1460, 0, 0, '+', 0, '+', 0, 0, 1e-5, 40, 1, 40, 1080, 1120, '-', '-', '-', 0, 0),
        ('gene_1', 1750, 1900, 0, 0, '+', 0, '+', 0, 0, 1e-5, 150, 1, 150, 1210, 1360, '-', '-', '-', 0, 0)
    ]
    exonize_obj.data_container.initialize_database()
    exonize_obj.database_interface.insert_matches(
        gene_args_tuple=mock_gene1,
        fragments_tuples_list=fragments_gene1
    )
    exonize_obj.database_interface.insert_matches(
        gene_args_tuple=mock_gene0,
        fragments_tuples_list=fragments_gene0
    )
    exonize_obj.event_classifier.data_container.gene_hierarchy_dictionary = gene_hierarchy_dictionary
    matches_list = exonize_obj.database_interface.query_raw_matches()
    exonize_obj.database_interface.insert_identity_and_dna_algns_columns(
        list_tuples=[(1, 1, '', '', i[0]) for i in matches_list]
    )
    exonize_obj.database_interface.insert_percent_query_column_to_fragments()
    exonize_obj.database_interface.create_filtered_full_length_events_view(
        query_overlap_threshold=exonize_obj.environment.query_coverage_threshold,
        evalue_threshold=exonize_obj.environment.evalue_threshold
            )
    exonize_obj.events_reconciliation()
    exonize_obj.transcript_interdependence_classification()
    return exonize_obj, results_db_path


def test_representative_cdss():
    exonize_obj, _ = create_exonize_test1()
    representative_cds_gene_0 = {
        P.open(1, 200),
        P.open(400, 500),
        P.open(600, 700),
        P.open(800, 1500),
        P.open(850, 950),
        P.open(1000, 1500),
        P.open(1600, 1700),
        P.open(1800, 2000),
        P.open(2400, 2500),
        P.open(2600, 3000),
    }

    representative_cds_gene_1 = {
        P.open(1, 200), P.open(250, 350),
        P.open(400, 500), P.open(550, 600),
        P.open(680, 780), P.open(840, 900),
        P.open(950, 1000), P.open(1080, 1120),
        P.open(1200, 1400), P.open(1420, 1460),
        P.open(1520, 1700), P.open(1520, 1560),
        P.open(1750, 1900), P.open(1750, 2100),
        P.open(1800, 2100), P.open(2300, 2400),
        P.open(2540, 2600), P.open(2800, 3000)
    }
    gene_0_rcs = exonize_obj.search_engine.get_candidate_cds_coordinates('gene_0')['candidates_cds_coordinates']
    gene_1_rcs = exonize_obj.search_engine.get_candidate_cds_coordinates('gene_1')['candidates_cds_coordinates']
    assert set(gene_0_rcs) == representative_cds_gene_0
    assert set(gene_1_rcs) == representative_cds_gene_1


def test_expansion():
    exonize_obj, results_db_path = create_exonize_test1()
    expansions = [
        ('gene_0', exonize_obj.environment.full, 1, 200, 1, 0),
        ('gene_0', exonize_obj.environment.partial_insertion, 1300, 1500, 1, 0),
        ('gene_0', exonize_obj.environment.full, 600, 700, 4, 1),
        ('gene_0', exonize_obj.environment.partial_insertion, 100, 200, 1, 1),
        ('gene_0', exonize_obj.environment.partial_insertion, 1400, 1500, 1, 1),
        ('gene_0', exonize_obj.environment.inter_boundary, 1750, 1850, 1, 1),
        ('gene_0', exonize_obj.environment.intronic, 2200, 2300, 1, 1),
        ('gene_0', exonize_obj.environment.full, 400, 500, 2, 2),
        ('gene_0', exonize_obj.environment.full, 850, 950, 2, 2),
        ('gene_1', exonize_obj.environment.full, 1, 200, 2, 0),
        ('gene_1', exonize_obj.environment.full, 1200, 1400, 2, 0),
        ('gene_1', exonize_obj.environment.full, 400, 500, 2, 1),
        ('gene_1', exonize_obj.environment.full, 2300, 2400, 2, 1),
        ('gene_1', exonize_obj.environment.full, 840, 900, 2, 2),
        ('gene_1', exonize_obj.environment.full, 2540, 2600, 2, 2),
        ('gene_1', exonize_obj.environment.full, 550, 600, 2, 3),
        ('gene_1', exonize_obj.environment.full, 950, 1000, 2, 3),
        ('gene_1', exonize_obj.environment.full, 250, 350, 2, 4),
        ('gene_1', exonize_obj.environment.full, 680, 780, 2, 4),
        ('gene_1', exonize_obj.environment.full, 1080, 1120, 2, 5),
        ('gene_1', exonize_obj.environment.full, 1420, 1460, 2, 5),
        ('gene_1', exonize_obj.environment.full, 1750, 1900, 1, 6),
        ('gene_1', exonize_obj.environment.partial_insertion, 1210, 1360, 1, 6)
    ]
    with sqlite3.connect(results_db_path) as db:
        cursor = db.cursor()
        cursor.execute(
            """
            SELECT GeneID, Mode, EventStart, EventEnd, EventDegree, ExpansionID FROM Expansions
            """
        )
        records = cursor.fetchall()
        events_gene_0 = len(set([record[-1] for record in records if record[0] == 'gene_0']))
        expected_events_gene_0 = len(set([i[-1] for i in [rec for rec in expansions if rec[0] == 'gene_0']]))
        events_gene_1 = len(set([record[-1] for record in records if record[0] == 'gene_1']))
        expected_events_gene_1 = len(set([i[-1] for i in [rec for rec in expansions if rec[0] == 'gene_1']]))
        assert set([i[:-1] for i in records]) == set([i[:-1] for i in expansions])
        assert events_gene_0 == expected_events_gene_0
        assert events_gene_1 == expected_events_gene_1


def test_matches_interdependence_counts():
    def sort_coordinates(a, b, c, d):
        query, target = sorted(
            [(a, b), (c, d)],
            key=lambda x: (x[0], x[1])
        )
        return query[0], query[1], target[0], target[1]

    _, results_db_path = create_exonize_test1()
    non_reciprocal_matches_count = [
        # gene_1
        (1, 200, 1200, 1400, 'FLEXIBLE'),
        (250, 350, 680, 780, 'OPTIONAL_FLEXIBLE'),
        (400, 500, 2300, 2400, 'OBLIGATE'),
        (550, 600, 950, 1000, 'OPTIONAL_OBLIGATE'),
        (840, 900, 2540, 2600, 'EXCLUSIVE'),
        (1080, 1120, 1420, 1460, 'OPTIONAL_EXCLUSIVE'),
        # gene_0
        (400, 500, 850, 950, 'FLEXIBLE')
    ]
    with sqlite3.connect(results_db_path) as db:
        cursor = db.cursor()
        cursor.execute(
            f"""
            SELECT
             QueryExonStart,
             QueryExonEnd,
             TargetStart,
             TargetEnd,
             Classification
            FROM Local_matches_non_reciprocal
            WHERE Mode="FULL";
            """
        )
        records = {
            (sort_coordinates(query_s, query_e, target_s, target_e), class_)
            for (query_s, query_e, target_s, target_e, class_) in cursor.fetchall()
        }
        expected_records = {
            (sort_coordinates(query_s, query_e, target_s, target_e), class_)
            for query_s, query_e, target_s, target_e, class_ in non_reciprocal_matches_count
        }
        assert records == expected_records


# there's missing a test for optional, where all > 0 and intersection is empty
def create_exonize_test2():
    results_db_path2 = Path("mock_results2.db")
    if results_db_path2.exists():
        os.remove("mock_results2.db")
    exonize_obj2 = Exonize(
        gff_file_path=Path('mock_gff.gff3'),
        genome_file_path=Path('mock_genome.fa'),
        gene_annot_feature='gene',
        cds_annot_feature='CDS',
        transcript_annot_feature='mRNA',
        evalue_threshold=0.01,
        self_hit_threshold=0.5,
        query_coverage_threshold=0.9,
        exon_clustering_overlap_threshold=0.91,
        targets_clustering_overlap_threshold=0.9,
        output_directory_path=Path("."),
        output_prefix="mock_specie2",
    )
    shutil.rmtree("mock_specie2_exonize", ignore_errors=True)
    exonize_obj2.environment.results_database_path = results_db_path2
    exonize_obj2.data_container.initialize_database()
    gene_hierarchy_dictionary_expansions_test = {
        'gene1': {
            'coordinate': P.open(0, 1500),
            'mRNAs': {
                'tran1': {
                    'structure': [
                        {'coordinate': P.open(0, 100), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(150, 250), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(300, 500), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(600, 700), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(900, 1000), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(1100, 1300), 'type': 'CDS', 'frame': 0},
                    ]
                },
                'tran2': {
                    'structure': [
                        {'coordinate': P.open(0, 100), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(150, 250), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(600, 700), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(900, 1000), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(1100, 1300), 'type': 'CDS', 'frame': 0},
                    ]
                },
                'tran3': {
                    'structure': [
                        {'coordinate': P.open(0, 100), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(600, 700), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(900, 1000), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(1100, 1300), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(1400, 1500), 'type': 'CDS', 'frame': 0},
                    ]
                }
            }
        },
        'gene2': {
            'coordinate': P.open(0, 4600),
            'mRNAs': {
                'tran1': {
                    'structure': [
                        {'coordinate': P.open(0, 100), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(150, 250), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(300, 500), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(1650, 1750), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(2100, 2200), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(2400, 2500), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(2900, 3000), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(3100, 3200), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(3400, 3500), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(4000, 4100), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(4200, 4300), 'type': 'CDS', 'frame': 0},

                    ]
                },
                'tran2': {
                    'structure': [
                        {'coordinate': P.open(0, 100), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(150, 250), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(2900, 3000), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(3100, 3200), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(3400, 3500), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(4000, 4100), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(4500, 4600), 'type': 'CDS', 'frame': 0},

                    ]
                },
                'tran3': {
                    'structure': [
                        {'coordinate': P.open(600, 700), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(900, 1000), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(1200, 1300), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(1400, 1500), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(1650, 1750), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(2700, 2800), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(4200, 4300), 'type': 'CDS', 'frame': 0},
                        {'coordinate': P.open(4500, 4600), 'type': 'CDS', 'frame': 0},
                    ]
                }
            }
        }}
    mock_gene1 = ['gene1', 'Y', '+',
                  len(gene_hierarchy_dictionary_expansions_test['gene1']['mRNAs']),
                  0, 1500]
    fragments_gene1 = [
        ('gene1', 0, 100, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 300, 400, '-', '-', '-', 0, 0),
        ('gene1', 0, 100, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 700, 800, '-', '-', '-', 0, 0),
        ('gene1', 0, 100, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 150, 250, '-', '-', '-', 0, 0),
        ('gene1', 150, 250, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 0, 100, '-', '-', '-', 0, 0),
        ('gene1', 150, 250, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 300, 400, '-', '-', '-', 0, 0),
        ('gene1', 150, 250, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 700, 800, '-', '-', '-', 0, 0),
        ('gene1', 300, 500, 0, 0, '+', 0, '+', 0, 0, 1e-5, 200, 0, 200, 1000, 1200, '-', '-', '-', 0, 0),
        ('gene1', 600, 700, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 900, 1000, '-', '-', '-', 0, 0),
        ('gene1', 600, 700, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 1100, 1200, '-', '-', '-', 0, 0),
        ('gene1', 900, 1000, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 600, 700, '-', '-', '-', 0, 0),
        ('gene1', 900, 1000, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 1100, 1200, '-', '-', '-', 0, 0),
    ]
    mock_gene2 = ['gene2', 'X', '+',
                  len(gene_hierarchy_dictionary_expansions_test['gene2']['mRNAs']),
                  0, 4600]
    fragments_gene2 = [
        ('gene2', 0, 100, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 150, 250, '-', '-', '-', 0, 0),
        ('gene2', 0, 100, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 600, 700, '-', '-', '-', 0, 0),
        ('gene2', 150, 250, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 0, 100, '-', '-', '-', 0, 0),
        ('gene2', 150, 250, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 600, 700, '-', '-', '-', 0, 0),
        ('gene2', 600, 700, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 150, 250, '-', '-', '-', 0, 0),
        ('gene2', 600, 700, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 0, 100, '-', '-', '-', 0, 0),
        ('gene2', 4000, 4100, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 4200, 4300, '-', '-', '-', 0, 0),
        ('gene2', 4000, 4100, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 4500, 4600, '-', '-', '-', 0, 0),
        ('gene2', 4200, 4300, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 4000, 4100, '-', '-', '-', 0, 0),
        ('gene2', 4200, 4300, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 4500, 4600, '-', '-', '-', 0, 0),
        ('gene2', 4500, 4600, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 4200, 4300, '-', '-', '-', 0, 0),
        ('gene2', 4500, 4600, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 4000, 4100, '-', '-', '-', 0, 0),
        ('gene2', 1200, 1300, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 1400, 1500, '-', '-', '-', 0, 0),
        ('gene2', 1200, 1300, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 1650, 1750, '-', '-', '-', 0, 0),
        ('gene2', 1400, 1500, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 1200, 1300, '-', '-', '-', 0, 0),
        ('gene2', 1400, 1500, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 1650, 1750, '-', '-', '-', 0, 0),
        ('gene2', 1650, 1750, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 1200, 1300, '-', '-', '-', 0, 0),
        ('gene2', 1650, 1750, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 1400, 1500, '-', '-', '-', 0, 0),
        ('gene2', 2100, 2200, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 2400, 2500, '-', '-', '-', 0, 0),
        ('gene2', 2100, 2200, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 2700, 2800, '-', '-', '-', 0, 0),
        ('gene2', 2900, 3000, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 3100, 3200, '-', '-', '-', 0, 0),
        ('gene2', 2900, 3000, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 3400, 3500, '-', '-', '-', 0, 0),
        ('gene2', 3400, 3500, 0, 0, '+', 0, '+', 0, 0, 1e-5, 100, 0, 100, 2900, 3000, '-', '-', '-', 0, 0)
    ]

    test_events = [
        # FLEXIBLE
        ('gene1', exonize_obj2.environment.full, 0, 100, 2, None, 0),
        ('gene1', exonize_obj2.environment.partial_insertion, 300, 400, 1, None, 0),
        ('gene1', exonize_obj2.environment.intronic, 700, 800, 1, None, 0),

        # FLEXIBLE
        ('gene1', exonize_obj2.environment.full, 0, 100, 3, None, 1),
        ('gene1', exonize_obj2.environment.full, 150, 250, 3, None, 1),
        ('gene1', exonize_obj2.environment.partial_insertion, 300, 400, 2, None, 1),
        ('gene1', exonize_obj2.environment.intronic, 700, 800, 2, None, 1),

        # OBLIGATE
        ('gene1', exonize_obj2.environment.full, 600, 700, 2, None, 2),
        ('gene1', exonize_obj2.environment.full, 900, 1000, 2, None, 2),
        ('gene1', exonize_obj2.environment.partial_insertion, 1100, 1200, 2, None, 2),

        # EXCLUSIVE
        ('gene2', exonize_obj2.environment.full, 0, 100, 2, None, 0),
        ('gene2', exonize_obj2.environment.full, 150, 250, 2, None, 0),
        ('gene2', exonize_obj2.environment.full, 600, 700, 2, None, 0),

        # EXCLUSIVE
        ('gene2', exonize_obj2.environment.full, 4000, 4100, 2, None, 1),
        ('gene2', exonize_obj2.environment.full, 4200, 4300, 2, None, 1),
        ('gene2', exonize_obj2.environment.full, 4500, 4600, 2, None, 1),

        # OPTIONAL - FLEXIBLE
        ('gene2', exonize_obj2.environment.full, 1200, 1300, 2, None, 2),
        ('gene2', exonize_obj2.environment.full, 1400, 1500, 2, None, 2),
        ('gene2', exonize_obj2.environment.full, 1650, 1750, 2, None, 2),

        # OPTIONAL - EXCLUSIVE
        ('gene2', exonize_obj2.environment.full, 2100, 2200, 2, None, 3),
        ('gene2', exonize_obj2.environment.full, 2400, 2500, 2, None, 3),
        ('gene2', exonize_obj2.environment.full, 2700, 2800, 2, None, 3),

        # OPTIONAL - OBLIGATE
        ('gene2', exonize_obj2.environment.full, 2900, 3000, 2, None, 4),
        ('gene2', exonize_obj2.environment.full, 3100, 3200, 2, None, 4),
        ('gene2', exonize_obj2.environment.full, 3400, 3500, 2, None, 4),

        # NON-CODING
        ('gene1', exonize_obj2.environment.full, 300, 500, 1, None, 3),
        ('gene1', exonize_obj2.environment.intronic, 1000, 1200, 1, None, 3),
    ]
    exonize_obj2.database_interface.insert_matches(
        gene_args_tuple=mock_gene1,
        fragments_tuples_list=fragments_gene1
    )
    exonize_obj2.database_interface.insert_matches(
        gene_args_tuple=mock_gene2,
        fragments_tuples_list=fragments_gene2
    )
    exonize_obj2.data_container.gene_hierarchy_dictionary = gene_hierarchy_dictionary_expansions_test
    matches_list2 = exonize_obj2.database_interface.query_raw_matches()
    exonize_obj2.database_interface.insert_identity_and_dna_algns_columns(
        list_tuples=[(1, 1, '', '', i[0]) for i in matches_list2]
    )
    exonize_obj2.database_interface.insert_percent_query_column_to_fragments()

    exonize_obj2.database_interface.create_filtered_full_length_events_view(
        query_overlap_threshold=exonize_obj2.environment.query_coverage_threshold,
        evalue_threshold=exonize_obj2.environment.evalue_threshold
    )
    exonize_obj2.database_interface.create_non_reciprocal_fragments_table()
    tblastx_full_matches_list = exonize_obj2.database_interface.query_full_length_events()
    exonize_obj2.full_matches_dictionary = exonize_obj2.event_reconciler.get_gene_events_dictionary(
        local_full_matches_list=tblastx_full_matches_list
    )
    genes_list = list(exonize_obj2.full_matches_dictionary.keys())
    for gene_id in genes_list:
        tblastx_records_set = exonize_obj2.full_matches_dictionary[gene_id]
        cds_candidates_dictionary = exonize_obj2.search_engine.get_candidate_cds_coordinates(
            gene_id=gene_id
        )
        (query_coordinates,
         targets_reference_coordinates_dictionary
         ) = exonize_obj2.event_reconciler.align_target_coordinates(
            gene_id=gene_id,
            local_records_set=tblastx_records_set,
            cds_candidates_dictionary=cds_candidates_dictionary
        )
        gene_graph = exonize_obj2.event_reconciler.create_events_multigraph(
            local_records_set=tblastx_records_set,
            global_records_set=set(),
            query_local_coordinates_set=query_coordinates,
            targets_reference_coordinates_dictionary=targets_reference_coordinates_dictionary
        )
        (gene_events_list,
         non_reciprocal_fragment_ids_list,
         full_events_list
         ) = exonize_obj2.event_reconciler.get_reconciled_graph_and_expansion_events_tuples(
            targets_reference_coordinates_dictionary=targets_reference_coordinates_dictionary,
            gene_id=gene_id,
            gene_graph=gene_graph
        )
        if full_events_list:
            expansions_dictionary = exonize_obj2.event_reconciler.build_expansion_dictionary(
                records=full_events_list
            )
            tuples_to_insert = exonize_obj2.event_reconciler.get_gene_full_events_tandemness_tuples(
                expansions_dictionary
            )
            exonize_obj2.database_interface.insert_expansion_table(
                list_tuples=gene_events_list,
                list_tuples_full=full_events_list,
                list_tuples_tandemness=tuples_to_insert,

            )
        exonize_obj2.database_interface.insert_in_non_reciprocal_fragments_table(
            fragment_ids_list=non_reciprocal_fragment_ids_list,
            gene_id=gene_id
        )

    exonize_obj2.transcript_interdependence_classification()
    return exonize_obj2, results_db_path2


def test_expansion_transcript_iterdependence_classification():
    _, results_db_path2 = create_exonize_test2()
    expected_expansions_classification = [
        ('gene1', 1, 3, 2, 4, 1, 1, 0, 'FLEXIBLE', ''),  # n x (k + 1) = 3 x (1 + 1) = 6
        ('gene1', 2, 3, 2, 6, 0, 0, 0, 'OBLIGATE', ''),
        ('gene2', 0, 3, 3, 0, 5, 4, 0, 'EXCLUSIVE',
         '_'.join([
             str(i)
             for i in (
                 P.open(600, 700),
                 tuple((P.open(0, 100), P.open(150, 250)))
             )
         ])),
        ('gene2', 1, 3, 3, 0, 6, 3, 0, 'EXCLUSIVE',
         '_'.join([str(i) for i in
                   (P.open(4500, 4600),
                    P.open(4200, 4300),
                    P.open(4000, 4100))
                   ])),
        ('gene2', 2, 3, 3, 3, 1, 2, 3, 'OPTIONAL_FLEXIBLE', ''),
        ('gene2', 3, 3, 3, 0, 3, 3, 3, 'OPTIONAL_EXCLUSIVE',
         '_'.join([str(i) for i
                   in (P.open(2700, 2800),
                       (P.open(2100, 2200),
                        P.open(2400, 2500)))
                   ])),
        ('gene2', 4, 3, 3, 6, 0, 0, 3, 'OPTIONAL_OBLIGATE', ''),
    ]
    with sqlite3.connect(results_db_path2) as db:
        cursor = db.cursor()
        cursor.execute(
            """
            SELECT
                GeneID,
                NumberTranscripts,
                NumberCodingEvents,
                All_,
                Present,
                Absent,
                Neither,
                Classification,
                ExclusiveEvents
            FROM Expansions_transcript_interdependence;
            """
        )
        records = cursor.fetchall()
        assert set(records) == set([(i[0], *i[2:]) for i in expected_expansions_classification])


def test_expansion_full_events():
    exonize_obj, results_db_path2 = create_exonize_test2()
    exp_full_events = [
        ('gene2', exonize_obj.environment.full, 0, 100, 2),
        ('gene2', exonize_obj.environment.full, 150, 250, 2),
        ('gene2', exonize_obj.environment.full, 600, 700, 2),

        ('gene2', exonize_obj.environment.full, 4000, 4100, 2),
        ('gene2', exonize_obj.environment.full, 4200, 4300, 2),
        ('gene2', exonize_obj.environment.full, 4500, 4600, 2),

        ('gene2', exonize_obj.environment.full, 1200, 1300, 2),
        ('gene2', exonize_obj.environment.full, 1400, 1500, 2),
        ('gene2', exonize_obj.environment.full, 1650, 1750, 2),

        ('gene2', exonize_obj.environment.full, 2100, 2200, 2),
        ('gene2', exonize_obj.environment.full, 2400, 2500, 1),
        ('gene2', exonize_obj.environment.full, 2700, 2800, 1),

        ('gene2', exonize_obj.environment.full, 2900, 3000, 2),
        ('gene2', exonize_obj.environment.full, 3100, 3200, 1),
        ('gene2', exonize_obj.environment.full, 3400, 3500, 1),

        ('gene1', exonize_obj.environment.full, 0, 100, 1),
        ('gene1', exonize_obj.environment.full, 150, 250, 1),

        ('gene1', exonize_obj.environment.full, 600, 700, 1),
        ('gene1', exonize_obj.environment.full, 900, 1000, 1),
    ]
    with sqlite3.connect(results_db_path2) as db:
        cursor = db.cursor()
        cursor.execute(
            """
            select GeneID, Mode, EventStart, EventEnd, EventDegree
            from Expansions_full;
            """
        )
        records = cursor.fetchall()
        assert set(records) == set(exp_full_events)


def test_expansion_full_events_tandemness():
    _, results_db_path2 = create_exonize_test2()
    exp_full_events_tandemness = [
        ('gene2', 0, 100, 150, 250, 1),
        ('gene2', 150, 250, 600, 700, 0),
        ('gene2', 4000, 4100, 4200, 4300, 1),
        ('gene2', 4200, 4300, 4500, 4600, 1),
        ('gene2', 1200, 1300, 1400, 1500, 1),
        ('gene2', 1400, 1500, 1650, 1750, 1),
        ('gene2', 2100, 2200, 2400, 2500, 1),
        ('gene2', 2400, 2500, 2700, 2800, 1),
        ('gene2', 2900, 3000, 3100, 3200, 1),
        ('gene2', 3100, 3200, 3400, 3500, 1),
        ('gene1', 0, 100, 150, 250, 1),
        ('gene1', 600, 700, 900, 1000, 1)
    ]
    with sqlite3.connect(results_db_path2) as db:
        cursor = db.cursor()
        cursor.execute(
            """
            SELECT 
                GeneID,
                PredecessorStart,
                PredecessorEnd,
                SuccessorStart,
                SuccessorEnd,
                TandemPair
            FROM Expansions_full_tandem;
            """
        )
        records = cursor.fetchall()
        assert set(records) == set(exp_full_events_tandemness)
