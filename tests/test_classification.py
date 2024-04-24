import os.path
import portion as P
import sqlite3
from exonize.exonize_handler import Exonize
import shutil

gene_hierarchy_dictionary = {
    'gene_1': {
        'coordinate': P.open(1, 3000),
        'chrom': 'Y',
        'strand': '+',
        'mRNAs': {
            'transcript_g1_1': {
                'coordinate': P.open(1, 2400),
                'strand': '+',
                'structure': [
                    {'id': 'cds1_g1_t1', 'coordinate': P.open(1, 200), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron1_g1_t1', 'coordinate': P.open(201, 249), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds2_g1_t1', 'coordinate': P.open(250, 350), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron2_g1_t1', 'coordinate': P.open(351, 399), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds3_g1_t1', 'coordinate': P.open(400, 500), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron3_g1_t1', 'coordinate': P.open(501, 549), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds4_g1_t1', 'coordinate': P.open(550, 600), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron4_g1_t1', 'coordinate': P.open(601, 679), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds5_g1_t1', 'coordinate': P.open(680, 780), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron5_g1_t1', 'coordinate': P.open(781, 839), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds6_g1_t1', 'coordinate': P.open(840, 900), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron6_g1_t1', 'coordinate': P.open(901, 949), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds7_g1_t1', 'coordinate': P.open(950, 1000), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron7_g1_t1', 'coordinate': P.open(1001, 1079), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds8_g1_t1', 'coordinate': P.open(1080, 1120), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron8_g1_t1', 'coordinate': P.open(1121, 1199), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds9_g1_t1', 'coordinate': P.open(1200, 1400), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron9_g1_t1', 'coordinate': P.open(1401, 2299), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds10_g1_t1', 'coordinate': P.open(2300, 2400), 'frame': 0, 'type': 'CDS'},
                ]
            },
            'transcript_g1_2': {
                'coordinate': P.open(1, 3000),
                'strand': '+',
                'structure': [
                    {'id': 'cds1_g1_t2', 'coordinate': P.open(1, 200), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron0_g1_t2', 'coordinate': P.open(201, 249), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds2_g1_t2', 'coordinate': P.open(250, 350), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron1_g1_t2', 'coordinate': P.open(351, 399), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds3_g1_t2', 'coordinate': P.open(400, 500), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron2_g1_t2', 'coordinate': P.open(501, 549), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds4_g1_t2', 'coordinate': P.open(550, 600), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron3_g1_t2', 'coordinate': P.open(601, 839), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds5_g1_t2', 'coordinate': P.open(840, 900), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron4_g1_t2', 'coordinate': P.open(901, 949), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds6_g1_t2', 'coordinate': P.open(950, 1000), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron5_g1_t2', 'coordinate': P.open(1001, 1199), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds7_g1_t2', 'coordinate': P.open(1200, 1400), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron6_g1_t2', 'coordinate': P.open(1401, 1419), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds8_g1_t2', 'coordinate': P.open(1420, 1460), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron7_g1_t2', 'coordinate': P.open(1461, 1519), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds9_g1_t2', 'coordinate': P.open(1520, 1560), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron8_g1_t2', 'coordinate': P.open(1561, 1799), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds10_g1_t2', 'coordinate': P.open(1800, 2100), 'frame': 0, 'type': 'CDS'},
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
mock_gene = ['gene_1', 'Y', '+', 3, 1, 3000, 1]
fragments = [
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
    ('gene_1', 1420, 1460, 0, 0, '+', 0, '+', 0, 0, 1e-5, 40, 1, 40, 1080, 1120, '-', '-', '-', 0, 0)
]

matches = [
    # FLEXIBLE
    (1, 1, 'gene_1', 'transcript_g1_1', 1, 200, 'cds1_g1_t1', 1, 200, 'FULL', 'cds9_g1_t1',
     1200, 1400, 1200, 1400, 0, 0, 0, 1, 1e-05),
    (2, 1, 'gene_1', 'transcript_g1_2', 1, 200, 'cds1_g1_t2', 1, 200, 'FULL', 'cds7_g1_t2',
     1200, 1400, 1200, 1400, 0, 0, 0, 1, 1e-05),
    (3, 1, 'gene_1', 'transcript_g1_3', 1, 200, 'cds1_g1_t3', 1, 200, 'DEACTIVATED', 'intron2_g1_t3',
     501, 1519, 1200, 1400, 0, 1, 0, 0, 1e-05),
    # FLEXIBLE - 1 RECIPROCAL
    (4, 2, 'gene_1', 'transcript_g1_1', 1200, 1400, 'cds9_g1_t1', 1, 200, 'FULL', 'cds1_g1_t1',
     1, 200, 1, 200, 0, 0, 0, 1, 1e-05),
    (5, 2, 'gene_1', 'transcript_g1_2', 1200, 1400, 'cds7_g1_t2', 1, 200, 'FULL', 'cds1_g1_t2',
     1, 200, 1, 200, 0, 0, 0, 1, 1e-05),
    (6, 2, 'gene_1', 'transcript_g1_3', 1200, 1400, '-', 1, 200, 'FULL', 'cds1_g1_t3',
     1, 200, 1, 200, 0, 0, 1, 0, 1e-05),
    # OBLIGATE
    (7, 3, 'gene_1', 'transcript_g1_1', 400, 500, 'cds3_g1_t1', 1, 100, 'FULL', 'cds10_g1_t1',
     2300, 2400, 2300, 2400, 0, 0, 0, 1, 1e-05),
    (8, 3, 'gene_1', 'transcript_g1_2', 400, 500, 'cds3_g1_t2', 1, 100, 'FULL', 'cds11_g1_t2',
     2300, 2400, 2300, 2400, 0, 0, 0, 1, 1e-05),
    (9, 3, 'gene_1', 'transcript_g1_3', 400, 500, 'cds2_g1_t3', 1, 100, 'FULL', 'cds5_g1_t3',
     2300, 2400, 2300, 2400, 0, 0, 0, 1, 1e-05),
    # OBLIGATE - 3 RECIPROCAL
    (10, 4, 'gene_1', 'transcript_g1_1', 2300, 2400, 'cds10_g1_t1', 1, 100, 'FULL', 'cds3_g1_t1',
     400, 500, 400, 500, 0, 0, 0, 1, 1e-05),
    (11, 4, 'gene_1', 'transcript_g1_2', 2300, 2400, 'cds11_g1_t2', 1, 100, 'FULL', 'cds3_g1_t2',
     400, 500, 400, 500, 0, 0, 0, 1, 1e-05),
    (12, 4, 'gene_1', 'transcript_g1_3', 2300, 2400, 'cds5_g1_t3', 1, 100, 'FULL', 'cds2_g1_t3',
     400, 500, 400, 500, 0, 0, 0, 1, 1e-05),
    # EXCLUSIVE
    (13, 5, 'gene_1', 'transcript_g1_1', 840, 900, 'cds6_g1_t1', 1, 60, 'OUT_OF_MRNA', '-', None, None,
     2540, 2600, 0, 1, 0, 0, 1e-05),
    (14, 5, 'gene_1', 'transcript_g1_2', 840, 900, 'cds5_g1_t2', 1, 60, 'DEACTIVATED', 'intron10_g1_t2',
     2401, 2799, 2540, 2600, 0, 1, 0, 0, 1e-05),
    (15, 5, 'gene_1', 'transcript_g1_3', 840, 900, '-', 1, 60, 'FULL', 'cds6_g1_t3',
     2540, 2600, 2540, 2600, 0, 0, 1, 0, 1e-05),
    # EXCLUSIVE - 5 RECIPROCAL
    (16, 6, 'gene_1', 'transcript_g1_1', 2540, 2600, '-', 1, 60, 'FULL', 'cds6_g1_t1',
     840, 900, 840, 900, 0, 0, 1, 0, 1e-05),
    (17, 6, 'gene_1', 'transcript_g1_2', 2540, 2600, '-', 1, 60, 'FULL', 'cds5_g1_t2',
     840, 900, 840, 900, 0, 0, 1, 0, 1e-05),
    (18, 6, 'gene_1', 'transcript_g1_3', 2540, 2600, 'cds6_g1_t3', 1, 60, 'DEACTIVATED', 'intron2_g1_t3',
     501, 1519, 840, 900, 0, 1, 0, 0, 1e-05),
    # OPTIONAL - OBLIGATE
    (19, 7, 'gene_1', 'transcript_g1_1', 550, 600, 'cds4_g1_t1', 1, 50, 'FULL', 'cds7_g1_t1',
     950, 1000, 950, 1000, 0, 0, 0, 1, 1e-05),
    (20, 7, 'gene_1', 'transcript_g1_2', 550, 600, 'cds4_g1_t2', 1, 50, 'FULL', 'cds6_g1_t2',
     950, 1000, 950, 1000, 0, 0, 0, 1, 1e-05),
    (21, 7, 'gene_1', 'transcript_g1_3', 550, 600, '-', 1, 50, 'NEITHER', 'intron2_g1_t3',
     501, 1519, 950, 1000, 1, 0, 0, 0, 1e-05),
    # OPTIONAL - OBLIGATE - 7 RECIPROCAL
    (22, 8, 'gene_1', 'transcript_g1_1', 950, 1000, 'cds7_g1_t1', 1, 50, 'FULL', 'cds4_g1_t1',
     550, 600, 550, 600, 0, 0, 0, 1, 1e-05),
    (23, 8, 'gene_1', 'transcript_g1_2', 950, 1000, 'cds6_g1_t2', 1, 50, 'FULL', 'cds4_g1_t2',
     550, 600, 550, 600, 0, 0, 0, 1, 1e-05),
    (24, 8, 'gene_1', 'transcript_g1_3', 950, 1000, '-', 1, 50, 'NEITHER', 'intron2_g1_t3',
     501, 1519, 550, 600, 1, 0, 0, 0, 1e-05),
    # OPTIONAL - FLEXIBLE
    (25, 9, 'gene_1', 'transcript_g1_1', 250, 350, 'cds2_g1_t1', 1, 100, 'FULL', 'cds5_g1_t1',
     680, 780, 680, 780, 0, 0, 0, 1, 1e-05),
    (26, 9, 'gene_1', 'transcript_g1_2', 250, 350, 'cds2_g1_t2', 1, 100, 'DEACTIVATED', 'intron3_g1_t2',
     601, 839, 680, 780, 0, 1, 0, 0, 1e-05),
    (27, 9, 'gene_1', 'transcript_g1_3', 250, 350, '-', 1, 100, 'NEITHER', 'intron2_g1_t3',
     501, 1519, 680, 780, 1, 0, 0, 0, 1e-05),
    # OPTIONAL - FLEXIBLE  - 9 RECIPROCAL
    (28, 10, 'gene_1', 'transcript_g1_1', 680, 780, 'cds5_g1_t1', 1, 100, 'FULL', 'cds2_g1_t1',
     250, 350, 250, 350, 0, 0, 0, 1, 1e-05),
    (29, 10, 'gene_1', 'transcript_g1_2', 680, 780, '-', 1, 100, 'FULL', 'cds2_g1_t2',
     250, 350, 250, 350, 0, 0, 1, 0, 1e-05),
    (30, 10, 'gene_1', 'transcript_g1_3', 680, 780, '-', 1, 100, 'NEITHER', 'intron1_g1_t3',
     201, 399, 250, 350, 1, 0, 0, 0, 1e-05),
    # OPTIONAL - EXCLUSIVE
    (31, 11, 'gene_1', 'transcript_g1_1', 1080, 1120, 'cds8_g1_t1', 1, 40, 'DEACTIVATED', 'intron9_g1_t1',
     1401, 2299, 1420, 1460, 0, 1, 0, 0, 1e-05),
    (32, 11, 'gene_1', 'transcript_g1_2', 1080, 1120, '-', 1, 40, 'FULL', 'cds8_g1_t2',
     1420, 1460, 1420, 1460, 0, 0, 1, 0, 1e-05),
    (33, 11, 'gene_1', 'transcript_g1_3', 1080, 1120, '-', 1, 40, 'NEITHER', 'intron2_g1_t3',
     501, 1519, 1420, 1460, 1, 0, 0, 0, 1e-05),
    # OPTIONAL - EXCLUSIVE - 11 RECIPROCAL
    (34, 12, 'gene_1', 'transcript_g1_1', 1420, 1460, '-', 1, 40, 'FULL', 'cds8_g1_t1',
     1080, 1120, 1080, 1120, 0, 0, 1, 0, 1e-05),
    (35, 12, 'gene_1', 'transcript_g1_2', 1420, 1460, 'cds8_g1_t2', 1, 40, 'DEACTIVATED', 'intron5_g1_t2',
     1001, 1199, 1080, 1120, 0, 1, 0, 0, 1e-05),
    (36, 12, 'gene_1', 'transcript_g1_3', 1420, 1460, '-', 1, 40, 'NEITHER', 'intron2_g1_t3',
     501, 1519, 1080, 1120, 1, 0, 0, 0, 1e-05)
]

if os.path.exists("mock_results.db"):
    os.remove("mock_results.db")

exonize_obj = Exonize(
        gff_file_path='mock_gff.gff3',
        genome_file_path='mock_genome.fa',
        specie_identifier="mock_specie",
        draw_event_multigraphs=False,
        enable_debug=False,
        soft_force=False,
        evalue_threshold=0.01,
        sleep_max_seconds=0,
        min_exon_length=30,
        self_hit_threshold=0.5,
        timeout_database=60,
        hard_force=False,
        cds_overlapping_threshold=0.9,
        query_overlapping_threshold=0.9,
        genome_pickled_file_path=".",
        output_directory_path="",
    )
shutil.rmtree("mock_specie_exonize", ignore_errors=True)
exonize_obj.database_interface.results_database_path = "mock_results.db"
exonize_obj.database_interface.connect_create_results_database()
exonize_obj.database_interface.insert_fragments(
    gene_args_tuple=mock_gene,
    fragments_tuples_list=fragments
)
exonize_obj.database_interface.insert_percent_query_column_to_fragments()
exonize_obj.database_interface.create_filtered_full_length_events_view(
    query_overlap_threshold=0.9
)
exonize_obj.event_classifier.data_container.gene_hierarchy_dictionary = gene_hierarchy_dictionary
exonize_obj.event_classifier.identify_full_length_duplications()
exonize_obj.event_classifier.insert_classified_tuples_in_results_database()
exonize_obj.database_interface.create_cumulative_counts_table()
exonize_obj.database_interface.create_exclusive_pairs_view()


def test_matches_transcript_classification():
    with sqlite3.connect("mock_results.db") as db:
        cursor = db.cursor()
        cursor.execute(
            """
            SELECT * FROM Matches_transcript_classification
            ORDER BY gene_id, fragment_id;
            """
        )
        records = cursor.fetchall()
    for my_record, exonize_record in zip(matches, records):
        assert my_record == exonize_record


def test_obligate_evnts():
    with sqlite3.connect("mock_results.db") as db:
        cursor = db.cursor()
        cursor.execute(
            """
            SELECT
            fragment_id,
            gene_id,
            transcript_id,
            query_cds_start,
            query_cds_end,
            query_cds_id,
            query_start,
            query_end,
            type,
            target_cds_id,
            target_cds_start,
            target_cds_end,
            target_start,
            target_end
            FROM Obligate_events
            ORDER BY gene_id, fragment_id
            """
        )
        records = cursor.fetchall()
    obligate_pairs = [record[1:-5] for record in matches if record[-2] == 1]

    for my_record, exonize_record in zip(obligate_pairs, records):
        assert my_record == exonize_record


def test_truncate_results():
    exonize_obj.database_interface.truncate_results_database()
    with sqlite3.connect("mock_results.db") as db:
        cursor = db.cursor()
        cursor.execute(
            """
            SELECT * FROM Truncation_events
            """
        )
        records = cursor.fetchall()
    assert len(records) == 0


def test_exclusive_events():
    with sqlite3.connect("mock_results.db") as db:
        cursor = db.cursor()
        cursor.execute(
            """
            SELECT DISTINCT fragment_id FROM Exclusive_pairs
            """
        )
        records = [fragment_id[0] for fragment_id in cursor.fetchall()]
    assert records == [5, 6]


os.remove("mock_results.db")
