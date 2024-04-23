import portion as P
from exonize.exonize_handler import Exonize

gene_hierarchy_dictionary = {
    'gene_1': {
        'coordinate': P.open(1, 3000),
        'chrom': 'Y',
        'strand': '+',
        'mRNAs': {
            'transcript_g1_1': {
                'coordinate': P.open(1, 3000),
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
                    {'id': 'intron9_g1_t1', 'coordinate': P.open(1401, 1799), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds10_g1_t1', 'coordinate': P.open(1800, 2100), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron10_g1_t1', 'coordinate': P.open(2101, 2299), 'frame': 0, 'type': 'intron'}
                ]
            },
            'transcript_g1_2': {
                'coordinate': P.open(1, 3000),
                'strand': '+',
                'structure': [
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
                    {'id': 'cds8_g1_t2', 'coordinate': P.open(1420, 1450), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron7_g1_t2', 'coordinate': P.open(1451, 1519), 'frame': 0, 'type': 'intron'},
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
                    {'id': 'transcript_g1_3', 'coordinate': P.open(1, 3000), 'frame': 0, 'type': 'mRNA'},
                    {'id': 'cds1_g1_t3', 'coordinate': P.open(400, 500), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron1_g1_t3', 'coordinate': P.open(501, 1519), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds2_g1_t3', 'coordinate': P.open(1520, 1700), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron2_g1_t3', 'coordinate': P.open(1701, 1799), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds3_g1_t3', 'coordinate': P.open(1800, 2100), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron3_g1_t3', 'coordinate': P.open(2101, 2299), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds4_g1_t3', 'coordinate': P.open(2300, 2400), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron4_g1_t3', 'coordinate': P.open(2401, 2539), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds5_g1_t3', 'coordinate': P.open(2540, 2600), 'frame': 0, 'type': 'CDS'},
                    {'id': 'intron5_g1_t3', 'coordinate': P.open(2601, 2799), 'frame': 0, 'type': 'intron'},
                    {'id': 'cds6_g1_t3', 'coordinate': P.open(2800, 3000), 'frame': 0, 'type': 'CDS'}
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
    ('gene_1', 1080, 1120, 0, 0, '+', 0, '+', 0, 0, 1e-5, 40, 1, 40, 1420, 1450, '-', '-', '-', 0, 0),
    ('gene_1', 1420, 1450, 0, 0, '+', 0, '+', 0, 0, 1e-5, 40, 1, 40, 1080, 1120, '-', '-', '-', 0, 0)
]


def main():
    exonize_obj = Exonize(
        gff_file_path='mock_gff.gff3',
        genome_file_path='mock_genome.fa',
        specie_identifier="human",
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
    exonize_obj.database_interface.results_database_path = "mock_results_db.db"
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


if __name__ == '__main__':
    main()
