import os
import portion as P
import sqlite3
from exonize.exonize_handler import Exonize
import shutil
from pathlib import Path
from collections import defaultdict


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
expansions = [
    ('gene_0', 'FULL', 1, 200, 1, 0),
    ('gene_0', 'INSERTION_EXCISION', 1300, 1500, 1, 0),
    ('gene_0', 'FULL', 600, 700, 4, 1),
    ('gene_0', 'INSERTION_EXCISION', 100, 200, 1, 1),
    ('gene_0', 'INSERTION_EXCISION', 1400, 1500, 1, 1),
    ('gene_0', 'TRUNCATION_ACQUISITION', 1750, 1850, 1, 1),
    ('gene_0', 'INACTIVE_UNANNOTATED', 2200, 2300, 1, 1),
    ('gene_0', 'FULL', 400, 500, 2, 2),
    ('gene_0', 'FULL', 850, 950, 2, 2),
    ('gene_1', 'FULL', 1, 200, 2, 0),
    ('gene_1', 'FULL', 1200, 1400, 2, 0),
    ('gene_1', 'FULL', 400, 500, 2, 1),
    ('gene_1', 'FULL', 2300, 2400, 2, 1),
    ('gene_1', 'FULL', 840, 900, 2, 2),
    ('gene_1', 'FULL', 2540, 2600, 2, 2),
    ('gene_1', 'FULL', 550, 600, 2, 3),
    ('gene_1', 'FULL', 950, 1000, 2, 3),
    ('gene_1', 'FULL', 250, 350, 2, 4),
    ('gene_1', 'FULL', 680, 780, 2, 4),
    ('gene_1', 'FULL', 1080, 1120, 2, 5),
    ('gene_1', 'FULL', 1420, 1460, 2, 5),
    ('gene_1', 'FULL', 1750, 1900, 1, 6),
    ('gene_1', 'INSERTION_EXCISION', 1210, 1360, 1, 6)
]

matches = [
    # FULL - FLEXIBLE
    (1, 'gene_1', 'transcript_g1_1', 1, 200, 'cds1_g1_t1', 'FULL', 'cds9_g1_t1',
     1200, 1400, 1200, 1400, 0, 0, 0, 1, 1e-05),
    (2, 'gene_1', 'transcript_g1_2', 1, 200, 'cds1_g1_t2', 'FULL', 'cds7_g1_t2',
     1200, 1400, 1200, 1400, 0, 0, 0, 1, 1e-05),
    (3, 'gene_1', 'transcript_g1_3', 1, 200, 'cds1_g1_t3', 'DEACTIVATED', 'intron2_g1_t3',
     501, 1519, 1200, 1400, 0, 1, 0, 0, 1e-05),
    # FULL - FLEXIBLE - 1 RECIPROCAL
    (4, 'gene_1', 'transcript_g1_1', 1200, 1400, 'cds9_g1_t1', 'FULL', 'cds1_g1_t1',
     1, 200, 1, 200, 0, 0, 0, 1, 1e-05),
    (5, 'gene_1', 'transcript_g1_2', 1200, 1400, 'cds7_g1_t2', 'FULL', 'cds1_g1_t2',
     1, 200, 1, 200, 0, 0, 0, 1, 1e-05),
    (6, 'gene_1', 'transcript_g1_3', 1200, 1400, '-', 'FULL', 'cds1_g1_t3',
     1, 200, 1, 200, 0, 0, 1, 0, 1e-05),
    # FULL - OBLIGATE
    (7, 'gene_1', 'transcript_g1_1', 400, 500, 'cds3_g1_t1', 'FULL', 'cds11_g1_t1',
     2300, 2400, 2300, 2400, 0, 0, 0, 1, 1e-05),
    (8,  'gene_1', 'transcript_g1_2', 400, 500, 'cds3_g1_t2', 'FULL', 'cds11_g1_t2',
     2300, 2400, 2300, 2400, 0, 0, 0, 1, 1e-05),
    (9,  'gene_1', 'transcript_g1_3', 400, 500, 'cds2_g1_t3', 'FULL', 'cds5_g1_t3',
     2300, 2400, 2300, 2400, 0, 0, 0, 1, 1e-05),
    # FULL - OBLIGATE - 3 RECIPROCAL
    (10, 'gene_1', 'transcript_g1_1', 2300, 2400, 'cds11_g1_t1', 'FULL', 'cds3_g1_t1',
     400, 500, 400, 500, 0, 0, 0, 1, 1e-05),
    (11, 'gene_1', 'transcript_g1_2', 2300, 2400, 'cds11_g1_t2', 'FULL', 'cds3_g1_t2',
     400, 500, 400, 500, 0, 0, 0, 1, 1e-05),
    (12, 'gene_1', 'transcript_g1_3', 2300, 2400, 'cds5_g1_t3', 'FULL', 'cds2_g1_t3',
     400, 500, 400, 500, 0, 0, 0, 1, 1e-05),
    # FULL - EXCLUSIVE
    (13, 'gene_1', 'transcript_g1_1', 840, 900, 'cds6_g1_t1', 'OUT_OF_MRNA', '-', None, None,
     2540, 2600, 0, 1, 0, 0, 1e-05),
    (14, 'gene_1', 'transcript_g1_2', 840, 900, 'cds5_g1_t2', 'DEACTIVATED', 'intron10_g1_t2',
     2401, 2799, 2540, 2600, 0, 1, 0, 0, 1e-05),
    (15, 'gene_1', 'transcript_g1_3', 840, 900, '-', 'FULL', 'cds6_g1_t3',
     2540, 2600, 2540, 2600, 0, 0, 1, 0, 1e-05),
    # FULL - EXCLUSIVE - 5 RECIPROCAL
    (16, 'gene_1', 'transcript_g1_1', 2540, 2600, '-', 'FULL', 'cds6_g1_t1',
     840, 900, 840, 900, 0, 0, 1, 0, 1e-05),
    (17, 'gene_1', 'transcript_g1_2', 2540, 2600, '-', 'FULL', 'cds5_g1_t2',
     840, 900, 840, 900, 0, 0, 1, 0, 1e-05),
    (18, 'gene_1', 'transcript_g1_3', 2540, 2600, 'cds6_g1_t3', 'DEACTIVATED', 'intron2_g1_t3',
     501, 1519, 840, 900, 0, 1, 0, 0, 1e-05),
    # FULL - OPTIONAL - OBLIGATE
    (19, 'gene_1', 'transcript_g1_1', 550, 600, 'cds4_g1_t1', 'FULL', 'cds7_g1_t1',
     950, 1000, 950, 1000, 0, 0, 0, 1, 1e-05),
    (20, 'gene_1', 'transcript_g1_2', 550, 600, 'cds4_g1_t2', 'FULL', 'cds6_g1_t2',
     950, 1000, 950, 1000, 0, 0, 0, 1, 1e-05),
    (21, 'gene_1', 'transcript_g1_3', 550, 600, '-', 'NEITHER', 'intron2_g1_t3',
     501, 1519, 950, 1000, 1, 0, 0, 0, 1e-05),
    # FULL - OPTIONAL - OBLIGATE - 7 RECIPROCAL
    (22, 'gene_1', 'transcript_g1_1', 950, 1000, 'cds7_g1_t1', 'FULL', 'cds4_g1_t1',
     550, 600, 550, 600, 0, 0, 0, 1, 1e-05),
    (23, 'gene_1', 'transcript_g1_2', 950, 1000, 'cds6_g1_t2', 'FULL', 'cds4_g1_t2',
     550, 600, 550, 600, 0, 0, 0, 1, 1e-05),
    (24, 'gene_1', 'transcript_g1_3', 950, 1000, '-', 'NEITHER', 'intron2_g1_t3',
     501, 1519, 550, 600, 1, 0, 0, 0, 1e-05),
    # FULL - OPTIONAL - FLEXIBLE
    (25, 'gene_1', 'transcript_g1_1', 250, 350, 'cds2_g1_t1', 'FULL', 'cds5_g1_t1',
     680, 780, 680, 780, 0, 0, 0, 1, 1e-05),
    (26, 'gene_1', 'transcript_g1_2', 250, 350, 'cds2_g1_t2', 'DEACTIVATED', 'intron3_g1_t2',
     601, 839, 680, 780, 0, 1, 0, 0, 1e-05),
    (27, 'gene_1', 'transcript_g1_3', 250, 350, '-', 'NEITHER', 'intron2_g1_t3',
     501, 1519, 680, 780, 1, 0, 0, 0, 1e-05),
    # FULL - OPTIONAL - FLEXIBLE  - 9 RECIPROCAL
    (28, 'gene_1', 'transcript_g1_1', 680, 780, 'cds5_g1_t1', 'FULL', 'cds2_g1_t1',
     250, 350, 250, 350, 0, 0, 0, 1, 1e-05),
    (29, 'gene_1', 'transcript_g1_2', 680, 780, '-', 'FULL', 'cds2_g1_t2',
     250, 350, 250, 350, 0, 0, 1, 0, 1e-05),
    (30, 'gene_1', 'transcript_g1_3', 680, 780, '-', 'NEITHER', 'intron1_g1_t3',
     201, 399, 250, 350, 1, 0, 0, 0, 1e-05),
    # FULL - OPTIONAL - EXCLUSIVE
    (31, 'gene_1', 'transcript_g1_1', 1080, 1120, 'cds8_g1_t1', 'DEACTIVATED', 'intron9_g1_t1',
     1401, 1749, 1420, 1460, 0, 1, 0, 0, 1e-05),
    (32, 'gene_1', 'transcript_g1_2', 1080, 1120, '-', 'FULL', 'cds8_g1_t2',
     1420, 1460, 1420, 1460, 0, 0, 1, 0, 1e-05),
    (33, 'gene_1', 'transcript_g1_3', 1080, 1120, '-', 'NEITHER', 'intron2_g1_t3',
     501, 1519, 1420, 1460, 1, 0, 0, 0, 1e-05),
    # FULL - OPTIONAL - EXCLUSIVE - 11 RECIPROCAL
    (34, 'gene_1', 'transcript_g1_1', 1420, 1460, '-', 'FULL', 'cds8_g1_t1',
     1080, 1120, 1080, 1120, 0, 0, 1, 0, 1e-05),
    (35, 'gene_1', 'transcript_g1_2', 1420, 1460, 'cds8_g1_t2', 'DEACTIVATED', 'intron5_g1_t2',
     1001, 1199, 1080, 1120, 0, 1, 0, 0, 1e-05),
    (36, 'gene_1', 'transcript_g1_3', 1420, 1460, '-', 'NEITHER', 'intron2_g1_t3',
     501, 1519, 1080, 1120, 1, 0, 0, 0, 1e-05),
    # INSERTION - OPTIONAL
    (37, 'gene_1', 'transcript_g1_1', 1750, 1900, 'cds10_g1_t1', 'INS_CDS', 'cds9_g1_t1',
     1200, 1400, 1210, 1360, 0, 0, 0, 1, 1e-05),
    (38, 'gene_1', 'transcript_g1_2', 1750, 1900, 'cds10_g1_t2', 'INS_CDS', 'cds7_g1_t2',
     1200, 1400, 1210, 1360, 0, 0, 1, 0, 1e-05),
    (39, 'gene_1', 'transcript_g1_3', 1750, 1900, '-', 'NEITHER', 'intron2_g1_t3',
     501, 1519, 1210, 1360, 1, 0, 0, 0, 1e-05)
]

non_reciprocal_matches_count = [
    # gene_1
    (1, 200, 1200, 1400, 'FLEXIBLE'),
    (250, 350, 680, 780, 'OPTIONAL_FLEXIBLE'),
    (400, 500, 2300, 2400, 'OBLIGATE'),
    (550, 600, 950, 1000, 'OPTIONAL_OBLIGATE'),
    (840, 900, 2540, 2600, 'EXCLUSIVE'),
    (1080, 1120, 1420, 1460, 'OPTIONAL_EXCLUSIVE'),
    (1210, 1360, 1750, 1900, 'OPTIONAL_OBLIGATE'),
    # gene_0
    (1, 200, 1300, 1500, 'FLEXIBLE'),
    (100, 200, 600, 700, 'EXCLUSIVE'),
    (600, 700, 1400, 1500, 'FLEXIBLE'),
    (400, 500, 850, 950, 'FLEXIBLE')
]


results_db_path = Path("mock_results.db")
if results_db_path.exists():
    os.remove("mock_results.db")

exonize_obj = Exonize(
        gff_file_path=Path('mock_gff.gff3'),
        genome_file_path=Path('mock_genome.fa'),
        output_prefix="mock_specie",
        draw_event_multigraphs=False,
        csv=False,
        enable_debug=False,
        soft_force=False,
        evalue_threshold=0.01,
        sleep_max_seconds=0,
        min_exon_length=30,
        self_hit_threshold=0.5,
        cpus_number=1,
        timeout_database=60,
        hard_force=False,
        cds_overlapping_threshold=0.9,
        query_overlapping_threshold=0.9,
        output_directory_path=Path("."),
    )
shutil.rmtree("mock_specie_exonize", ignore_errors=True)
exonize_obj.database_interface.results_database_path = results_db_path
exonize_obj.database_interface.connect_create_results_database()
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
    query_overlap_threshold=exonize_obj.query_overlapping_threshold,
    evalue_threshold=exonize_obj.evalue_threshold
        )
exonize_obj.events_reconciliation()
exonize_obj.transcript_interdependence_classification()


def test_representative_cdss():
    gene_0_rcs = exonize_obj.blast_engine.get_candidate_cds_coordinates('gene_0')['candidates_cds_coordinates']
    gene_1_rcs = exonize_obj.blast_engine.get_candidate_cds_coordinates('gene_1')['candidates_cds_coordinates']
    assert set(gene_0_rcs) == representative_cds_gene_0
    assert set(gene_1_rcs) == representative_cds_gene_1


def test_expansion():
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

    with sqlite3.connect(results_db_path) as db:
        cursor = db.cursor()
        cursor.execute(
            """
            SELECT
             QueryExonStart,
             QueryExonEnd,
             TargetStart,
             TargetEnd,
             Classification
            FROM Matches_full_length_non_reciprocal
            WHERE Mode="FULL" or Mode="INSERTION_EXCISION";
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


gene_hierarchy_dictionary_expansions_test = {'gene1': {
    'mRNAs': {
        'tran1': {
            'structure': [
                {'coordinate': P.open(0, 100), 'type': 'CDS'},
                {'coordinate': P.open(150, 250), 'type': 'CDS'},
                {'coordinate': P.open(300, 500), 'type': 'CDS'},
                {'coordinate': P.open(600, 700), 'type': 'CDS'},
                {'coordinate': P.open(900, 1000), 'type': 'CDS'},
                {'coordinate': P.open(1100, 1300), 'type': 'CDS'},
            ]
        },
        'tran2': {
            'structure': [
                {'coordinate': P.open(0, 100), 'type': 'CDS'},
                {'coordinate': P.open(150, 250), 'type': 'CDS'},
                {'coordinate': P.open(600, 700), 'type': 'CDS'},
                {'coordinate': P.open(900, 1000), 'type': 'CDS'},
                {'coordinate': P.open(1100, 1300), 'type': 'CDS'},
            ]
        },
        'tran3': {
            'structure': [
                {'coordinate': P.open(0, 100), 'type': 'CDS'},
                {'coordinate': P.open(600, 700), 'type': 'CDS'},
                {'coordinate': P.open(900, 1000), 'type': 'CDS'},
                {'coordinate': P.open(1100, 1300), 'type': 'CDS'},
                {'coordinate': P.open(1400, 1500), 'type': 'CDS'},
            ]
        }
    }
},
    'gene2': {
        'mRNAs': {
            'tran1': {
                'structure': [
                    {'coordinate': P.open(0, 100), 'type': 'CDS'},
                    {'coordinate': P.open(150, 250), 'type': 'CDS'},
                    {'coordinate': P.open(300, 500), 'type': 'CDS'},
                    {'coordinate': P.open(1650, 1750), 'type': 'CDS'},
                    {'coordinate': P.open(2100, 2200), 'type': 'CDS'},
                    {'coordinate': P.open(2400, 2500), 'type': 'CDS'},
                    {'coordinate': P.open(2900, 3000), 'type': 'CDS'},
                    {'coordinate': P.open(3100, 3200), 'type': 'CDS'},
                    {'coordinate': P.open(3400, 3500), 'type': 'CDS'},
                    {'coordinate': P.open(4000, 4100), 'type': 'CDS'},
                    {'coordinate': P.open(4200, 4300), 'type': 'CDS'},

                ]
            },
            'tran2': {
                'structure': [
                    {'coordinate': P.open(0, 100), 'type': 'CDS'},
                    {'coordinate': P.open(150, 250), 'type': 'CDS'},
                    {'coordinate': P.open(2900, 3000), 'type': 'CDS'},
                    {'coordinate': P.open(3100, 3200), 'type': 'CDS'},
                    {'coordinate': P.open(3400, 3500), 'type': 'CDS'},
                    {'coordinate': P.open(4000, 4100), 'type': 'CDS'},
                    {'coordinate': P.open(4500, 4600), 'type': 'CDS'},

                ]
            },
            'tran3': {
                'structure': [
                    {'coordinate': P.open(600, 700), 'type': 'CDS'},
                    {'coordinate': P.open(900, 1000), 'type': 'CDS'},
                    {'coordinate': P.open(1200, 1300), 'type': 'CDS'},
                    {'coordinate': P.open(1400, 1500), 'type': 'CDS'},
                    {'coordinate': P.open(1650, 1750), 'type': 'CDS'},
                    {'coordinate': P.open(2700, 2800), 'type': 'CDS'},
                    {'coordinate': P.open(4200, 4300), 'type': 'CDS'},
                    {'coordinate': P.open(4500, 4600), 'type': 'CDS'},
                ]
            }
        }
    }}

test_events = [
    # FLEXIBLE
    (1, 'gene1', "FULL", 0, 100, 2, None, 0),
    (2, 'gene1', "INSERTION_EXCISION", 300, 400, 1, None, 0),
    (3, 'gene1', "INACTIVE_UNANNOTATED", 700, 800, 1, None, 0),

    # FLEXIBLE
    (4, 'gene1', "FULL", 0, 100, 3, None, 1),
    (5, 'gene1', "FULL", 150, 250, 3, None, 1),
    (6, 'gene1', "INSERTION_EXCISION", 300, 400, 2, None, 1),
    (7, 'gene1', "INACTIVE_UNANNOTATED", 700, 800, 2, None, 1),

    # OBLIGATE
    (8, 'gene1', "FULL", 600, 700, 2, None, 2),
    (9, 'gene1', "FULL", 900, 1000, 2, None, 2),
    (10, 'gene1', "INSERTION_EXCISION", 1100, 1200, 2, None, 2),

    # EXCLUSIVE
    (11, 'gene2', "FULL", 0, 100, 2, None, 0),
    (12, 'gene2', "FULL", 150, 250, 2, None, 0),
    (13, 'gene2', "FULL", 600, 700, 2, None, 0),

    # EXCLUSIVE
    (14, 'gene2', "FULL", 4000, 4100, 2, None, 1),
    (15, 'gene2', "FULL", 4200, 4300, 2, None, 1),
    (16, 'gene2', "FULL", 4500, 4600, 2, None, 1),

    # OPTIONAL - FLEXIBLE
    (17, 'gene2', "FULL", 1200, 1300, 2, None, 2),
    (18, 'gene2', "FULL", 1400, 1500, 2, None, 2),
    (19, 'gene2', "FULL", 1650, 1750, 2, None, 2),

    # OPTIONAL - EXCLUSIVE
    (20, 'gene2', "FULL", 2100, 2200, 2, None, 3),
    (21, 'gene2', "FULL", 2400, 2500, 2, None, 3),
    (22, 'gene2', "FULL", 2700, 2800, 2, None, 3),

    # OPTIONAL - OBLIGATE
    (23, 'gene2', "FULL", 2900, 3000, 2, None, 4),
    (24, 'gene2', "FULL", 3100, 3200, 2, None, 4),
    (25, 'gene2', "FULL", 3400, 3500, 2, None, 4),

    # NON-CODING
    (26, 'gene1', "FULL", 300, 500, 1, None, 3),
    (27, 'gene1', "INACTIVE_UNANNOTATED", 1000, 1200, 1, None, 3),

]

expected_expansions_classification = [
    ('gene1', 0, 3, 2, 2, 2, 2, 0, 'FLEXIBLE', ''),  # n x (k + 1) = 3 x (1 + 1) = 6
    ('gene1', 1, 3, 3, 3, 3, 3, 0, 'FLEXIBLE', ''),
    ('gene1', 2, 3, 3, 9, 0, 0, 0, 'OBLIGATE', ''),
    ('gene2', 0, 3, 3, 0, 5, 4, 0, 'EXCLUSIVE',
     '_'.join([str(i) for i in
               (P.open(600, 700), tuple((P.open(0, 100), P.open(150, 250))))
               ])),
    ('gene2', 1, 3, 3, 0, 6, 3, 0, 'EXCLUSIVE',
     '_'.join([str(i) for i in
               (P.open(4500, 4600), P.open(4200, 4300), P.open(4000, 4100))
               ])),
    ('gene2', 2, 3, 3, 3, 1, 2, 3, 'OPTIONAL_FLEXIBLE', ''),
    ('gene2', 3, 3, 3, 0, 3, 3, 3, 'OPTIONAL_EXCLUSIVE',
     '_'.join([str(i) for i
               in (P.open(2700, 2800), (P.open(2100, 2200), P.open(2400, 2500)))
               ])),
    ('gene2', 4, 3, 3, 6, 0, 0, 3, 'OPTIONAL_OBLIGATE', ''),
]
# there's missing a test for optional, where all > 0 and intersection is empty


def test_expansion_transcript_iterdependence_classification():
    exonize_obj.data_container.gene_hierarchy_dictionary = gene_hierarchy_dictionary_expansions_test
    expansions_dict = defaultdict(lambda: defaultdict(lambda: list()))
    for event in test_events:
        matchid, geneid, mode, start, end, degree, clusterid, expansionid = event
        if mode in ['FULL', 'INSERTION_EXCISION']:
            expansions_dict[geneid][expansionid].append(P.open(start, end))
    expansion_interdependence_tuples = exonize_obj.event_classifier.classify_expansion_interdependence(
        expansions_dictionary=expansions_dict
    )
    assert expansion_interdependence_tuples == expected_expansions_classification
