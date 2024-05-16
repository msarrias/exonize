#!/bin/bash

GFF_FILE_PATH=Homo_sapiens.GRCh38.111.chromosome.Y.gff3
GENOME_FILE_PATH=Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
SPECIE_ID=Homo_sapiens_chrom_Y
MIN_SEQ_LENGTH=30
OUTPUT_DIRECTORY_PATH=.
exonize $GFF_FILE_PATH \
        $GENOME_FILE_PATH \
        --output_prefix $SPECIE_ID \
        --min-exon-length $MIN_SEQ_LENGTH \
        --multigraphs \
        --hard-force \
        --output-directory-path $OUTPUT_DIRECTORY_PATH
