#!/bin/bash

GFF_FILE_PATH=Homo_sapiens.GRCh38.111.chromosome.Y.gff3
GENOME_FILE_PATH=Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
SPECIE_ID=Homo_sapiens_chrom_Y
EVALUE=1e-3
MIN_SEQ_LENGTH=30
exonize $GFF_FILE_PATH \
        $GENOME_FILE_PATH \
        --output_prefix $SPECIE_ID \
        --min-exon-length $MIN_SEQ_LENGTH \
        --evalue-threshold $EVALUE \
        --multigraphs \
        --csv \
        --cpus_number 5 \

