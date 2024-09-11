#!/bin/bash

GFF_FILE_PATH=Homo_sapiens.GRCh38.111.chromosome.Y.gff3
GENOME_FILE_PATH=Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
SPECIE_ID=Homo_sapiens_chrom_Y
EVALUE=1e-3
CDS_ANNOT='CDS'
GENE_ANNOT='gene'
TRANS_ANNOT='mRNA'
MIN_SEQ_LENGTH=30
exonize $GFF_FILE_PATH \
        $GENOME_FILE_PATH \
        --gene-annot-feature $GENE_ANNOT \
        --cds-annot-feature $CDS_ANNOT \
        --transcript-annot-feature $TRANS_ANNOT \
        --output_prefix $SPECIE_ID \
        --min-exon-length $MIN_SEQ_LENGTH \
        --evalue-threshold $EVALUE \
        --csv \

