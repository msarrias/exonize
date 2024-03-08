#!/bin/bash

GFF_FILE_PATH=data/bonobo/bonobo_first1081lines.gff
GENOME_FILE_PATH=data/bonobo/bonobo.fa
SPECIE_ID=bonobo
MIN_SEQ_LENGTH=20
OUTPUT_DIRECTORY_PATH=/Users/marinaherrerasarrias/Desktop
exonize $GFF_FILE_PATH \
        $GENOME_FILE_PATH \
        $SPECIE_ID \
        -el $MIN_SEQ_LENGTH \
        --hard-force \
        --multigraphs \
        --output-directory-path $OUTPUT_DIRECTORY_PATH