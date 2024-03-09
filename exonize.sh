#!/bin/bash

GFF_FILE_PATH=data/bonobo/bonobo_first1081lines.gff
GENOME_FILE_PATH=data/bonobo/Pan_paniscus.panpan1.1.dna_rm.toplevel.fa
SPECIE_ID=bonobo
MIN_SEQ_LENGTH=20
OUTPUT_DIRECTORY_PATH=.
exonize $GFF_FILE_PATH \
        $GENOME_FILE_PATH \
        $SPECIE_ID \
        -el $MIN_SEQ_LENGTH \
        --hard-force \
        --multigraphs \
        --output-directory-path $OUTPUT_DIRECTORY_PATH