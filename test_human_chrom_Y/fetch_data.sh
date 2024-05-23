#!/bin/bash

curl -O https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
curl -O https://ftp.ensembl.org/pub/release-111/gff3/homo_sapiens/Homo_sapiens.GRCh38.112.chromosome.Y.gff3.gz
gunzip *gff3.gz