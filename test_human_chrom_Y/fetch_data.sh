#!/bin/bash

curl -O https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
curl -O https://ftp.ensembl.org/pub/release-114/gff3/homo_sapiens/Homo_sapiens.GRCh38.114.chromosome.Y.gff3.gz
gunzip Homo_sapiens.GRCh38.114.chromosome.Y.gff3.gz