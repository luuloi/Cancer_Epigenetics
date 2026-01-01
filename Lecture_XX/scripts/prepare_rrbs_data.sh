#!/usr/bin/env bash
set -euo pipefail

# Create data directories
cd "dna_meth_analysis"
mkdir -p "./data"
mkdir -p "./data/tmp"

# Download hg19 reference genome
wget -O \
    "./refs/hg19.fa.gz" \
    "https://hgdownload.gi.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
## Index the reference genome (in bgzipped format)
samtools faidx "./refs/hg19.fa.gz"

# Bismark Genome Preparation
bismark_genome_preparation \
    --verbose \
    --bowtie2 \
    ./refs/

# Download raw RRBS fastq files from SRA
cut -d "," -f 2 data/samplesheet.csv | sed '1d' | xargs -n 1 -P 2 -I {} \
    fasterq-dump {} \
        --split-files \
        --threads 4 \
        --progress \
        --outdir "./data" \
        --temp "./data/tmp/{}"

# gzip the fastq files
pigz -p 6 ./data/*.fastq
