#!/usr/bin/env bash
set -euo pipefail

# Bismark alignment for RRBS data
mkdir -p ./bismark/

for r1 in ./trimmed/*/*_val_1.fq.gz; do
    # Extract sample ID from filename and set up paths for R2
    sample_basename="$(basename "${r1}")"
    sample_id="${sample_basename%_1_val_1.fq.gz}"
    r2="./trimmed/${sample_id}/${sample_id}_2_val_2.fq.gz"

    # Create output directory for each sample
    outdir="./bismark/${sample_id}"
    mkdir -p "${outdir}"

    # Run Bismark alignment
    echo "Aligning ${sample_id} with Bismark..."
    bismark \
        --genome ./refs/ \
        --fastq \
        --bowtie2 \
        --output_dir "${outdir}" \
        --parallel 4 \
        -N 0 \
        --gzip \
        -1 "${r1}" -2 "${r2}"
done

# Sort and index the BAM files
for bam in ./bismark/*/*_pe.bam; do
    sorted_bam="${bam%.bam}.sorted.bam"
    echo "Sorting ${bam} to ${sorted_bam}..."
    samtools sort -o "${sorted_bam}" "${bam}"
    echo "Indexing ${sorted_bam}..."
    samtools index "${sorted_bam}"
done

# Keep proper pairs on chr16, then downsample to 10% (pair-safe)
for sorted_bam in ./bismark/*/*_pe.sorted.bam; do
    sample_basename="$(basename "${sorted_bam}")"
    sample_id="${sample_basename%_1_val_1_bismark_bt2_pe.sorted.bam}"
    outdir="./bismark/${sample_id}"

    echo "Processing ${sample_id} for chr16 proper pairs downsampling..."

    # Steps:
    # 1) Subset to chr16, keep only proper pairs and drop secondary/supplementary alignments
    # 2) Queryname-sort so pairs are together
    # 3) Downsample to 10% (keeps pairs together)
    # 4) Coordinate-sort + index

    ## 1) Subset to chr16, keep only proper pairs and drop secondary/supplementary alignments
    chr16_proper_bam="${outdir}/${sample_id}_pe_chr16_proper.bam"
    samtools view \
        -@ 8 \
        -b \
        -f 2 \
        -F 0x900 \
        "${sorted_bam}" chr16 > "${chr16_proper_bam}"

    ## 2) Queryname-sort using Picard (Picard-compatible queryname ordering)
    picard SortSam \
        -I "${chr16_proper_bam}" \
        -O "${outdir}/${sample_id}_pe_chr16_proper_qname.bam" \
        --SORT_ORDER queryname

    ## 3) Downsample to 10% (keeps pairs together)
    picard DownsampleSam \
        -I "${outdir}/${sample_id}_pe_chr16_proper_qname.bam" \
        -O "${outdir}/${sample_id}_pe_chr16_proper_10pct_qname.bam" \
        -P 0.10 \
        --STRATEGY Chained \
        -R 123

    # 4) Coordinate-sort + index
    samtools sort \
        -@ 8 \
        -o "${outdir}/${sample_id}_pe_chr16_proper_10pct.sorted.bam" \
        "${outdir}/${sample_id}_pe_chr16_proper_10pct_qname.bam"
    samtools index -@ 8 "${outdir}/${sample_id}_pe_chr16_proper_10pct.sorted.bam"

done

# Convert BAM to paired-FASTQ file
for final_bam in ./bismark/*/*_pe_chr16_proper_10pct.sorted.bam; do
    sample_basename="$(basename "${final_bam}")"
    sample_id="${sample_basename%_pe_chr16_proper_10pct.sorted.bam}"
    outdir="./bismark/${sample_id}"

    echo "Converting ${sample_id} BAM to paired FASTQ..."

    samtools fastq -@ 8 \
        -1 "${outdir}/${sample_id}_chr16_10pct_R1.fastq.gz" \
        -2 "${outdir}/${sample_id}_chr16_10pct_R2.fastq.gz" \
        -0 /dev/null -s /dev/null \
        -n \
        "${final_bam}"
done

# Move subset FASTQ files to a separate directory
mv ./bismark/*/*_chr16_10pct_*.fastq.gz "./raw_subset"