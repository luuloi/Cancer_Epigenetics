#!/usr/bin/env bash
set -euo pipefail

mkdir -p ./trimmed/

# trim_galore \
#     --rrbs \
#     --paired \
#     --illumina \
#     --gzip \
#     --fastqc \
#     --output_dir ./trimmed/ \
#     --cores 4 \
#     "./data/SRR5311114_1.fastq.gz" "./data/SRR5311114_2.fastq.gz"

for r1 in ./data/*_1.fastq.gz; do
    # Extract sample ID from filename and set up paths for R2
    sample_basename="$(basename "${r1}")"
    sample_id="${sample_basename%_1.fastq.gz}"
    r2="./data/${sample_id}_2.fastq.gz"

    # Create output directory for each trimmed sample
    outdir="./trimmed/${sample_id}"
    mkdir -p "${outdir}"

    # Run Trim Galore for RRBS data
    echo "Trimming ${sample_id}..."
    trim_galore \
        --rrbs \
        --paired \
        --illumina \
        --gzip \
        --fastqc \
        --output_dir "${outdir}" \
        --cores 4 \
        "$r1" "$r2"
done

echo "All done."