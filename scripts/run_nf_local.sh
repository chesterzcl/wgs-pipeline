#!/usr/bin/env bash
set -euo pipefail

# ========== USER PARAMETERS ==========
PROFILE="docker"
MAIN_NF="wgs-pipeline/main.nf"

# Absolute paths recommended
SAMPLESHEET="wgs-pipeline/samplesheet"
FASTA="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2/primary_dedup_chr_masked_hp_sealed.fa"
OUTDIR="results_test"

# Optional flags
RESUME="-resume"
REPORTS="-with-report ${OUTDIR}/run_report.html -with-timeline ${OUTDIR}/run_timeline.html"
# ====================================

# Create output directory
mkdir -p "${OUTDIR}"

# Auto-overwrite old report files
rm -f "${OUTDIR}/run_report.html" "${OUTDIR}/run_timeline.html"

# Make sure FASTA index exists (bcftools / samtools need .fai)
if [[ ! -f "${FASTA}.fai" ]]; then
  echo "Indexing FASTA..."
  docker run --rm -v "$(dirname "${FASTA}"):/ref" quay.io/biocontainers/samtools:1.20--h50ea8bc_1 \
    samtools faidx "/ref/$(basename "${FASTA}")"
fi

# Run Nextflow pipeline
./nextflow run "${MAIN_NF}" \
  -profile "${PROFILE}" \
  --samplesheet "${SAMPLESHEET}" \
  --fasta "${FASTA}" \
  --outdir "${OUTDIR}" \
  ${RESUME} \
  ${REPORTS} \

