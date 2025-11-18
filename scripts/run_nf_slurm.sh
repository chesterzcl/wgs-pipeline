#!/bin/bash
#SBATCH --job-name=wgs_pipeline
#SBATCH --time=29-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=long

# ============================================================
# 1. Environment setup
# ============================================================
set -euo pipefail

module load Java

# ============================================================
# 2. User parameters
# ============================================================
DIR="/home/zl436/palmer_scratch/pipe"
PROFILE="slurm"
MAIN_NF=${DIR}/"wgs-pipeline/main.nf"
SAMPLESHEET=${DIR}/"wgs-pipeline/samplesheet"
FASTA="/home/zl436/project/Contigs_dedup_scfd_3ddna_mc_hp/primary_dedup_chr_masked_hp_sealed.fa"
OUTDIR="${DIR}/results_wgs_cv3"

RESUME="-resume"
REPORTS="-with-report ${OUTDIR}/run_report.html -with-timeline ${OUTDIR}/run_timeline.html"

mkdir -p "${OUTDIR}" logs

# ============================================================
# 3. Run pipeline
# ============================================================
./nextflow run "${MAIN_NF}" \
  -profile "${PROFILE}" \
  --samplesheet "${SAMPLESHEET}" \
  --fasta "${FASTA}" \
  --outdir "${OUTDIR}" \
  ${RESUME} \
  ${REPORTS}