# WGS-Pipeline

A reproducible, modular, and fully automated whole-genome sequencing (WGS) variant-calling workflow built with [Nextflow](https://www.nextflow.io/).  
This pipeline performs high-throughput joint genotyping of multiple samples using GATK best-practice tools, with built-in support for HPC environments (Slurm) and customizable configurations for any reference genome or organism.

---

## Overview

The WGS-Pipeline implements a complete variant discovery workflow, including quality control, read preprocessing, alignment, variant calling, joint genotyping, filtering, and functional annotation.  
It is suitable for large-scale population genomic analyses or cohort-based studies involving hundreds of WGS samples.

---

## Key Features

- **End-to-end processing:** from raw FASTQ files to annotated and filtered VCFs.
- **Fully modular design:** each step is implemented as an independent Nextflow module located under `modules/`.
- **Scalable and reproducible:**
  - Supports execution on local or HPC environments.
  - Built-in profiles for Slurm clusters and local testing.
  - Automatically manages job submission, retries, and resource usage.
- **Custom genome support:**
  - Accepts user-provided reference assemblies and custom [snpEff](https://pcingola.github.io/SnpEff/) databases.
  - Compatible with non-model organisms and draft assemblies.
- **Optimized performance:**
  - Configurable scatter-gather intervals for parallel execution.
  - Separate workflows for SNPs and INDELs, combined at the end.
- **Robust fault tolerance:**
  - Includes automatic retries and queue throttling to prevent job submission failures on busy clusters.

---

## Prerequisites

- **Software**
  - [Nextflow](https://www.nextflow.io/) (v21.04.0 or later)
  - Java 11 or higher
  - [Docker](https://www.docker.com/) or [Apptainer/Singularity](https://apptainer.org/)
  - GATK (bundled via container)
  - bcftools, tabix, and snpEff (bundled or custom)
- **Input files**
  - Reference genome in FASTA format, with:
    - Corresponding index file (`.fai`)
    - Sequence dictionary (`.dict`)
  - Pre-built snpEff database for the target genome
  - Paired-end FASTQ files or a tab-delimited sample sheet (`samplesheet.tsv`)
- **Compute environment**
  - Sufficient disk space for BAM and GVCF intermediates (hundreds of GBs per sample)
  - HPC with Slurm scheduler for large datasets

---

## Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/chesterzcl/wgs-pipeline.git
cd wgs-pipeline

### 2. Configure the pipeline parameters

Edit nextflow.config and other configs files within conf folder.

### 3. Launch with preconfigured scripts

For local testing, use scripts/run_nf_local.sh

For execution on Slurm managed computing clusters, use scripts/run_nf_slurm.sh


