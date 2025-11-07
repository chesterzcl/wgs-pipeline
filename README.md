# WGS-Pipeline

A fully automated whole-genome sequencing (WGS) variant-calling workflow built with [Nextflow](https://www.nextflow.io/) for joint-genotype analysis of multiple samples. Designed for high-throughput use (100s of samples) on HPC/Slurm clusters, yet usable locally.

## Key features

- End-to-end processing: QC → trimming → alignment → post-processing → variant calling → joint genotyping → hard filters → annotation & stats.  
- Modular pipeline built under `modules/` for each major step: e.g., BWA mapping, GATK HaplotypeCaller, GenomicsDB, joint genotyping, snpEff annotation, bcftools stats.  
- Built for custom assemblies/organisms: supports user-provided reference FASTA + annotation database (e.g., your oyster genome).  
- Optimized for HPC:  
  - Slurm executor profile included (`-profile slurm`) with partition/QoS routing and job-queue control.  
  - Local mode (`-profile local`) for testing or small deployments.  
- Strong grouping of streams: SNPs vs INDELs handled and annotated separately, then merged.

## Prerequisites

- [Nextflow](https://www.nextflow.io/) installed (v21+ recommended).  
- Java 11 or higher.  
- On HPC: Slurm cluster with sufficient permissions; Apptainer/Singularity or Docker support.  
- Input files:  
  - Reference genome FASTA + index (.fai), sequence dictionary (.dict).  
  - Annotation for snpEff (custom genome built).  
  - Samplesheet listing input FASTQ/reads (see `samplesheet/`).  
- Sufficient storage and memory: expect large intermediate files (BAM, GVCF) when processing hundreds of samples.

## Quickstart

1. Clone the repo:
   ```bash
   git clone https://github.com/chesterzcl/wgs-pipeline.git
   cd wgs-pipeline

2. Configure the pipeline parameters

   Edit nextflow.config and other configs files within conf folder.

3. Launch with preconfigured scripts

   For local testing, use scripts/run_nf_local.sh

   For execution on Slurm managed computing clusters, use scripts/run_nf_slurm.sh


