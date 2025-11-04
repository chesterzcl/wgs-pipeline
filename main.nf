nextflow.enable.dsl = 2

// ---------- Params ----------
params.samplesheet = params.samplesheet ?: 'samplesheet'   // your file has no .csv
params.outdir      = params.outdir ?: 'results'
params.lib_default = params.lib_default ?: 'lib1'

// ---------- Includes (can be top-level) ----------
include { FASTQC as FASTQC_RAW }  from './modules/qc/fastqc.nf'
include { FASTQC as FASTQC_TRIM } from './modules/qc/fastqc.nf'
include { FASTP  }                from './modules/preprocess/fastp.nf'
include { BWA_MEM }           from './modules/align/bwa_mem.nf'
include { SAMTOOLS_POSTPROC } from './modules/postprocess/samtools_postproc.nf'
include { HAPLOTYPECALLER }   from './modules/call/haplotypecaller_gvcf.nf'
include { GENOTYPEGVCFS }     from './modules/call/genotypegvcfs.nf'
include { HARD_FILTER }       from './modules/call/hardfilter.nf'
include { SNPEFF as SNPEFF_SNPS }   from './modules/annotate/snpeff.nf'
include { SNPEFF as SNPEFF_INDELS } from './modules/annotate/snpeff.nf'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_SNPS }    from './modules/annotate/bcftools_stats.nf'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_INDELS }    from './modules/annotate/bcftools_stats.nf'


// ---------- Workflow (all process invocations must be inside) ----------
workflow {

  // Samplesheet → READS
  Channel
    .fromPath(params.samplesheet, checkIfExists: true)
    .ifEmpty { error "Samplesheet not found or empty: ${params.samplesheet}" }
    .splitCsv(header: true)
    .map { row ->
      tuple( row.sample as String, [ file(row.fastq1), file(row.fastq2) ] )
    }
    .set { READS }   // (sample, [R1, R2])

  // Add library + read-group
  def defaultLib = params.lib_default

  READS
    .map { sample, fqs ->
      def lib = defaultLib
      def rg  = "@RG\\tID:${sample}\\tSM:${sample}\\tLB:${lib}\\tPL:ILLUMINA"
      tuple(sample, lib, fqs, rg)
    }
    .set { READS_LIB_RG }

  READS.view()
  READS_LIB_RG.view()

  // QC raw, trim, QC trimmed
  QC_RAW  = READS_LIB_RG.map { s, lib, fqs, rg -> tuple(s, lib, fqs) } | FASTQC_RAW
  TRIMMED = FASTP(READS_LIB_RG)
  QC_TRIM = TRIMMED.map      { s, lib, fqs, rg -> tuple(s, lib, fqs) } | FASTQC_TRIM

  // Align → postprocess (fixmate → sort → rmdup → MAPQ≥30 → index)
  ALIGNED_SAM = BWA_MEM(TRIMMED)
  POSTPROC    = SAMTOOLS_POSTPROC(ALIGNED_SAM)

  // No BQSR — feed postprocessed BAMs directly to HC
  FINAL_BAM = POSTPROC

  // Per-sample GVCF, then joint genotyping
  GVCFS     = HAPLOTYPECALLER(FINAL_BAM)
  ALL_GVCFS = GVCFS.map { sample,gvcf->gvcf }.collect()
  JOINT_RAW = GENOTYPEGVCFS(ALL_GVCFS)

  // Hard filter → two VCFs (SNPs & INDELs)
  HARD = HARD_FILTER(JOINT_RAW)

  SNPS   = HARD.snps
  INDELS = HARD.indels

  // Annotate (pick the single labeled stream `.vcf` from each)
  ANN_SNPS   = SNPEFF_SNPS(SNPS).vcf
  ANN_INDELS = SNPEFF_INDELS(INDELS).vcf

  // stats + multiqc (concat works fine here too)
  VCF_STATS_SNPS = BCFTOOLS_STATS_SNPS(ANN_SNPS).stats
  VCF_STATS_INDELS = BCFTOOLS_STATS_INDELS(ANN_INDELS).stats

  // Emits
  emit:
    filtered_vcf   = Channel.concat(SNPS, INDELS)
    annotated_vcf  = Channel.concat(ANN_SNPS, ANN_INDELS)
    vcf_stats      = Channel.concat(VCF_STATS_SNPS, VCF_STATS_INDELS)
}