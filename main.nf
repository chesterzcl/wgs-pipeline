nextflow.enable.dsl = 2

// ---------- Params ----------
params.samplesheet = params.samplesheet ?: 'samplesheet'   // your file has no .csv
params.outdir      = params.outdir ?: 'results'
params.lib_default = params.lib_default ?: 'lib1'

// ---------- Includes (can be top-level) ----------
include { FAIDX }       from './modules/ref/faidx.nf'
include { CREATE_DICT } from './modules/ref/create_dict.nf'
include { BWA_INDEX }   from './modules/ref/bwa_index.nf' 
include { PARTITION_GENOME } from './modules/ref/partition_genome.nf'
include { FASTQC as FASTQC_RAW }  from './modules/qc/fastqc.nf'
include { FASTQC as FASTQC_TRIM } from './modules/qc/fastqc.nf'
include { FASTP  }                from './modules/preprocess/fastp.nf'
include { BWA_MEM }           from './modules/align/bwa_mem.nf'
include { SAMTOOLS_POSTPROC } from './modules/postprocess/samtools_postproc.nf'
include { HAPLOTYPECALLER }   from './modules/call/haplotypecaller_gvcf.nf'
include { GENOMICSDB }     from './modules/call/genomicsdb.nf'
include { GENOTYPEGVCFS }     from './modules/call/genotypegvcfs.nf'
include { HARD_FILTER }       from './modules/call/hardfilter.nf'
include { ORDER_VCFS as ORDER_VCFS_SNPS } from './modules/call/order_vcfs.nf'
include { ORDER_VCFS as ORDER_VCFS_INDELS } from './modules/call/order_vcfs.nf'
include { SNPEFF as SNPEFF_SNPS }   from './modules/annotate/snpeff.nf'
include { SNPEFF as SNPEFF_INDELS } from './modules/annotate/snpeff.nf'
include { CONCAT_VCFS as CONCAT_VCFS_SNPS } from './modules/call/concat_vcfs.nf'
include { CONCAT_VCFS as CONCAT_VCFS_INDELS } from './modules/call/concat_vcfs.nf'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_SNPS }    from './modules/annotate/bcftools_stats.nf'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_INDELS }    from './modules/annotate/bcftools_stats.nf'


// ---------- Workflow (all process invocations must be inside) ----------
workflow {

  REF_FASTA = Channel.value( file(params.fasta) )
  REF_FAI   = FAIDX(REF_FASTA).fai
  REF_DICT  = CREATE_DICT(REF_FASTA).dict
  REF_BWA   = BWA_INDEX(REF_FASTA).index 
  INTERVAL_FILE = PARTITION_GENOME(REF_FAI).list
  INTVL_CH = INTERVAL_FILE.splitText().map{ it.trim() }.filter{ it }
  
  // Samplesheet READS
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


  // QC raw, trim, QC trimmed
  QC_RAW  = READS_LIB_RG.map { s, lib, fqs, rg -> tuple(s, lib, fqs) } | FASTQC_RAW
  TRIMMED = FASTP(READS_LIB_RG)
  QC_TRIM = TRIMMED.map      { s, lib, fqs, rg -> tuple(s, lib, fqs) } | FASTQC_TRIM

  // Align postprocess (fixmate/sort/rmdup/MAPQ30/index)
  ALIGNED_SAM   = BWA_MEM(TRIMMED)
  POSTPROC    = SAMTOOLS_POSTPROC(ALIGNED_SAM)
  GVCFS          = HAPLOTYPECALLER(POSTPROC)

  // Per-sample GVCF, then joint genotyping
  SAMPLE_MAP = GVCFS.collectFile { id, gvcf ->
    [ "${params.cohort_name}_map.tsv", "${id}\t${gvcf}\n" ]
  } 
  
  SAMPLE_MAP_VAL = SAMPLE_MAP.first()
  GDB_WS = GENOMICSDB(SAMPLE_MAP_VAL, INTVL_CH)
  JOINT_RAW = GENOTYPEGVCFS(GDB_WS)
  
  HARD = HARD_FILTER(JOINT_RAW)

  SNPS   = HARD.snps
  INDELS = HARD.indels
  
  SNPS_LIST_FILE = SNPS.toList().map { rows ->
     def f = file("snps.list")
     f.text = rows.collect { ivl, vcf -> "${ivl}\t${vcf}" }.join('\n') + '\n'
     f
  } 

  INDELS_LIST_FILE = INDELS.toList().map { rows ->
     def f = file("indels.list")
     f.text = rows.collect { ivl, vcf -> "${ivl}\t${vcf}" }.join('\n') + '\n'
     f
  }

  SNPS_ORDERED_LIST = ORDER_VCFS_SNPS( INTERVAL_FILE , SNPS_LIST_FILE )
  INDELS_ORDERED_LIST = ORDER_VCFS_INDELS( INTERVAL_FILE , INDELS_LIST_FILE )

  MERGED_SNPS   = CONCAT_VCFS_SNPS( SNPS_ORDERED_LIST.map   { lf -> tuple('snps',   lf) } )
  MERGED_INDELS = CONCAT_VCFS_INDELS( INDELS_ORDERED_LIST.map { lf -> tuple('indels', lf) } )

  ANN_SNPS   = SNPEFF_SNPS( MERGED_SNPS ).vcf
  ANN_INDELS = SNPEFF_INDELS( MERGED_INDELS ).vcf
  
  ANN_SNPS_VCF   = ANN_SNPS.map   { type, vcf, tbi -> tuple(type, vcf) }
  ANN_INDELS_VCF = ANN_INDELS.map { type, vcf, tbi -> tuple(type, vcf) }

  // stats + multiqc (concat works fine here too)
  VCF_STATS_SNPS = BCFTOOLS_STATS_SNPS(ANN_SNPS_VCF).stats
  VCF_STATS_INDELS = BCFTOOLS_STATS_INDELS(ANN_INDELS_VCF).stats

  // Emits
  //emit:
    //FILTERED_VCF
    //annotated_vcf  = Channel.concat(ANN_SNPS, ANN_INDELS)
    //vcf_stats      = Channel.concat(VCF_STATS_SNPS, VCF_STATS_INDELS)
}
