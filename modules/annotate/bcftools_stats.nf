process BCFTOOLS_STATS {
  tag "${type}:${vcf.simpleName}"
  publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.vcfstats.txt'

  input:
    tuple val(type), path(vcf)

  output:
    tuple val(type), path("${vcf.simpleName}.${type}.vcfstats.txt"), emit: stats

  script:
  """
  bcftools stats -F ${params.fasta} ${vcf} > ${vcf.simpleName}.${type}.vcfstats.txt
  """
}
