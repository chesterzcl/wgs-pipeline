process BCFTOOLS_STATS {
  tag "${vcf.simpleName}"
  publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.vcfstats.txt'

  input:
    path vcf

  output:
    path "${vcf.simpleName}.vcfstats.txt", emit: stats

  script:
  """
  bcftools stats -F ${params.fasta} ${vcf} > ${vcf.simpleName}.vcfstats.txt
  """
}