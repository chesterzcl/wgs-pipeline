process HAPLOTYPECALLER {
  tag "${sample}"
  publishDir "${params.outdir}/gvcf", mode: 'copy', pattern: '*.g.vcf.gz*'

  input:
    tuple val(sample), path(bam), path(bai)

  output:
    tuple val(sample), path("${sample}.g.vcf.gz")

  script:
  """
  gatk HaplotypeCaller \
    -R ${params.fasta} \
    -I ${bam} \
    -O ${sample}.g.vcf.gz \
    -ERC GVCF
  """
}
