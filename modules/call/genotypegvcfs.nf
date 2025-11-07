process GENOTYPEGVCFS {
  tag "${interval}"
  publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.raw.vcf.gz*'

  input:
    tuple val(interval), path(workspace)

  output:
    tuple val(interval),
          path("${params.cohort_name}.${interval}.raw.vcf.gz"),
          path("${params.cohort_name}.${interval}.raw.vcf.gz.tbi")

  script:
  """
  gatk GenotypeGVCFs \
    -R ${params.fasta} \
    -V gendb://${workspace} \
    -L ${interval} \
    -O ${params.cohort_name}.${interval}.raw.vcf.gz

  gatk IndexFeatureFile -I ${params.cohort_name}.${interval}.raw.vcf.gz

  """
}
