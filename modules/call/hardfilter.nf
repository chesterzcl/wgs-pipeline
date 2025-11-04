process HARD_FILTER {
  tag "GATK_HardFilter"
  publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.filt.vcf.gz*'

  input:
    path vcf

  output:
    path 'snps.filt.vcf.gz'   , emit: snps
    path 'indels.filt.vcf.gz' , emit: indels

  script:
  """
  gatk SelectVariants -R ${params.fasta} -V ${vcf} -select-type SNP   -O snps.vcf.gz
  gatk SelectVariants -R ${params.fasta} -V ${vcf} -select-type INDEL -O indels.vcf.gz

  gatk VariantFiltration -R ${params.fasta} -V snps.vcf.gz -O snps.filt.vcf.gz \
    --filter-name QD2        --filter-expression 'QD < 2.0' \
    --filter-name FS60       --filter-expression 'FS > 60.0' \
    --filter-name MQ40       --filter-expression 'MQ < 40.0' \
    --filter-name SOR3       --filter-expression 'SOR > 3.0' \
    --filter-name MQRankSum  --filter-expression 'MQRankSum < -12.5' \
    --filter-name ReadPosRS  --filter-expression 'ReadPosRankSum < -8.0'

  gatk VariantFiltration -R ${params.fasta} -V indels.vcf.gz -O indels.filt.vcf.gz \
    --filter-name QD2        --filter-expression 'QD < 2.0' \
    --filter-name FS200      --filter-expression 'FS > 200.0' \
    --filter-name SOR10      --filter-expression 'SOR > 10.0' \
    --filter-name ReadPosRS  --filter-expression 'ReadPosRankSum < -20.0'

  tabix -p vcf snps.filt.vcf.gz
  tabix -p vcf indels.filt.vcf.gz
  """
}