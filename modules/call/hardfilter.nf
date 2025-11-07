process HARD_FILTER {
  tag "GATK_HardFilter"
  publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*_bi*_filtered.vcf*'

  input:
    tuple val(interval), path(vcf),path(vcf_tbi)

  output:
    tuple val(interval), path("${vcf.baseName}_biSNP_filtered.vcf"),   emit: snps
    tuple val(interval), path("${vcf.baseName}_biINDEL_filtered.vcf"), emit: indels  

  script:
  """
  gatk SelectVariants \
    --restrict-alleles-to BIALLELIC \
    --select-type-to-include SNP \
    -R ${params.fasta} \
    -V ${vcf} \
    -O ${vcf.baseName}_biSNP.vcf

  gatk IndexFeatureFile \
    -I ${vcf.baseName}_biSNP.vcf

  gatk VariantFiltration \
    -V ${vcf.baseName}_biSNP.vcf \
    -O ${vcf.baseName}_biSNP_marked.vcf \
    -filter "QD < 2.0"                --filter-name "QD2" \
    -filter "QUAL < 30.0"             --filter-name "QUAL30" \
    -filter "SOR > 3.0"               --filter-name "SOR3" \
    -filter "FS > 60.0"               --filter-name "FS60" \
    -filter "MQ < 40.0"               --filter-name "MQ40" \
    -filter "MQRankSum < -12.5"       --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0"   --filter-name "ReadPosRankSum-8"

  gatk IndexFeatureFile \
    -I ${vcf.baseName}_biSNP_marked.vcf

  gatk SelectVariants \
    --restrict-alleles-to BIALLELIC \
    --exclude-filtered \
    --select-type-to-include SNP \
    -R ${params.fasta} \
    -V ${vcf.baseName}_biSNP_marked.vcf \
    -O ${vcf.baseName}_biSNP_filtered.vcf

  gatk IndexFeatureFile \
    -I ${vcf.baseName}_biSNP_filtered.vcf

  rm -f ${vcf.baseName}_biSNP.vcf ${vcf.baseName}_biSNP_marked.vcf

  gatk SelectVariants \
    --restrict-alleles-to BIALLELIC \
    --select-type-to-include INDEL \
    -R ${params.fasta} \
    -V ${vcf} \
    -O ${vcf.baseName}_biINDEL.vcf

  gatk IndexFeatureFile \
    -I ${vcf.baseName}_biINDEL.vcf

  gatk VariantFiltration \
    -V ${vcf.baseName}_biINDEL.vcf \
    -O ${vcf.baseName}_biINDEL_marked.vcf \
    -filter "QD < 2.0"                --filter-name "QD2" \
    -filter "FS > 200.0"              --filter-name "FS200" \
    -filter "SOR > 10.0"              --filter-name "SOR10" \
    -filter "ReadPosRankSum < -20.0"  --filter-name "ReadPosRankSum-20"

  gatk IndexFeatureFile \
    -I ${vcf.baseName}_biINDEL_marked.vcf

  gatk SelectVariants \
    --restrict-alleles-to BIALLELIC \
    --exclude-filtered \
    --select-type-to-include INDEL \
    -R ${params.fasta} \
    -V ${vcf.baseName}_biINDEL_marked.vcf \
    -O ${vcf.baseName}_biINDEL_filtered.vcf

  gatk IndexFeatureFile \
    -I ${vcf.baseName}_biINDEL_filtered.vcf

  rm -f ${vcf.baseName}_biINDEL.vcf ${vcf.baseName}_biINDEL_marked.vcf

   """
}
