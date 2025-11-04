process SNPEFF {
  tag "snpEff"
  publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.snpeff.vcf.gz*'

  input:
    path vcf

  output:
    path "${vcf.simpleName}.snpeff.vcf.gz", emit: vcf

  script:
  """
  java -Xmx40g -jar /home/zl436/project/tools/snpEff/snpEff.jar ann -v \
    CV_final ${vcf} \
    > ${vcf.simpleName}.snpeff.vcf \
    -s ${vcf.simpleName}.snpeff.vcf.html

  bgzip -c ${vcf.simpleName}.snpeff.vcf > ${vcf.simpleName}.snpeff.vcf.gz
  tabix -p vcf ${vcf.simpleName}.snpeff.vcf.gz
  """
}