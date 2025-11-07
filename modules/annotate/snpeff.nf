process SNPEFF {
  tag "${type}:${vcf.simpleName}"
  publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.snpeff.vcf.gz*'

  input:
    tuple val(type), path(vcf)     // or (type, vcf, tbi) if you pass the index too

  output:
    tuple val(type),
          path("*.snpeff.vcf.gz"),
          path("*.snpeff.vcf.gz.tbi"), emit: vcf

  script:
  """
  BN=${vcf.simpleName}
  OUT="\${BN}.snpeff.vcf.gz"

  java -Xmx4g -jar /opt/snpEff/snpEff.jar -v ${params.snpeff_genome} ${vcf} \
    | bgzip -c > "\$OUT"
  tabix -p vcf "\$OUT"
  """
}
