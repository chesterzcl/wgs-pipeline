process GENOTYPEGVCFS {
  tag "JointGenotyping"
  publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.raw.vcf.gz*'

  input:
    // This expects a *list* of GVCF files (from .collect() in main.nf)
    path gvcfs

  output:
    path "joint.raw.vcf.gz"

  script:
  // Build "-V file1 -V file2 ..." safely in Groovy
  def vopts = gvcfs.collect { f -> "-V ${f}" }.join(' ')
  """
  gatk GenotypeGVCFs \
    -R ${params.fasta} \
    ${vopts} \
    -O joint.raw.vcf.gz

  tabix -p vcf joint.raw.vcf.gz
  """
}