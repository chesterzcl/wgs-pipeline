process FASTQC {
  tag "${sample}"
  publishDir "${params.outdir}/qc/fastqc", mode: 'copy', pattern: '*_fastqc.*'

  input:
    tuple val(sample), val(lib), path(fqs)

  output:
    path "*_fastqc.zip",  emit: zip
    path "*_fastqc.html", emit: html

  script:
  // Join the list of input FASTQs into a space-separated string
  def inputs = (fqs instanceof List) ? fqs.join(' ') : fqs.toString()
  """
  fastqc -o . ${inputs}
  """
}