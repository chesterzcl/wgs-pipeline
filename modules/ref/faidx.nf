process FAIDX {
  tag "${fasta.simpleName}"
  publishDir "${params.outdir}/ref", mode: 'copy', pattern: '*.fai'

  input:
    path fasta

  output:
    path "${fasta}.fai", emit: fai

  """
  samtools faidx ${fasta}
  """
}
