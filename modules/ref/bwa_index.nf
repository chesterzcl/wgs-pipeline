process BWA_INDEX {
  tag "${fasta.simpleName}"
  publishDir "${params.outdir}/ref/bwa", mode: 'copy', pattern: "${fasta.simpleName}.*"

  input:
    path fasta

  output:
    path "${fasta}.bwt", emit: index

  """
  ln -s ${fasta} reference.fa
  bwa index reference.fa
  for ext in amb ann bwt pac sa; do
    mv reference.fa.\${ext} ${fasta}.\${ext}
  done
  """
}
