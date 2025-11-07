process CREATE_DICT {
  tag "${fasta.simpleName}"
  publishDir "${params.outdir}/ref", mode: 'copy', pattern: '*.dict'

  input:
    path fasta

  output:
    path "${fasta.simpleName}.dict", emit: dict

  """
  gatk CreateSequenceDictionary \
    -R ${fasta} \
    -O ${fasta.simpleName}.dict
  """
}

