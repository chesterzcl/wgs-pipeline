process BWA_MEM {
  tag "${sample}"
  publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.unsorted.sam'

  input:
    tuple val(sample), val(library), path(fqs), val(rg)

  output:
    tuple val(sample), path("${sample}.unsorted.sam")

  script:
  """
  bwa mem -t ${task.cpus} -R '${rg}' ${params.fasta} ${fqs[0]} ${fqs[1]} > ${sample}.unsorted.sam
  """
}