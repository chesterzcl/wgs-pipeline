process SAMTOOLS_POSTPROC {
  tag "${sample}"
  publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.q30.bam*'

  input:
    tuple val(sample), path(samfile)

  output:
    tuple val(sample), path("${sample}.q30.bam"), path("${sample}.q30.bam.bai")

  script:
  """
  samtools view -@ ${task.cpus} -bS ${samfile} -o ${sample}.unsorted.bam
  samtools sort -@ ${task.cpus} -n -o ${sample}.name.bam ${sample}.unsorted.bam
  rm -f ${samfile} ${sample}.unsorted.bam

  samtools fixmate -@ ${task.cpus} -m ${sample}.name.bam ${sample}.fixmate.bam
  rm -f ${sample}.name.bam

  samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam ${sample}.fixmate.bam
  rm -f ${sample}.fixmate.bam

  samtools markdup -@ ${task.cpus} -r ${sample}.sorted.bam ${sample}.dedup.bam
  rm -f ${sample}.sorted.bam

  samtools view -@ ${task.cpus} -b -q 30 -o ${sample}.q30.bam ${sample}.dedup.bam
  rm -f ${sample}.dedup.bam

  samtools index -@ ${task.cpus} ${sample}.q30.bam
  """
}