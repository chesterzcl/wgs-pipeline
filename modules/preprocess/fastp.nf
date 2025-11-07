process FASTP {
  tag "${sample}"
  publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: '*.trim.fastq.gz'
  
  input:
    tuple val(sample), val(library), path(fqs), val(rg)

  output:
    tuple val(sample), val(library), path("*.trim.fastq.gz"), val(rg)

  script:
  def r1 = fqs[0]
  def r2 = fqs[1]

  """
  fastp \
    -i ${r1} -I ${r2} \
    -o ${sample}.${library}.R1.trim.fastq.gz \
    -O ${sample}.${library}.R2.trim.fastq.gz \
    --detect_adapter_for_pe \
    --trim_front1 10 --trim_front2 10 \
    --cut_right --cut_right_window_size 10 --cut_right_mean_quality 30 --length_required 50 \
    --trim_poly_g \
    --thread ${task.cpus} \
    --html ${sample}.${library}.fastp.html \
    --json ${sample}.${library}.fastp.json
  """
}
