process GENOMICSDB {
  tag "${params.cohort_name}:${interval}"
  publishDir "${params.outdir}/gdb", mode: 'copy', pattern: 'gendb_*'

  input:
    path sample_map
    val  interval

  output:
    tuple val(interval), path("gendb_${params.cohort_name}_${interval}")

  script:
  """
  set -euo pipefail

  gatk GenomicsDBImport \
    --sample-name-map ${sample_map} \
    --genomicsdb-workspace-path gendb_${params.cohort_name}_${interval} \
    -L ${interval} \
    --reader-threads ${task.cpus}
  """
}
