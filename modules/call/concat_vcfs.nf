process CONCAT_VCFS {
  tag "${type}"
  publishDir "${params.outdir}/vcf", mode: 'copy',
             pattern: "${params.cohort_name}.${type}.filtered.merged.vcf.gz*"

  input:
    tuple val(type), path(listfile)

  output:
    tuple val(type),
          path("${params.cohort_name}.${type}.filtered.merged.vcf.gz"),
          path("${params.cohort_name}.${type}.filtered.merged.vcf.gz.tbi")

  script:
  """
  OUT="${params.cohort_name}.${type}.filtered.merged.vcf.gz"
  [ -s "${listfile}" ] || { echo "[ERROR] empty listfile"; exit 1; }

  echo "[INFO] First 15 inputs (pre-ordered):"
  nl -ba "${listfile}" | head -n 15

  # Build -I args *without* command substitution
  mapfile -t FILES < "${listfile}"
  ARGS=()
  for f in "\${FILES[@]}"; do
    [[ -n "\$f" ]] && ARGS+=("-I" "\$f")
  done
  echo "gatk GatherVcfs \${ARGS[*]} -O \$OUT"

  gatk GatherVcfs "\${ARGS[@]}" -O "\$OUT"
  gatk IndexFeatureFile -I "\$OUT"
  """
}
