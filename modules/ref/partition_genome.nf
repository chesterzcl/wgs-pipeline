process PARTITION_GENOME {
  tag "W=${params.len_interval_mb}Mb"
  publishDir "${params.outdir}/intervals", mode: 'copy', pattern: "intervals.fixed.${params.len_interval_mb}Mb.txt"

  input:
    path fai

  output:
    path "intervals.fixed.${params.len_interval_mb}Mb.txt", emit: list

  """
  set -euo pipefail

  if [[ ! -s "${fai}" ]]; then
    echo "ERROR: FAI not found or empty: ${fai}" >&2
    exit 1
  fi

  W=\$(( ${params.len_interval_mb} * 1000000 ))
  if [[ "\$W" -lt 1 ]]; then
    echo "ERROR: params.len_interval_mb must be >= 1" >&2
    exit 1
  fi

  awk -v W="\$W" '{
    chr=\$1; len=\$2;
    for (start=1; start<=len; start+=W) {
      end = start + W - 1;
      if (end > len) end = len;
      printf "%s:%d-%d\\n", chr, start, end;
    }
  }' "${fai}" > "intervals.fixed.${params.len_interval_mb}Mb.txt"
  """
}

