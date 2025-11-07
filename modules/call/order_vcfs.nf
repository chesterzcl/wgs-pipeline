process ORDER_VCFS {
  tag "order_vcfs"
  publishDir "${params.outdir}/vcf", mode: 'copy', pattern: 'ordered.snps.list'

  input:
    path order_file
    path pairs_file   

  output:
    path "ordered.snps.list"

  script:
  """
  set -euo pipefail

  awk 'NF{gsub(/\\r/,""); print}' "${order_file}" > order.tmp
  awk 'NF{gsub(/\\r/,""); print}' "${pairs_file}" > pairs.tmp
  
  wc -l order.tmp | awk '{print "[INFO] intervals in order.tmp : " \$1}'
  wc -l pairs.tmp | awk '{print "[INFO] pairs in pairs.tmp     : " \$1}'

  bad_lines=\$(awk 'BEGIN{FS="\\t"} NF!=2{c++} END{print c+0}' pairs.tmp)
  if [[ "\${bad_lines}" -gt 0 ]]; then
    echo "[WARN] \${bad_lines} line(s) in pairs.tmp do not have exactly 2 TAB-separated fields." >&2
  fi

  awk 'BEGIN{FS=OFS="\\t"}
       NR==FNR { a[\$1]=\$2; next }   # pairs.tmp -> a[interval]=path
       { k=\$1; if (k in a) print a[k] }
      ' pairs.tmp order.tmp > ordered.snps.list

  wc -l ordered.snps.list | awk '{print "[INFO] written lines in ordered.snps.list : " \$1}'
  echo "[INFO] head of ordered.snps.list:"
  head -n 10 ordered.snps.list || true

  """
}
