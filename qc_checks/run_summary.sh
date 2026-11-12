#!/usr/bin/env bash
# Per-file run/flow cell/barcode summary from read headers, one TSV per final folder.
F=/home/genomics/Desktop/species_id/anopheles-panel/species_id_final_fastqs
for sp in "$@"; do
  OUT=/home/genomics/Desktop/species_id/anopheles-panel/run_summary_$sp.tsv
  printf 'file\treads\trun\tflow_cell_id\trunid\tbarcode\tbasecall_model\n' > $OUT
  for f in $F/$sp/*; do
    case "$f" in *.tsv) continue;; esac
    awk -v f="$(basename $f)" 'NR%4==1{r=b=fc=id=m="NA"; for(i=2;i<=NF;i++){split($i,a,"="); if(a[1]=="protocol_group_id")r=a[2]; else if(a[1]=="barcode")b=a[2]; else if(a[1]=="flow_cell_id")fc=a[2]; else if(a[1]=="runid")id=a[2]; else if(a[1]=="basecall_model_version_id")m=a[2]} k=r"\t"fc"\t"id"\t"b"\t"m; c[k]++} END{for(k in c) print f"\t"c[k]"\t"k}' "$f"
  done >> $OUT
done
