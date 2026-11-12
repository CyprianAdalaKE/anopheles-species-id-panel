#!/usr/bin/env bash
# Copy the 374 read files into one upload folder under their submission names (sample_id.fastq.gz) and verify MD5s.
# Usage: bash stage_upload_files.sh <upload_dir>   (needs ~ the size of the verified folder; source files are not modified)
set -euo pipefail
OUT=${1:?upload dir}; mkdir -p "$OUT"; cd "$(dirname "$0")"
tail -n +2 files_manifest.tsv | while IFS=$'\t' read sid sp qc reason run bc reads src up; do cp -n "$src" "$OUT/$up"; done
awk 'NR==FNR{m[$2]=$1; next} FNR>1{print m[$8]"  "$9}' md5_source_files.txt files_manifest.tsv > "$OUT/md5sums.txt"
(cd "$OUT" && md5sum -c md5sums.txt | grep -v ': OK$' || true)
echo "staged $(ls "$OUT"/*.fastq.gz | wc -l) files in $OUT"
