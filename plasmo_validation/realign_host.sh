#!/bin/bash
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate anopheles-panel
cd "$(dirname "$0")"

FASTQ_DIR=../infectivity
PARALLEL_JOBS=8

align_one() {
    bc="$1"
    fq="$FASTQ_DIR/combined.${bc}.fastq"
    [ -s "$fq" ] || { echo "SKIP $bc"; return; }
    minimap2 -ax map-ont -t 2 agamp4_plus_18s.mmi "$fq" 2>>logs/realign_host.log \
      | samtools sort -@ 2 -o "bams/${bc}.host18s.sorted.bam" - 2>>logs/realign_host.log
    samtools index "bams/${bc}.host18s.sorted.bam"
    echo "DONE $bc" >> logs/realign_host.log
}
export -f align_one
export FASTQ_DIR

cat barcodes.txt | xargs -P $PARALLEL_JOBS -I{} bash -c 'align_one "$@"' _ {}
echo "ALL DONE" >> logs/realign_host.log
