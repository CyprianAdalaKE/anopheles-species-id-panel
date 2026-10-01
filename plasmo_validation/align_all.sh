#!/bin/bash
# Align every infectivity-experiment sample to (a) the AgamP4 host genome and
# (b) the combined 4-species Plasmodium reference, in parallel (xargs -P).
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate anopheles-panel
cd "$(dirname "$0")"

FASTQ_DIR=../infectivity
THREADS_PER_JOB=2
PARALLEL_JOBS=8

align_one() {
    bc="$1"
    fq="$FASTQ_DIR/combined.${bc}.fastq"
    if [ ! -s "$fq" ]; then
        echo "SKIP $bc (missing/empty fastq)" >> logs/align.log
        return
    fi
    if [ ! -s "bams/${bc}.host.sorted.bam" ]; then
        minimap2 -ax map-ont -t $THREADS_PER_JOB agamp4.mmi "$fq" 2>>logs/align.log \
          | samtools sort -@ $THREADS_PER_JOB -o "bams/${bc}.host.sorted.bam" - 2>>logs/align.log
        samtools index "bams/${bc}.host.sorted.bam"
    fi
    if [ ! -s "bams/${bc}.plasmo.sorted.bam" ]; then
        minimap2 -ax map-ont -t $THREADS_PER_JOB plasmo_combined.mmi "$fq" 2>>logs/align.log \
          | samtools sort -@ $THREADS_PER_JOB -o "bams/${bc}.plasmo.sorted.bam" - 2>>logs/align.log
        samtools index "bams/${bc}.plasmo.sorted.bam"
    fi
    echo "DONE $bc" >> logs/align.log
}
export -f align_one
export FASTQ_DIR THREADS_PER_JOB

cat barcodes.txt | xargs -P $PARALLEL_JOBS -I{} bash -c 'align_one "$@"' _ {}
echo "ALL DONE" >> logs/align.log
