#!/usr/bin/env bash
# Bring all per-sample analysis inputs in line with the verified set after the 2026-11-11 QC (349 samples):
#  - excluded samples: move their k2 / bracken900 / per-read / BAM / mosdepth outputs to *_excluded_qc/ folders
#  - stephensi barcode15 -> barcode16 (file renamed 2026-11-11) in BAM / mosdepth names
#  - 9 rebuilt SLCA013 files (coluzzii 03-08, arabiensis 59/67/68): move old outputs aside and regenerate
#  - regenerate file lists (accuracy6, bracken900, reviewer samples.txt) from the verified folder
set -euo pipefail
P=/home/genomics/Desktop/species_id/anopheles-panel; G=/home/genomics/Desktop/species_id/gamcol_validation
V=/media/genomics/Sequence_Bkup/species_id_verified_fastqs; DB=/home/genomics/species_id/species_id_infectivity_db
RA=$P/analysis_verified/reviewer_analyses; A6=$G/accuracy6; B9=$RA/bracken900
EXQ=/media/genomics/Sequence_Bkup/species_id_excluded_background_2026-11-11/EXCLUDED.tsv
source /home/genomics/miniforge3/etc/profile.d/conda.sh; conda activate anopheles-panel

aside() { # $1 = dir, $2 = glob; move matches into $1_excluded_qc (or $1/excluded_qc for nested dirs)
  local d=$1 g=$2 dest; dest=${3:-$1/excluded_qc}; mkdir -p "$dest"
  for f in $d/$g; do [ -e "$f" ] && mv "$f" "$dest/"; done; return 0; }

STEPS13_DONE=$P/run_checks/.rerun_steps1-3.done
if [ ! -e $STEPS13_DONE ]; then
# 1. excluded samples (stems as used by each output type)
tail -n +2 $EXQ | while IFS=$'\t' read sp f reason; do
  stem=${f%.fastq.gz}; s=${stem#combined.}; s1=$stem; [ "$sp" = coluzzii ] || [ "$sp" = coustani ] && s1=$stem.fastq
  aside $A6/k2 "${sp}__${stem}[._]*"; aside $A6/k2 "${sp}__${s1}[._]*"
  aside $B9/k2 "${sp}__${stem}[._]*"; aside $B9/k2 "${sp}__${s1}[._]*"
  aside $RA/per_read "${sp}__${stem}.*"
  aside $P/aligned_verified/$sp "$s.sorted.bam*" $P/aligned_verified_excluded_qc/$sp
  aside $P/coverage_results_verified "${sp}_$s.*" $P/coverage_results_verified_excluded_qc
done

# 2. stephensi barcode15 -> barcode16 names (BAM, mosdepth)
for f in $P/aligned_verified/stephensi/barcode15.sorted.bam*; do [ -e "$f" ] && mv "$f" "${f/barcode15/barcode16}"; done
for f in $P/coverage_results_verified/stephensi_barcode15.*; do [ -e "$f" ] && mv "$f" "${f/barcode15/barcode16}"; done

# 3. rebuilt SLCA013 files: move old outputs aside
for x in coluzzii:03 coluzzii:04 coluzzii:05 coluzzii:06 coluzzii:07 coluzzii:08 arabiensis:59 arabiensis:67 arabiensis:68; do
  sp=${x%%:*}; n=${x##*:}; stem=combined.barcode$n; s1=$stem; [ $sp = coluzzii ] && s1=$stem.fastq
  aside $A6/k2 "${sp}__${s1}[._]*" $A6/k2_replaced_BatchCOL; aside $B9/k2 "${sp}__${s1}[._]*" $B9/k2_replaced_BatchCOL
  aside $RA/per_read "${sp}__${stem}.*" $RA/per_read_replaced_BatchCOL
  aside $P/aligned_verified/$sp "barcode$n.sorted.bam*" $P/aligned_verified_replaced_BatchCOL/$sp
  aside $P/coverage_results_verified "${sp}_barcode$n.*" $P/coverage_results_verified_replaced_BatchCOL
done

touch $STEPS13_DONE; fi   # steps 1-3 move files and must run only once

# 4. file lists from the verified folder (coluzzii/coustani keep the historical ".fastq" stem in accuracy6/bracken900)
for L in $A6/files.tsv $B9/files.tsv; do
  for f in $V/*/*.fastq.gz; do sp=$(basename $(dirname $f)); s=$(basename $f .fastq.gz)
    [ $sp = coluzzii ] || [ $sp = coustani ] && s=$s.fastq; printf '%s\t%s\n' $sp $s; done | sort -k1,1 -k2,2V > $L
done
for f in $V/*/*.fastq.gz; do echo "$(basename $(dirname $f)) $(basename $f)"; done | sort -k1,1 -k2,2V > $RA/samples.txt
echo "file lists: $(wc -l < $A6/files.tsv) samples"

# 5. regenerate missing outputs: Kraken2 + Bracken (100 and 900), per-read subsample, alignment + mosdepth
while IFS=$'\t' read sp s; do
  s=${s%.fastq}; f=$V/$sp/$s.fastq.gz; o=$A6/k2/${sp}__$s; o9=$B9/k2/${sp}__$s   # k2 outputs use the stem without .fastq
  if [ ! -s $o.kreport ]; then
    kraken2 --db $DB --threads 16 --gzip-compressed --report $o.kreport --output /dev/null $f 2> $o.k2.log
    python3 -I /home/genomics/miniforge3/bin/est_abundance.py -i $o.kreport -o $o.bracken -k $DB/database100mers.kmer_distrib -l S -t 10 > $o.bracken.log 2>&1 || true
  fi
  if [ ! -s $o9.bracken ]; then
    ln -sf $o.kreport $o9.kreport   # bracken900 reuses the main Kraken2 report (as before)
    python3 -I /home/genomics/miniforge3/bin/est_abundance.py -i $o9.kreport -o $o9.bracken -k $DB/database900mers.kmer_distrib -l S -t 10 > $o9.bracken.log 2>&1 || true
  fi
done < $A6/files.tsv
while read sp f; do bash $RA/run_sample.sh $sp $f; done < $RA/samples.txt
bash $P/align_coverage_verified.sh
echo "=== rerun_after_qc done"
