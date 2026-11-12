#!/usr/bin/env bash
# Kraken2 + Bracken on the verified coluzzii and coustani files (n=30 each, species_id_verified_fastqs),
# same DB and Bracken call as run_k2_verified.sh (accuracy4); the four-species outputs are reused from accuracy4.
set -eu
W=/home/genomics/Desktop/species_id/gamcol_validation/accuracy6
DB=/home/genomics/species_id/species_id_infectivity_db
V=/media/genomics/Sequence_Bkup/species_id_verified_fastqs
for sp in coluzzii coustani; do
  for f in $V/$sp/*.fastq.gz; do
    b=$(basename $f .fastq.gz); o=$W/k2/${sp}__$b
    [ -s "$o.bracken" ] && continue
    kraken2 --db $DB --threads 16 --gzip-compressed --report "$o.kreport" --output /dev/null "$f" 2>"$o.k2.log"
    python3 -I /home/genomics/miniforge3/bin/est_abundance.py -i "$o.kreport" -o "$o.bracken" \
      -k $DB/database100mers.kmer_distrib -l S -t 10 > "$o.bracken.log" 2>&1 || true
  done
done
echo done
