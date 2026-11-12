#!/usr/bin/env bash
# Per sample: subsample 5,000 reads (seqtk, seed 100); align to 60 species-specific amplicon references (all hits kept);
# Kraken2 per-read classification of the same reads (same DB as the main analysis).
set -eu
sp=$1; f=$2; R=/home/genomics/Desktop/species_id/anopheles-panel/analysis_verified/reviewer_analyses
V=/media/genomics/Sequence_Bkup/species_id_verified_fastqs; DB=/home/genomics/species_id/species_id_infectivity_db
o=$R/per_read/${sp}__${f%.fastq.gz}; [ -s $o.paf.gz ] && [ -s $o.k2.gz ] && exit 0
seqtk sample -s100 $V/$sp/$f 5000 > $o.fq
minimap2 -x map-ont -c -N 100 -p 0 --secondary=yes -t 2 $R/amplicons_ref.fa $o.fq 2>/dev/null | cut -f1-12 | gzip > $o.paf.gz
kraken2 --db $DB --memory-mapping --threads 2 --output /dev/stdout $o.fq 2>/dev/null | cut -f2,3 | gzip > $o.k2.gz
rm $o.fq
