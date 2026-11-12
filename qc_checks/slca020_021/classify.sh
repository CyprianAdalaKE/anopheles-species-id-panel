#!/usr/bin/env bash
# One Kraken2 pass over all subsampled reads (same DB as the accuracy analysis); per-read output + combined report.
set -eu
DB=/home/genomics/species_id/species_id_infectivity_db
kraken2 --db $DB --threads 18 --report subsample.kreport --output subsample.k2out subsample.fastq 2> kraken2.log
echo done
