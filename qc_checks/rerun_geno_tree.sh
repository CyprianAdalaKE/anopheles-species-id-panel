#!/usr/bin/env bash
# Genotyping -> PCA/LDA/Fst/diagnostic SNPs -> SNP alignment -> IQ-TREE (-nm 5000, same seed/model/threads as before), then read stats.
set -uo pipefail
P=/home/genomics/Desktop/species_id/anopheles-panel
source /home/genomics/miniforge3/etc/profile.d/conda.sh; conda activate anopheles-panel
cd $P/analysis_verified/pca && bash run_genotyping.sh 2>&1 | tail -1
for s in check_missingness.py pca.py pca2.py lda_complex.py find_diagnostic_snps.py locus_discriminatory_power.py; do echo "== $s"; python3 $s 2>&1 | grep -v -i warn | tail -25; done
cd $P/analysis_verified/phylo_tree && python3 build_alignment.py | tail -1
seqkit stats -a -T -j 8 /media/genomics/Sequence_Bkup/species_id_verified_fastqs/*/*.fastq.gz > $P/verified_stats/read_stats.tsv && echo "read_stats rows: $(($(wc -l < $P/verified_stats/read_stats.tsv)-1))"
echo "=== geno done"
/home/genomics/miniforge3/envs/phylo_env/bin/iqtree -s snp_alignment.fasta -st DNA -m GTR+G -B 1000 -nm 5000 -seed 669192 -T 4 -pre snp_tree_nm5000 -redo > iqtree_nm5000_run.log 2>&1
echo "=== iqtree done"
