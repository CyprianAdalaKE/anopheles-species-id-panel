#!/bin/bash
# Joint diploid SNP calling across all samples at the 10 nuclear target loci.
# Commands recovered verbatim from the ##bcftoolsCommand header of the original
# vcf/raw_calls.vcf.gz (bcftools 1.24). Run build_manifest.py first; bam_list.txt
# is the bam column of manifest.csv.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p vcf results logs

tail -n +2 manifest.csv | cut -d, -f4 > bam_list.txt

bcftools mpileup -R ../nomadsSID.amplicons.bed -f ../VectorBase-67_AgambiaePEST_Genome.fasta \
    -B -I -a AD,DP -q 20 -Q 15 -d 60 -Ou -b bam_list.txt 2> logs/mpileup_full.log \
  | bcftools call -mv --ploidy 2 -Oz -o vcf/raw_calls.vcf.gz
bcftools index -t vcf/raw_calls.vcf.gz

# Next: filter_and_matrix.sh, then pca.py / pca2.py / locus_discriminatory_power.py
