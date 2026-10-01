#!/bin/bash
# Filter joint calls to biallelic SNPs with reasonable MAF/missingness/depth,
# then dump a genotype matrix (CHROM POS REF ALT + one GT column per sample).
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate anopheles-panel
cd "$(dirname "$0")"

IN=vcf/raw_calls.vcf.gz
FILT=vcf/filtered.vcf.gz

# Biallelic SNPs only; minimum call rate 85%; MAF > 0.05; minimum depth 8 per genotype
bcftools view -m2 -M2 -v snps "$IN" \
  | bcftools +fill-tags -- -t AF \
  | bcftools view -e 'F_MISSING > 0.15 || (AF < 0.05 && AF > 0) || MAF < 0.05' \
  -Oz -o "$FILT"
bcftools index -t "$FILT"

echo "Sites before filtering:"
zcat "$IN" | grep -vc "^#"
echo "Sites after filtering:"
zcat "$FILT" | grep -vc "^#"

# Sample order
bcftools query -l "$FILT" > results/sample_order.txt

# Genotype matrix: one row per site, columns = CHROM POS REF ALT then GT per sample
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$FILT" > results/genotype_matrix.tsv

echo "Wrote results/genotype_matrix.tsv and results/sample_order.txt"
