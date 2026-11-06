#!/usr/bin/env bash
# Joint genotyping of verified-sample BAMs at the 10 nuclear loci (same commands as pca_validation, recorded in its VCF header).
set -euo pipefail
source /home/genomics/miniforge3/etc/profile.d/conda.sh; conda activate anopheles-panel
cd /home/genomics/Desktop/species_id/anopheles-panel/analysis_verified/pca
ls /home/genomics/Desktop/species_id/anopheles-panel/aligned_verified/*/*.sorted.bam > bam_list.txt
python3 -I - <<'PY'
import csv
rows=[]
for b in open('bam_list.txt').read().split():
    sp=b.split('/')[-2]; bc=b.split('/')[-1].split('.')[0]
    rows.append(dict(sample_id=f"{sp}_{bc}",species=sp,barcode=bc,bam=b))
w=csv.DictWriter(open('manifest.csv','w',newline=''),fieldnames=['sample_id','species','barcode','bam']); w.writeheader(); w.writerows(rows)
print(len(rows),'BAMs')
PY
bcftools mpileup -R /home/genomics/Desktop/species_id/anopheles-panel/nomadsSID.amplicons.bed -f /home/genomics/Desktop/species_id/anopheles-panel/VectorBase-67_AgambiaePEST_Genome.fasta -B -I -a AD,DP -q 20 -Q 15 -d 60 -Ou -b bam_list.txt 2>logs/mpileup.log \
  | bcftools call -mv --ploidy 2 -Oz -o vcf/raw_calls.vcf.gz 2>logs/call.log
bcftools index -t vcf/raw_calls.vcf.gz
bcftools view -m2 -M2 -v snps vcf/raw_calls.vcf.gz | bcftools +fill-tags -- -t AF \
  | bcftools view -e 'F_MISSING > 0.15 || (AF < 0.05 && AF > 0) || MAF < 0.05' -Oz -o vcf/filtered.vcf.gz
bcftools index -t vcf/filtered.vcf.gz
bcftools query -l vcf/filtered.vcf.gz > results/sample_order.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' vcf/filtered.vcf.gz > results/genotype_matrix.tsv
echo "sites raw: $(zcat vcf/raw_calls.vcf.gz | grep -vc '^#')  filtered: $(wc -l < results/genotype_matrix.tsv)  samples: $(wc -l < results/sample_order.txt)"
