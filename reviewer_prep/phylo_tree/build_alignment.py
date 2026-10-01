#!/usr/bin/env python3
"""
Build a SNP-only pseudo-sequence alignment (one row per sample, one column
per biallelic SNP) from the panel's genotype matrix, for ML tree building
with IQ-TREE using ascertainment-bias correction (standard practice for
SNP-only alignments, since invariant sites were already excluded during
variant calling/filtering).
"""
import csv
import numpy as np

with open("../../pca_validation/manifest.csv") as f:
    manifest = {row["sample_id"]: row["species"] for row in csv.DictReader(f)}
with open("../../pca_validation/results/sample_order.txt") as f:
    sample_order = [l.strip() for l in f if l.strip()]


def path_to_sample_id(p):
    parts = p.split("/")
    return f"{parts[-2]}_{parts[-1].split('.')[0]}"


sample_ids = [path_to_sample_id(s) for s in sample_order]
species = [manifest[s] for s in sample_ids]

IUPAC = {
    frozenset("AA"): "A", frozenset("CC"): "C", frozenset("GG"): "G", frozenset("TT"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GC"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M",
}


def encode(gt, ref, alt):
    if gt in (".", "./.", ".|."):
        return "N"
    sep = "/" if "/" in gt else "|"
    a, b = gt.split(sep)
    if a == "." or b == ".":
        return "N"
    bases = (ref if a == "0" else alt) + (ref if b == "0" else alt)
    return IUPAC.get(frozenset(bases), "N")


sites = []
matrix = []  # will be n_sites x n_samples, transpose later
with open("../../pca_validation/results/genotype_matrix.tsv") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt = parts[0:4]
        if len(ref) != 1 or len(alt) != 1:
            continue  # keep SNPs only
        gts = parts[4:]
        row = [encode(g, ref, alt) for g in gts]
        matrix.append(row)
        sites.append(f"{chrom}:{pos}")

matrix = np.array(matrix)  # sites x samples
n_sites, n_samples = matrix.shape
print(f"{n_sites} SNPs x {n_samples} samples")

# per-sample missingness QC (drop the known coluzzii low-coverage batch)
miss_per_sample = (matrix == "N").mean(axis=0)
keep = miss_per_sample <= 0.5
print(f"Dropping {(~keep).sum()} samples with >50% missing genotypes")

kept_matrix = matrix[:, keep]

# drop sites that are invariant (or all-N) within the retained sample set --
# required for IQ-TREE's +ASC ascertainment-bias correction
def is_variable(row):
    bases = set(row) - {"N"}
    # a het IUPAC code alone still counts as variation from a homozygous site;
    # require at least 2 distinct non-N symbols to call a site variable
    return len(bases) >= 2

site_variable = np.array([is_variable(kept_matrix[i]) for i in range(n_sites)])
print(f"Dropping {(~site_variable).sum()} sites invariant within the retained {keep.sum()} samples")
kept_matrix = kept_matrix[site_variable]
n_sites_final = kept_matrix.shape[0]

with open("snp_alignment.fasta", "w") as f:
    kept_i = 0
    for i in range(n_samples):
        if not keep[i]:
            continue
        seq = "".join(kept_matrix[:, kept_i])
        kept_i += 1
        label = f"{sample_ids[i]}_{species[i]}".replace(" ", "_")
        f.write(f">{label}\n{seq}\n")

with open("sample_species_map.tsv", "w") as f:
    for i in range(n_samples):
        if not keep[i]:
            continue
        label = f"{sample_ids[i]}_{species[i]}".replace(" ", "_")
        f.write(f"{label}\t{species[i]}\n")

print(f"Wrote snp_alignment.fasta ({keep.sum()} samples x {n_sites_final} sites)")
