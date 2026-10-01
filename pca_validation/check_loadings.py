import csv
import numpy as np

with open("manifest.csv") as f:
    manifest = {row["sample_id"]: row["species"] for row in csv.DictReader(f)}
with open("results/sample_order.txt") as f:
    sample_order = [l.strip() for l in f if l.strip()]

def path_to_sample_id(p):
    parts = p.split("/")
    return f"{parts[-2]}_{parts[-1].split('.')[0]}"

sample_ids = [path_to_sample_id(s) for s in sample_order]
species_list = [manifest[s] for s in sample_ids]

def gt_to_dosage(gt):
    if gt in (".", "./.", ".|."):
        return np.nan
    sep = "/" if "/" in gt else "|"
    a, b = gt.split(sep)
    if a == "." or b == ".":
        return np.nan
    return int(a) + int(b)

sites = []
rows = []
with open("results/genotype_matrix.tsv") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt = parts[0:4]
        sites.append(f"{chrom}:{pos}_{ref}>{alt}")
        rows.append([gt_to_dosage(g) for g in parts[4:]])
G_all = np.array(rows, dtype=float).T
sites = np.array(sites)

sample_ids = np.array(sample_ids)
species_list = np.array(species_list)
complex_mask = np.isin(species_list, ["gambiae", "coluzzii", "arabiensis"])
G = G_all[complex_mask]
sp = species_list[complex_mask]
sid = sample_ids[complex_mask]

n, m = G.shape
p = np.nanmean(G, axis=0) / 2.0
keep = (p > 0) & (p < 1) & ~np.isnan(p)
Gk = G[:, keep]
pk = p[keep]
sitesk = sites[keep]
impute = np.tile(2*pk, (n,1))
Gk = np.where(np.isnan(Gk), impute, Gk)
denom = np.sqrt(pk*(1-pk))
M = (Gk - 2*pk) / denom
U, S, Vt = np.linalg.svd(M, full_matrices=False)

pc1_loadings = Vt[0]
order = np.argsort(-np.abs(pc1_loadings))
print("Top 10 SNPs loading onto PC1 (complex-only run):")
for i in order[:10]:
    print(f"  {sitesk[i]:35s} loading={pc1_loadings[i]:+.3f}  MAF={pk[i]:.3f}")

# how much does the single top SNP explain of PC1 variance?
top = order[0]
print(f"\nTop SNP squared loading: {pc1_loadings[top]**2:.4f} (out of total 1.0 across {len(pc1_loadings)} SNPs)")
print(f"Sum of top-5 squared loadings: {sum(pc1_loadings[order[:5]]**2):.4f}")
