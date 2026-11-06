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

rows = []
with open("results/genotype_matrix.tsv") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        gts = parts[4:]
        rows.append([gt_to_dosage(g) for g in gts])
G = np.array(rows, dtype=float).T  # samples x sites

miss = np.isnan(G).mean(axis=1)
import collections
by_sp = collections.defaultdict(list)
for sid, sp, m in zip(sample_ids, species_list, miss):
    by_sp[sp].append(m)

print("Per-species mean missingness (fraction of 1022 SNP sites uncalled):")
for sp, vals in by_sp.items():
    vals = np.array(vals)
    print(f"  {sp}: mean={vals.mean():.2f}  median={np.median(vals):.2f}  >50% missing: {(vals>0.5).sum()}/{len(vals)}")

# Flag specific coluzzii samples with high missingness (candidates for the odd PC3 cluster)
print("\ncoluzzii samples with >40% missingness (barcode, missingness):")
for sid, sp, m in zip(sample_ids, species_list, miss):
    if sp == "coluzzii" and m > 0.4:
        print(f"  {sid}: {m:.2f}")
