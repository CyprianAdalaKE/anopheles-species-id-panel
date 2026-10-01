import csv
import numpy as np

# reload full genotype matrix + missingness
with open("manifest.csv") as f:
    manifest = {row["sample_id"]: row["species"] for row in csv.DictReader(f)}
with open("results/sample_order.txt") as f:
    sample_order = [l.strip() for l in f if l.strip()]

def path_to_sample_id(p):
    parts = p.split("/")
    return f"{parts[-2]}_{parts[-1].split('.')[0]}"

sample_ids = [path_to_sample_id(s) for s in sample_order]

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
G = np.array(rows, dtype=float).T
miss = dict(zip(sample_ids, np.isnan(G).mean(axis=1)))

# load complex-only PCA scores
tail = []
with open("results/pca_scores_complex_only.tsv") as f:
    header = f.readline()
    for line in f:
        parts = line.rstrip("\n").split("\t")
        sid, sp, pc1 = parts[0], parts[1], float(parts[2])
        tail.append((sid, sp, pc1, miss.get(sid, np.nan)))

tail.sort(key=lambda x: x[2])
print("10 most negative PC1 (the 'tail') -- sample, species, PC1, missingness:")
for sid, sp, pc1, m in tail[:15]:
    print(f"  {sid:25s} {sp:12s} PC1={pc1:8.2f}  missing={m:.3f}")

print("\nMain cluster (PC1 > -10) missingness stats:")
main_miss = [m for sid, sp, pc1, m in tail if pc1 > -10]
print(f"  n={len(main_miss)}  mean_missing={np.mean(main_miss):.3f}  max={np.max(main_miss):.3f}")

print("\nTail (PC1 <= -10) missingness stats:")
tail_miss = [m for sid, sp, pc1, m in tail if pc1 <= -10]
print(f"  n={len(tail_miss)}  mean_missing={np.mean(tail_miss):.3f}  max={np.max(tail_miss):.3f}  min={np.min(tail_miss):.3f}")
