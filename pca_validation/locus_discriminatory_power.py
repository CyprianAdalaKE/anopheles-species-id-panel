#!/usr/bin/env python3
"""
Per-locus discriminatory power for the 10-locus panel.

Primary metric: Hudson's Fst (ratio-of-averages multi-SNP estimator,
Bhatia et al. 2013), computed per locus for every pairwise species
comparison -- especially the three within-gambiae-complex pairs, which
is where the panel currently struggles.

Secondary metric: per-locus leave-one-out nearest-centroid classification
accuracy for the 3-way gambiae/coluzzii/arabiensis problem, as an
easy-to-interpret complement to Fst.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = "results"
MANIFEST = "manifest.csv"
BED = "../nomadsSID.amplicons.bed"

SPECIES_ORDER = ["gambiae", "coluzzii", "arabiensis", "funestus", "stephensi", "coustani"]
COMPLEX = ["gambiae", "coluzzii", "arabiensis"]


def gt_to_dosage(gt):
    if gt in (".", "./.", ".|."):
        return np.nan
    sep = "/" if "/" in gt else "|"
    a, b = gt.split(sep)
    if a == "." or b == ".":
        return np.nan
    return int(a) + int(b)


def load_data():
    with open(MANIFEST) as f:
        manifest = {row["sample_id"]: row["species"] for row in csv.DictReader(f)}
    with open(f"{RESULTS}/sample_order.txt") as f:
        sample_order = [l.strip() for l in f if l.strip()]

    def path_to_sample_id(p):
        parts = p.split("/")
        return f"{parts[-2]}_{parts[-1].split('.')[0]}"

    sample_ids = np.array([path_to_sample_id(s) for s in sample_order])
    species = np.array([manifest[s] for s in sample_ids])

    chroms, poss = [], []
    rows = []
    with open(f"{RESULTS}/genotype_matrix.tsv") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chroms.append(parts[0])
            poss.append(int(parts[1]))
            rows.append([gt_to_dosage(g) for g in parts[4:]])
    G = np.array(rows, dtype=float).T  # samples x SNPs
    chroms = np.array(chroms)
    poss = np.array(poss)

    # drop the 6 high-missingness coluzzii samples (same QC as before)
    miss = np.isnan(G).mean(axis=1)
    keep = miss <= 0.5
    return sample_ids[keep], species[keep], G[keep], chroms, poss


def load_bed():
    loci = []
    with open(BED) as f:
        for line in f:
            c, s, e, name = line.strip().split("\t")
            loci.append((name, c, int(s), int(e)))
    return loci


def assign_snps_to_loci(chroms, poss, loci):
    locus_of_snp = np.array([None] * len(chroms), dtype=object)
    for name, c, s, e in loci:
        mask = (chroms == c) & (poss >= s) & (poss <= e)
        locus_of_snp[mask] = name
    return locus_of_snp


def hudson_fst_locus(G, species, locus_of_snp, locus_name, spA, spB):
    idxA = species == spA
    idxB = species == spB
    site_idx = np.where(locus_of_snp == locus_name)[0]
    if len(site_idx) == 0:
        return np.nan, 0
    num_sum, den_sum, n_used = 0.0, 0.0, 0
    for j in site_idx:
        gA = G[idxA, j]
        gB = G[idxB, j]
        gA = gA[~np.isnan(gA)]
        gB = gB[~np.isnan(gB)]
        nA, nB = len(gA), len(gB)
        if nA < 2 or nB < 2:
            continue
        p1 = gA.mean() / 2.0
        p2 = gB.mean() / 2.0
        if p1 == p2 == 0 or p1 == p2 == 1:
            continue
        num = (p1 - p2) ** 2 - p1 * (1 - p1) / (2 * nA - 1) - p2 * (1 - p2) / (2 * nB - 1)
        den = p1 * (1 - p2) + p2 * (1 - p1)
        if den <= 0:
            continue
        num_sum += num
        den_sum += den
        n_used += 1
    if den_sum == 0:
        return np.nan, n_used
    return num_sum / den_sum, n_used


def nearest_centroid_loo(G_locus, species, classes=COMPLEX):
    mask = np.isin(species, classes)
    X = G_locus[mask]
    y = species[mask]
    # mean-impute per-column using overall mean (recomputed per fold would be better,
    # but with high depth/low missingness within these loci this simplification is fine)
    col_mean = np.nanmean(X, axis=0)
    col_mean = np.where(np.isnan(col_mean), 0, col_mean)
    Xi = np.where(np.isnan(X), col_mean, X)

    correct = 0
    n = Xi.shape[0]
    for i in range(n):
        train_mask = np.ones(n, dtype=bool)
        train_mask[i] = False
        centroids = {}
        for c in classes:
            cm = y[train_mask] == c
            if cm.sum() == 0:
                continue
            centroids[c] = Xi[train_mask][cm].mean(axis=0)
        dists = {c: np.linalg.norm(Xi[i] - v) for c, v in centroids.items()}
        pred = min(dists, key=dists.get)
        if pred == y[i]:
            correct += 1
    return correct / n, n


def main():
    sample_ids, species, G, chroms, poss = load_data()
    loci = load_bed()
    locus_of_snp = assign_snps_to_loci(chroms, poss, loci)
    locus_names = [name for name, *_ in loci]

    print(f"Samples: {len(sample_ids)}, SNPs: {G.shape[1]}")
    for name in locus_names:
        n_snps = (locus_of_snp == name).sum()
        print(f"  {name}: {n_snps} SNPs")

    complex_pairs = [("gambiae", "coluzzii"), ("gambiae", "arabiensis"), ("coluzzii", "arabiensis")]
    all_pairs = [(a, b) for i, a in enumerate(SPECIES_ORDER) for b in SPECIES_ORDER[i+1:]]

    # ---- Fst table: loci x all pairs ----
    fst_table = {}
    for name in locus_names:
        for a, b in all_pairs:
            fst, n_used = hudson_fst_locus(G, species, locus_of_snp, name, a, b)
            fst_table[(name, a, b)] = fst

    with open(f"{RESULTS}/locus_fst_table.tsv", "w") as f:
        f.write("locus\t" + "\t".join(f"{a}-{b}" for a, b in all_pairs) + "\tmean_complex_Fst\n")
        for name in locus_names:
            vals = [fst_table[(name, a, b)] for a, b in all_pairs]
            complex_vals = [fst_table[(name, a, b)] for a, b in complex_pairs]
            complex_mean = np.nanmean(complex_vals)
            f.write(name + "\t" + "\t".join(f"{v:.4f}" if not np.isnan(v) else "NA" for v in vals)
                    + f"\t{complex_mean:.4f}\n")

    print("\n=== Mean within-gambiae-complex Fst per locus (higher = more discriminatory) ===")
    ranking = []
    for name in locus_names:
        complex_vals = [fst_table[(name, a, b)] for a, b in complex_pairs]
        ranking.append((name, np.nanmean(complex_vals)))
    ranking.sort(key=lambda x: -x[1])
    for name, val in ranking:
        print(f"  {name}: {val:.4f}")

    # ---- Nearest-centroid leave-one-out accuracy per locus (complex trio) ----
    print("\n=== Per-locus 3-way (gambiae/coluzzii/arabiensis) LOO nearest-centroid accuracy ===")
    acc_ranking = []
    for name in locus_names:
        site_idx = np.where(locus_of_snp == name)[0]
        if len(site_idx) == 0:
            continue
        acc, n = nearest_centroid_loo(G[:, site_idx], species)
        acc_ranking.append((name, acc, n))
    acc_ranking.sort(key=lambda x: -x[1])
    for name, acc, n in acc_ranking:
        print(f"  {name}: {acc*100:.1f}% (n={n})")

    with open(f"{RESULTS}/locus_accuracy_ranking.tsv", "w") as f:
        f.write("locus\tloo_accuracy_complex_trio\tn_samples\n")
        for name, acc, n in acc_ranking:
            f.write(f"{name}\t{acc:.4f}\t{n}\n")

    # heatmap of Fst: loci x complex pairs
    fig, ax = plt.subplots(figsize=(6, 5))
    mat = np.array([[fst_table[(name, a, b)] for a, b in complex_pairs] for name in locus_names])
    im = ax.imshow(mat, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(complex_pairs)))
    ax.set_xticklabels([f"{a[:3]}-{b[:3]}" for a, b in complex_pairs])
    ax.set_yticks(range(len(locus_names)))
    ax.set_yticklabels(locus_names)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            ax.text(j, i, f"{v:.2f}" if not np.isnan(v) else "NA", ha="center", va="center",
                    color="white" if (not np.isnan(v) and v < mat[~np.isnan(mat)].max()*0.6) else "black", fontsize=8)
    plt.colorbar(im, label="Hudson Fst")
    ax.set_title("Per-locus Fst: gambiae-complex pairs")
    plt.tight_layout()
    plt.savefig(f"{RESULTS}/locus_fst_heatmap.png", dpi=200)
    print(f"\nSaved: {RESULTS}/locus_fst_table.tsv, locus_accuracy_ranking.tsv, locus_fst_heatmap.png")


if __name__ == "__main__":
    main()
