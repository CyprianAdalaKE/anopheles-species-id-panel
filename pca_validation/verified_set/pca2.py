#!/usr/bin/env python3
"""
Round 2: (a) full-panel PCA excluding high-missingness samples (the flagged
low-coverage coluzzii batch), (b) gambiae-complex-only PCA to isolate
within-complex genotype signal that PC1/PC2 of the full run is swamped by
stephensi/funestus divergence.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = "results"
MANIFEST = "manifest.csv"
MISSINGNESS_THRESH = 0.5

SPECIES_COLORS = {
    "gambiae":    "#1b9e77",
    "coluzzii":   "#d95f02",
    "arabiensis": "#7570b3",
    "funestus":   "#e7298a",
    "stephensi":  "#66a61e",
    "coustani":   "#e6ab02",
}


def gt_to_dosage(gt):
    if gt in (".", "./.", ".|."):
        return np.nan
    sep = "/" if "/" in gt else "|"
    a, b = gt.split(sep)
    if a == "." or b == ".":
        return np.nan
    return int(a) + int(b)


def load_all():
    with open(MANIFEST) as f:
        manifest = {row["sample_id"]: row["species"] for row in csv.DictReader(f)}
    with open(f"{RESULTS}/sample_order.txt") as f:
        sample_order = [l.strip() for l in f if l.strip()]

    def path_to_sample_id(p):
        parts = p.split("/")
        return f"{parts[-2]}_{parts[-1].split('.')[0]}"

    sample_ids = [path_to_sample_id(s) for s in sample_order]
    species_list = [manifest[s] for s in sample_ids]

    rows = []
    with open(f"{RESULTS}/genotype_matrix.tsv") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            gts = parts[4:]
            rows.append([gt_to_dosage(g) for g in gts])
    G = np.array(rows, dtype=float).T  # samples x sites
    return np.array(sample_ids), np.array(species_list), G


def patterson_pca(G, n_components=5):
    n, m = G.shape
    p = np.nanmean(G, axis=0) / 2.0
    keep = (p > 0) & (p < 1) & ~np.isnan(p)
    G = G[:, keep]
    p = p[keep]
    impute = np.tile(2 * p, (n, 1))
    mask = np.isnan(G)
    G = np.where(mask, impute, G)
    denom = np.sqrt(p * (1 - p))
    M = (G - 2 * p) / denom
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    var = S ** 2
    evr = var / var.sum()
    scores = U * S
    return scores[:, :n_components], evr[:n_components], keep.sum()


def plot_pca(scores, evr, species_list, title, outpath, pcs=(0, 1)):
    fig, ax = plt.subplots(figsize=(8, 7))
    for sp, color in SPECIES_COLORS.items():
        idx = [i for i, s in enumerate(species_list) if s == sp]
        if not idx:
            continue
        ax.scatter(scores[idx, pcs[0]], scores[idx, pcs[1]], s=22, alpha=0.75, label=sp, color=color)
    ax.set_xlabel(f"PC{pcs[0]+1} ({evr[pcs[0]]*100:.1f}%)")
    ax.set_ylabel(f"PC{pcs[1]+1} ({evr[pcs[1]]*100:.1f}%)")
    ax.set_title(title)
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def main():
    sample_ids, species_list, G = load_all()
    miss = np.isnan(G).mean(axis=1)
    keep_mask = miss <= MISSINGNESS_THRESH
    dropped = sample_ids[~keep_mask]
    print(f"Dropping {len(dropped)} samples with >{MISSINGNESS_THRESH*100:.0f}% missing genotypes:")
    for d in dropped:
        print(f"  {d}")

    sid2 = sample_ids[keep_mask]
    sp2 = species_list[keep_mask]
    G2 = G[keep_mask]

    # (a) Full panel, high-missingness samples removed
    scores, evr, n_kept = patterson_pca(G2, n_components=5)
    print(f"\n[Full panel, filtered] {G2.shape[0]} samples x {n_kept} SNPs")
    print("Explained variance ratio (PC1-5):", np.round(evr, 4))
    plot_pca(scores, evr, sp2, "Genotype PCA (high-missingness samples removed): PC1 vs PC2",
              f"{RESULTS}/pca_species_filtered.png", pcs=(0, 1))

    with open(f"{RESULTS}/pca_scores_filtered.tsv", "w") as f:
        f.write("sample_id\tspecies\t" + "\t".join(f"PC{i+1}" for i in range(scores.shape[1])) + "\n")
        for sid, sp, row in zip(sid2, sp2, scores):
            f.write(sid + "\t" + sp + "\t" + "\t".join(f"{v:.4f}" for v in row) + "\n")

    # (b) gambiae-complex only: gambiae, coluzzii, arabiensis
    complex_mask = np.isin(sp2, ["gambiae", "coluzzii", "arabiensis"])
    Gc = G2[complex_mask]
    spc = sp2[complex_mask]
    sidc = sid2[complex_mask]

    scores_c, evr_c, n_kept_c = patterson_pca(Gc, n_components=5)
    print(f"\n[gambiae-complex only] {Gc.shape[0]} samples x {n_kept_c} SNPs (re-filtered for polymorphism within complex)")
    print("Explained variance ratio (PC1-5):", np.round(evr_c, 4))
    plot_pca(scores_c, evr_c, spc, "Genotype PCA — gambiae complex only (gambiae/coluzzii/arabiensis): PC1 vs PC2",
              f"{RESULTS}/pca_complex_only.png", pcs=(0, 1))
    plot_pca(scores_c, evr_c, spc, "Genotype PCA — gambiae complex only: PC2 vs PC3",
              f"{RESULTS}/pca_complex_only_pc23.png", pcs=(1, 2))

    with open(f"{RESULTS}/pca_scores_complex_only.tsv", "w") as f:
        f.write("sample_id\tspecies\t" + "\t".join(f"PC{i+1}" for i in range(scores_c.shape[1])) + "\n")
        for sid, sp, row in zip(sidc, spc, scores_c):
            f.write(sid + "\t" + sp + "\t" + "\t".join(f"{v:.4f}" for v in row) + "\n")

    print(f"\nSaved: pca_species_filtered.png, pca_complex_only.png, pca_complex_only_pc23.png")


if __name__ == "__main__":
    main()
