#!/usr/bin/env python3
"""
Genotype-based PCA validation for the Anopheles amplicon panel.

Independent of Kraken2/Bracken: builds a diploid dosage genotype matrix
(0/1/2) from joint bcftools SNP calls at the 10 nuclear target loci,
Patterson-scales it (as in EIGENSOFT/smartpca), and runs PCA via SVD.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = "results"
MANIFEST = "manifest.csv"

SPECIES_COLORS = {
    "gambiae":    "#1b9e77",
    "coluzzii":   "#d95f02",
    "arabiensis": "#7570b3",
    "funestus":   "#e7298a",
    "stephensi":  "#66a61e",
    "coustani":   "#e6ab02",
}


def gt_to_dosage(gt):
    # bcftools GT strings: "0/0", "0/1", "1/1", "./." (missing)
    if gt in (".", "./.", ".|."):
        return np.nan
    sep = "/" if "/" in gt else "|"
    a, b = gt.split(sep)
    if a == "." or b == ".":
        return np.nan
    return int(a) + int(b)


def load_genotype_matrix(path, n_samples):
    sites = []
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = parts[0:4]
            gts = parts[4:4 + n_samples]
            dosages = [gt_to_dosage(g) for g in gts]
            sites.append(f"{chrom}:{pos}_{ref}>{alt}")
            rows.append(dosages)
    G = np.array(rows, dtype=float)  # sites x samples
    return G.T, sites  # samples x sites


def patterson_pca(G, n_components=5):
    """G: samples x SNPs, values in {0,1,2,NaN}. Returns scores, explained_var_ratio, kept_snp_mask."""
    n, m = G.shape
    p = np.nanmean(G, axis=0) / 2.0  # allele frequency per SNP

    # drop monomorphic / all-missing sites
    keep = (p > 0) & (p < 1) & ~np.isnan(p)
    G = G[:, keep]
    p = p[keep]

    # mean-impute missing genotypes to 2p
    impute = np.tile(2 * p, (n, 1))
    mask = np.isnan(G)
    G = np.where(mask, impute, G)

    # Patterson centering + scaling
    denom = np.sqrt(p * (1 - p))
    M = (G - 2 * p) / denom

    # SVD (samples x SNPs), PCs are U*S
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    var = S ** 2
    evr = var / var.sum()
    scores = U * S
    return scores[:, :n_components], evr[:n_components], keep.sum()


def main():
    with open(MANIFEST) as f:
        manifest = {row["sample_id"]: row["species"] for row in csv.DictReader(f)}

    with open(f"{RESULTS}/sample_order.txt") as f:
        sample_order = [l.strip() for l in f if l.strip()]

    # bcftools query -l returns the BAM file paths used as sample names by default;
    # map back to species via the basename pattern <species>/<barcode>.sorted.bam
    def path_to_sample_id(p):
        parts = p.rstrip("\n").split("/")
        species = parts[-2]
        barcode = parts[-1].split(".")[0]
        return f"{species}_{barcode}"

    sample_ids = [path_to_sample_id(s) for s in sample_order]
    species_list = [manifest[s] for s in sample_ids]

    G, sites = load_genotype_matrix(f"{RESULTS}/genotype_matrix.tsv", len(sample_ids))
    print(f"Loaded genotype matrix: {G.shape[0]} samples x {G.shape[1]} sites")

    scores, evr, n_kept = patterson_pca(G, n_components=5)
    print(f"SNPs retained after monomorphic filter: {n_kept}")
    print("Explained variance ratio (PC1-5):", np.round(evr, 4))

    # Save scores table
    with open(f"{RESULTS}/pca_scores.tsv", "w") as f:
        f.write("sample_id\tspecies\t" + "\t".join(f"PC{i+1}" for i in range(scores.shape[1])) + "\n")
        for sid, sp, row in zip(sample_ids, species_list, scores):
            f.write(sid + "\t" + sp + "\t" + "\t".join(f"{v:.4f}" for v in row) + "\n")

    # Plot PC1 vs PC2
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for sp, color in SPECIES_COLORS.items():
        idx = [i for i, s in enumerate(species_list) if s == sp]
        if not idx:
            continue
        axes[0].scatter(scores[idx, 0], scores[idx, 1], s=18, alpha=0.75, label=sp, color=color)
        axes[1].scatter(scores[idx, 0], scores[idx, 2], s=18, alpha=0.75, label=sp, color=color)

    axes[0].set_xlabel(f"PC1 ({evr[0]*100:.1f}%)")
    axes[0].set_ylabel(f"PC2 ({evr[1]*100:.1f}%)")
    axes[0].set_title("Genotype PCA: PC1 vs PC2")
    axes[0].legend(fontsize=8)

    axes[1].set_xlabel(f"PC1 ({evr[0]*100:.1f}%)")
    axes[1].set_ylabel(f"PC3 ({evr[2]*100:.1f}%)")
    axes[1].set_title("Genotype PCA: PC1 vs PC3")
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(f"{RESULTS}/pca_species.png", dpi=200)
    print(f"Saved {RESULTS}/pca_species.png and {RESULTS}/pca_scores.tsv")


if __name__ == "__main__":
    main()
