"""Figures 11, 12 (genotype PCA) and 14 (per-locus within-complex Fst)."""
import os, csv, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from style import SP, COL, INK, INK2, GRID, it
R = os.path.join(os.path.dirname(HERE), "analysis_verified/pca/results")
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2})
EVR_ALL, EVR_CX = (49.1, 22.2), (26.7, 10.3)

def pca_plot(f, species, evr, out, n_snps):
    d = pd.read_csv(f"{R}/{f}", sep="\t")
    fig, ax = plt.subplots(figsize=(6.6, 5.2), constrained_layout=True)
    for s in species:
        x = d[d.species == s]
        ax.scatter(x.PC1, x.PC2, s=22, color=COL[s], edgecolor="white", lw=0.6, alpha=0.9, label=f"{it(s)} (n={len(x)})", zorder=3)
    ax.set_xlabel(f"PC1 ({evr[0]}%)"); ax.set_ylabel(f"PC2 ({evr[1]}%)")
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="best", title=f"{len(d)} samples, {n_snps} SNPs", title_fontsize=8)
    fig.savefig(f"{HERE}/{out}", dpi=300, facecolor="white"); return d

pca_plot("pca_scores_filtered.tsv", SP, EVR_ALL, "figure11_pca_all.png", 980)
pca_plot("pca_scores_complex_only.tsv", ["arabiensis", "coluzzii", "gambiae"], EVR_CX, "figure12_pca_complex.png", 893)

f = pd.read_csv(f"{R}/locus_fst_table.tsv", sep="\t").set_index("locus")
pairs = [("gambiae-coluzzii", "gambiae–coluzzii"), ("gambiae-arabiensis", "gambiae–arabiensis"), ("coluzzii-arabiensis", "coluzzii–arabiensis")]
M = f[[p for p, _ in pairs]].clip(lower=0); M = M.loc[M.mean(axis=1).sort_values(ascending=False).index]
cmap = LinearSegmentedColormap.from_list("b", ["#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
fig, ax = plt.subplots(figsize=(5.4, 5.6), constrained_layout=True)
im = ax.imshow(M.to_numpy(), cmap=cmap, vmin=0, vmax=0.15, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M.iat[i, j]; ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=8, color="white" if v > 0.08 else INK)
ax.set_xticks(range(3), [lab for _, lab in pairs], fontsize=8.5); ax.set_yticks(range(M.shape[0]), M.index, fontsize=8.5)
ax.set_xticks(np.arange(-.5, 3.5), minor=True); ax.set_yticks(np.arange(-.5, M.shape[0] + .5), minor=True)
ax.grid(which="minor", color="white", lw=2); ax.tick_params(which="both", length=0)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=ax, shrink=0.7); cb.set_label("Hudson's Fst (negative values shown as 0)"); cb.outline.set_visible(False)
fig.savefig(f"{HERE}/figure14_fst_complex.png", dpi=300, facecolor="white")
print(M.round(4))
