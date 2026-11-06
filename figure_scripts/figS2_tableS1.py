"""Supplementary Table S1 (samples below 500x mean depth, by species) and Figure S2 (log10 depth heatmap)."""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram
HERE = os.path.dirname(os.path.abspath(__file__))
w = pd.read_csv(f"{HERE}/sample_by_locus_depth.csv"); L = [c for c in w.columns if c.startswith("AGAP")]
w["mean"] = w[L].mean(axis=1)
t = w.groupby("species").agg(n_below=("mean", lambda x: int((x < 500).sum())), n_total=("mean", "size"))
t["pct"] = (100 * t.n_below / t.n_total).round(1); t = t.sort_values("pct", ascending=False)
t.to_csv(f"{HERE}/tableS1_below500.tsv", sep="\t"); print(t)
SP = ["arabiensis", "coluzzii", "coustani", "funestus", "gambiae", "stephensi"]
sys.path.insert(0, HERE); from style import COL
M = np.log10(w[L].to_numpy() + 1)
order_l = dendrogram(linkage(M.T, "average", metric="correlation"), no_plot=True)["leaves"]
w["_sp"] = pd.Categorical(w.species, SP); idx = w.sort_values(["_sp", "mean"], ascending=[True, False]).index
cmap = LinearSegmentedColormap.from_list("b", ["#fcfcfb", "#d9d8d4", "#a8a7a2", "#7a7975", "#52514e", "#1f1f1e"])
fig = plt.figure(figsize=(7, 9), constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[0.25, 6, 0.25])
ax0, ax, cax = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
ax0.imshow(np.array([[matplotlib.colors.to_rgb(COL[s])] for s in w.loc[idx, "species"]]), aspect="auto")
ax0.set_xticks([]); ax0.set_yticks([])
im = ax.imshow(M[idx][:, order_l], aspect="auto", cmap=cmap, interpolation="nearest")
ax.set_xticks(range(len(L)), [L[i] for i in order_l], rotation=45, ha="right", fontsize=8); ax.set_yticks([])
ax0.set_ylabel(f"Samples (n={len(w)}), grouped by species")
cb = fig.colorbar(im, cax=cax); cb.set_label("log10(mean depth + 1)")
handles = [matplotlib.patches.Patch(color=COL[s], label=rf"$\it{{An.\ {s}}}$") for s in SP]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False, fontsize=8)
fig.savefig(f"{HERE}/figureS2_depth_heatmap.png", dpi=300, facecolor="white")
