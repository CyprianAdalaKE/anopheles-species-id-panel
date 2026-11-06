"""Figures 8 and 9 regenerated (2026-11-04) from species_id_verified_fastqs (374 samples).
Reads mosdepth regions files from coverage_results_verified/ (align_coverage_verified.sh) and writes
sample_by_locus_depth.csv, figure7_depth_overall.{png,pdf} and
figure8_depth_by_species.{png,pdf}. Same encoding as replot_depth_and_sens.py
(log-scale boxplot of per-sample mean depth with individual points); figure 9 shares one y-axis across species."""
import glob, gzip, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
COV = "/home/genomics/Desktop/species_id/anopheles-panel/coverage_results_verified"
SP = ["arabiensis", "coluzzii", "coustani", "funestus", "gambiae", "stephensi"]

rec = []
for f in glob.glob(f"{COV}/*.regions.bed.gz"):
    name = os.path.basename(f).removesuffix(".regions.bed.gz")
    sp, sample = name.split("_", 1)
    for line in gzip.open(f, "rt"):
        c = line.rstrip("\n").split("\t")
        rec.append((sp, f"{sp}_{sample}", c[3], float(c[4])))
long = pd.DataFrame(rec, columns=["species", "sample_label", "locus", "depth"])
wide = long.pivot_table(index=["sample_label", "species"], columns="locus", values="depth").reset_index()
wide.to_csv(f"{HERE}/sample_by_locus_depth.csv", index=False)
order = long.groupby("locus")["depth"].mean().sort_values(ascending=False).index.tolist()
print(wide.groupby("species").size().to_dict(), "samples;", len(order), "loci")

BOX, PT, INK, INK2, GRID = "#2a78d6", "#0b0b0b", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2,
                     "xtick.color": INK, "ytick.color": INK2})
rng = np.random.default_rng(1)
FLOOR = 0.1  # zero-depth samples are drawn at the floor of the log axis


def boxes(ax, d):
    data = [np.maximum(d.loc[d.locus == l, "depth"].to_numpy(), FLOOR) for l in order]
    ax.boxplot(data, positions=range(len(order)), widths=0.6, showfliers=False, patch_artist=True, whis=(5, 95),
               boxprops=dict(facecolor="#cde2fb", edgecolor=BOX, lw=1), medianprops=dict(color=BOX, lw=2),
               whiskerprops=dict(color=BOX, lw=1), capprops=dict(color=BOX, lw=1))
    for i, v in enumerate(data):
        ax.scatter(i + rng.uniform(-0.22, 0.22, len(v)), v, s=3, color=PT, alpha=0.25, lw=0, zorder=3)
    ax.set_yscale("log")
    ax.set_xticks(range(len(order)), order, rotation=45, ha="right", fontsize=8)
    ax.yaxis.grid(True, color=GRID, lw=0.8, which="major"); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)


fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
boxes(ax, long)
ax.set_ylabel("Mean depth per sample (reads, log scale)", color=INK); ax.set_xlabel("Locus", color=INK)
pass
for ext in ("png", "pdf"):
    fig.savefig(f"{HERE}/figure7_depth_overall.{ext}", dpi=300, facecolor="white")

fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharey=True, constrained_layout=True)
for ax, sp in zip(axes.flat, SP):
    d = long[long.species == sp]
    boxes(ax, d)
    ax.set_title(rf"$\it{{An.\ {sp}}}$ (n={d.sample_label.nunique()})", loc="left", fontsize=10, color=INK)
for ax in axes[:, 0]: ax.set_ylabel("Mean depth (reads, log scale)", color=INK)
for ext in ("png", "pdf"):
    fig.savefig(f"{HERE}/figure8_depth_by_species.{ext}", dpi=300, facecolor="white")
print("zero-depth sample-locus cells:", int((long.depth == 0).sum()))
