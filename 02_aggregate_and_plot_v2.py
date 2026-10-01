#!/usr/bin/env python3
"""
02_aggregate_and_plot.py (v2 - species-aware)

Reads all mosdepth *.regions.bed.gz files from coverage_results/, parses out
species and sample name from each filename, and builds:
  1. sample_by_locus_depth.csv       -- full matrix with species labels
  2. coverage_heatmap.png            -- samples (grouped/colored by species) x loci
  3. coverage_barchart_overall.png   -- mean depth per locus, across all samples
  4. coverage_barchart_by_species.png -- mean depth per locus, faceted by species
     (this is the figure most likely to reveal species-specific dropout)

Requires: pandas, seaborn, matplotlib, numpy
    pip install pandas seaborn matplotlib numpy --break-system-packages
"""

import glob
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ---- EDIT THESE ----
PANEL_DIR = os.environ.get("PANEL_DIR", os.path.dirname(os.path.abspath(__file__)))
COVERAGE_DIR = os.path.join(PANEL_DIR, "coverage_results")
OUT_DIR = os.path.join(PANEL_DIR, "figures")
# ----------------------

os.makedirs(OUT_DIR, exist_ok=True)

# Known species prefixes (used to split "species_sample" filenames correctly,
# since some sample names could themselves contain underscores)
KNOWN_SPECIES = ["gambiae", "arabiensis", "coluzzii", "funestus", "coustani", "stephensi"]


def parse_filename(path):
    """Extract species and normalized sample name from a mosdepth output filename."""
    fname = os.path.basename(path).replace(".regions.bed.gz", "")

    species = None
    sample = fname
    for sp in KNOWN_SPECIES:
        if fname.startswith(sp + "_"):
            species = sp
            sample = fname[len(sp) + 1:]
            break

    if species is None:
        # Fallback: split on first underscore
        parts = fname.split("_", 1)
        species = parts[0]
        sample = parts[1] if len(parts) > 1 else fname

    # Normalize known suffix quirks (e.g. "barcode59.fastq1" -> "barcode59")
    sample = re.sub(r"\.fastq\d*$", "", sample)

    return species, sample


# ---- Load all files ----
records = []
files = glob.glob(os.path.join(COVERAGE_DIR, "*.regions.bed.gz"))

if not files:
    raise FileNotFoundError(f"No coverage files found in {COVERAGE_DIR}")

for f in files:
    species, sample = parse_filename(f)
    df = pd.read_csv(
        f, sep="\t", header=None,
        names=["chrom", "start", "end", "locus", "mean_depth"]
    )
    df["species"] = species
    df["sample"] = sample
    records.append(df[["species", "sample", "locus", "mean_depth"]])

long_df = pd.concat(records, ignore_index=True)

# Combine species+sample into a unique row label (some barcode numbers repeat across species)
long_df["sample_label"] = long_df["species"] + "_" + long_df["sample"]

print(f"Loaded {long_df['sample_label'].nunique()} samples across "
      f"{long_df['species'].nunique()} species and {long_df['locus'].nunique()} loci")
print(long_df.groupby("species")["sample"].nunique())

# ---- Pivot into sample x locus matrix (keep species mapping separately) ----
species_map = long_df.drop_duplicates("sample_label").set_index("sample_label")["species"]
matrix = long_df.pivot(index="sample_label", columns="locus", values="mean_depth")

# Sort rows by species so the heatmap groups naturally
matrix = matrix.loc[species_map.sort_values().index]

# Save full matrix with species column for supplementary data
matrix_with_species = matrix.copy()
matrix_with_species.insert(0, "species", species_map)
matrix_with_species.to_csv(os.path.join(OUT_DIR, "sample_by_locus_depth.csv"))
print(f"\nSaved matrix: {OUT_DIR}/sample_by_locus_depth.csv")

# =====================================================================
# FIGURE 1: Heatmap, rows grouped and color-annotated by species
# =====================================================================
species_sorted = species_map.loc[matrix.index]
palette = sns.color_palette("Set2", n_colors=species_sorted.nunique())
species_colors = dict(zip(sorted(species_sorted.unique()), palette))
row_colors = species_sorted.map(species_colors)

log_matrix = np.log10(matrix + 1)

g = sns.clustermap(
    log_matrix,
    row_cluster=False,      # keep species-grouped order, don't reorder by similarity
    col_cluster=True,       # ok to cluster loci by co-behavior
    row_colors=row_colors,
    cmap="viridis",
    figsize=(max(10, 0.6 * matrix.shape[1]), max(10, 0.15 * matrix.shape[0])),
    cbar_kws={"label": "log10(mean depth + 1)"},
    linewidths=0.2,
)
g.ax_heatmap.set_xlabel("Locus")
g.ax_heatmap.set_ylabel("Sample (grouped by species)")

# Build a legend for species colors
from matplotlib.patches import Patch
handles = [Patch(facecolor=species_colors[sp], label=sp) for sp in sorted(species_colors)]
g.ax_heatmap.legend(
    handles=handles, title="Species",
    bbox_to_anchor=(1.25, 1), loc="upper left", borderaxespad=0.
)
g.savefig(os.path.join(OUT_DIR, "coverage_heatmap.png"), dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT_DIR}/coverage_heatmap.png")

# =====================================================================
# FIGURE 2: Overall mean depth per locus (across all samples/species)
# =====================================================================
locus_summary = matrix.mean(axis=0).sort_values(ascending=False)
locus_sd = matrix.std(axis=0).reindex(locus_summary.index)

plt.figure(figsize=(max(8, 0.6 * len(locus_summary)), 6))
plt.bar(
    locus_summary.index, locus_summary.values,
    yerr=locus_sd.values, capsize=4,
    color=sns.color_palette("crest", len(locus_summary))
)
plt.ylabel("Mean depth (reads), \u00b1 SD across all samples")
plt.xlabel("Locus")
plt.title("Mean per-locus depth across all samples (all species combined)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "coverage_barchart_overall.png"), dpi=300)
plt.close()
print(f"Saved: {OUT_DIR}/coverage_barchart_overall.png")

# =====================================================================
# FIGURE 3: Mean depth per locus, faceted by species
# This reveals whether low coverage is locus-specific, species-specific,
# or both -- directly useful for a Results statement about which loci/
# species combinations are less reliable.
# =====================================================================
species_locus_mean = long_df.groupby(["species", "locus"])["mean_depth"].mean().reset_index()

locus_order = matrix.columns.tolist()  # keep consistent locus order across panels

g2 = sns.catplot(
    data=species_locus_mean, x="locus", y="mean_depth", col="species",
    order=locus_order,
    kind="bar", col_wrap=3, height=3.5, aspect=1.3,
    palette="crest", hue="locus", hue_order=locus_order, legend=False,
    sharey=False
)
for ax in g2.axes.flat:
    ax.tick_params(axis="x", labelrotation=45)
    for label in ax.get_xticklabels():
        label.set_ha("right")
g2.set_axis_labels("Locus", "Mean depth")
g2.fig.suptitle("Mean per-locus depth by species", y=1.02)
g2.savefig(os.path.join(OUT_DIR, "coverage_barchart_by_species.png"), dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT_DIR}/coverage_barchart_by_species.png")

# ---- Summary stats for your Results text ----
print("\n--- Summary for Results text ---")
print(f"Total samples: {matrix.shape[0]}")
print(f"Total loci: {matrix.shape[1]}")
print(f"Overall mean depth: {matrix.values.mean():.1f}")
print(f"Highest-depth locus: {locus_summary.index[0]} ({locus_summary.iloc[0]:.1f}x)")
print(f"Lowest-depth locus: {locus_summary.index[-1]} ({locus_summary.iloc[-1]:.1f}x)")
print(f"Fold difference (max/min): {locus_summary.iloc[0] / locus_summary.iloc[-1]:.1f}x")
print("\nPer-species mean depth (averaged across all loci):")
print(long_df.groupby("species")["mean_depth"].mean().sort_values(ascending=False).round(1))
