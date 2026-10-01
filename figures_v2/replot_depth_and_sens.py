#!/usr/bin/env python3
"""
Replace Figures 8, 9 (per-locus depth, currently mean +/- SD capped at zero)
with log-scale boxplots (median/IQR + individual points, matching the
convention already used elsewhere in the manuscript e.g. Figure 16/Table 4),
and Figure 12 (sensitivity/specificity) with Wilson 95% CI error bars added.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

OUT = "."

# ---- QC exclusions matching the manuscript's 478 true-biological-sample set ----
EXCLUDE = {
    ("coluzzii", "barcode80"),
    ("coustani", "barcode84"),
    ("coustani", "barcode95"),
}

df = pd.read_csv("sample_by_locus_depth.csv")
df["barcode"] = df["sample_label"].str.extract(r"(barcode\d+)$")
mask = ~df.apply(lambda r: (r["species"], r["barcode"]) in EXCLUDE, axis=1)
df = df[mask].drop(columns=["barcode"])
print(f"Samples after QC exclusion: {len(df)}")

locus_cols = [c for c in df.columns if c.startswith("AGAP")]
long_df = df.melt(id_vars=["sample_label", "species"], value_vars=locus_cols,
                   var_name="locus", value_name="depth")

locus_order = df[locus_cols].mean().sort_values(ascending=False).index.tolist()

# =====================================================================
# Figure 8 replacement: overall per-locus depth, log-scale boxplot
# =====================================================================
plt.figure(figsize=(10, 6))
ax = sns.boxplot(data=long_df, x="locus", y="depth", order=locus_order,
                  palette="crest", showfliers=False, linewidth=1)
sns.stripplot(data=long_df, x="locus", y="depth", order=locus_order,
              color="black", alpha=0.25, size=2, jitter=0.25, ax=ax)
ax.set_yscale("log")
ax.set_ylabel("Sequencing depth (reads, log scale)")
ax.set_xlabel("Locus")
ax.set_title("Per-locus sequencing depth across all samples (all species combined)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(f"{OUT}/figure8_depth_overall_logscale.png", dpi=300)
plt.close()
print("Saved figure8_depth_overall_logscale.png")

# =====================================================================
# Figure 9 replacement: per-locus depth by species, log-scale boxplot
# =====================================================================
g = sns.catplot(
    data=long_df, x="locus", y="depth", col="species",
    order=locus_order, kind="box", col_wrap=3, height=3.5, aspect=1.3,
    palette="crest", showfliers=False, sharey=False,
)
for ax in g.axes.flat:
    ax.set_yscale("log")
    ax.tick_params(axis="x", labelrotation=45)
    for label in ax.get_xticklabels():
        label.set_ha("right")
g.set_axis_labels("Locus", "Depth (reads, log scale)")
g.fig.suptitle("Per-locus sequencing depth by species (log scale)", y=1.02)
g.savefig(f"{OUT}/figure9_depth_by_species_logscale.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved figure9_depth_by_species_logscale.png")

# =====================================================================
# Figure 12 replacement: sensitivity/specificity with Wilson 95% CI error bars
# =====================================================================
with open(f"{OUT}/wilson_ci_table3.json") as f:
    wilson = json.load(f)

species_order = ["arabiensis", "coluzzii", "coustani", "funestus", "gambiae", "stephensi"]
sens_vals, sens_err_lo, sens_err_hi = [], [], []
spec_vals, spec_err_lo, spec_err_hi = [], [], []
for sp in species_order:
    n, sens, sens_lo, sens_hi, spec, spec_lo, spec_hi, *_ = wilson[sp]
    sens_vals.append(sens); sens_err_lo.append(sens - sens_lo); sens_err_hi.append(sens_hi - sens)
    spec_vals.append(spec); spec_err_lo.append(spec - spec_lo); spec_err_hi.append(spec_hi - spec)

x = np.arange(len(species_order))
width = 0.35
fig, ax = plt.subplots(figsize=(12, 7))
ax.bar(x - width/2, sens_vals, width, yerr=[sens_err_lo, sens_err_hi], capsize=4,
       label="Sensitivity", color="#3b6a8f")
ax.bar(x + width/2, spec_vals, width, yerr=[spec_err_lo, spec_err_hi], capsize=4,
       label="Specificity", color="#4fa07a")
ax.set_xticks(x)
ax.set_xticklabels([f"Anopheles {sp}" for sp in species_order], rotation=45, ha="right")
ax.set_ylabel("Score (%)")
ax.set_ylim(0, 105)
ax.set_title("Sensitivity and Specificity per Species (with Wilson 95% CI)")
ax.legend(title="Metric", bbox_to_anchor=(1.02, 1), loc="upper left")
ax.grid(axis="y", linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig(f"{OUT}/figure12_sens_spec_with_CI.png", dpi=300)
plt.close()
print("Saved figure12_sens_spec_with_CI.png")
