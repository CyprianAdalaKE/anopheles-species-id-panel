#!/usr/bin/env python3
"""
Confusion matrix, per-species sensitivity/specificity (with Wilson 95% CIs) and
overall accuracy of panel-based species identification (Figure 9, Table 4, Figure 10).

Input: CSV with one row per sample and two columns:
  panel                 -- species called by the panel (Kraken2/Bracken)
  reference_species_id  -- reference (morphological/molecular) species identification

Portable version of crosstab_original_colab.ipynb (same one-vs-rest calculation).

Usage: python confusion_matrix_sens_spec.py species_id_joined_k2_Bracken_updated_Jul26_latest.csv [out_dir]
"""
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from statsmodels.stats.proportion import proportion_confint

in_csv = sys.argv[1]
out_dir = sys.argv[2] if len(sys.argv) > 2 else "."
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv(in_csv)
crosstab_df = pd.crosstab(df["reference_species_id"], df["panel"])
crosstab_df.to_csv(os.path.join(out_dir, "confusion_matrix.tsv"), sep="\t")
print(crosstab_df)

# Confusion matrix heatmap (Figure 9)
plt.figure(figsize=(12, 10))
sns.heatmap(crosstab_df.fillna(0).astype(int), annot=True, fmt="d", cmap="viridis", linewidths=.5)
plt.title("Heatmap of Reference Species vs Panel Counts")
plt.xlabel("Panel")
plt.ylabel("Reference Species ID")
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "crosstab_heatmap.png"), dpi=300)
plt.close()

# One-vs-rest sensitivity and specificity per species (Table 4)
total_sum = crosstab_df.sum().sum()
rows = []
for species in crosstab_df.index:
    tp = crosstab_df.loc[species, species] if species in crosstab_df.columns else 0
    fn = crosstab_df.loc[species, :].sum() - tp
    fp = (crosstab_df.loc[:, species].sum() if species in crosstab_df.columns else 0) - tp
    tn = total_sum - tp - fn - fp
    sens_lo, sens_hi = proportion_confint(tp, tp + fn, method="wilson")
    spec_lo, spec_hi = proportion_confint(tn, tn + fp, method="wilson")
    rows.append({
        "species": species, "TP": tp, "FN": fn, "FP": fp, "TN": tn,
        "sensitivity_pct": 100 * tp / (tp + fn),
        "sensitivity_ci_low": 100 * sens_lo, "sensitivity_ci_high": 100 * sens_hi,
        "specificity_pct": 100 * tn / (tn + fp),
        "specificity_ci_low": 100 * spec_lo, "specificity_ci_high": 100 * spec_hi,
    })
results = pd.DataFrame(rows)
results.to_csv(os.path.join(out_dir, "sensitivity_specificity.tsv"), sep="\t", index=False, float_format="%.1f")
print(results.round(1).to_string(index=False))

correct = sum(crosstab_df.loc[s, s] for s in crosstab_df.index if s in crosstab_df.columns)
print(f"\nOverall accuracy: {correct}/{total_sum} = {100 * correct / total_sum:.1f}%")

# Sensitivity/specificity bar chart (as in the original notebook; the CI version
# used in the manuscript is figures_v2/replot_depth_and_sens.py)
melted = results.melt(id_vars="species", value_vars=["sensitivity_pct", "specificity_pct"],
                      var_name="Metric", value_name="Value")
melted["Metric"] = melted["Metric"].map({"sensitivity_pct": "Sensitivity", "specificity_pct": "Specificity"})
plt.figure(figsize=(12, 7))
sns.barplot(x="species", y="Value", hue="Metric", data=melted, palette="viridis")
plt.title("Sensitivity and Specificity per Species")
plt.xlabel("Species")
plt.ylabel("Score (%)")
plt.xticks(rotation=45, ha="right")
plt.ylim(0, 100)
plt.legend(title="Metric", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "sensitivity_specificity_bar_graph.png"), dpi=300)
