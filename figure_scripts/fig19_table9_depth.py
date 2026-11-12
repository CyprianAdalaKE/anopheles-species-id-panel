"""Table 9 / Figure 19: pooled top-call accuracy vs target depth (48 samples, 8 per species; Wilson 95% CI)."""
import csv, math, os, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), "analysis_verified/depth_accuracy/depth_experiment_results.tsv")
def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z*z/n; c = (p + z*z/(2*n)) / d; h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return 100*(c-h), 100*min(1, c+h)
r = [x for x in csv.DictReader(open(SRC), delimiter="\t") if x["species"] != "ALL DONE"]
by = collections.defaultdict(list)
for x in r: by[int(x["target_depth"])].append(x["correct"] == "TRUE")
D = sorted(by); acc = [100*sum(by[d])/len(by[d]) for d in D]; ci = [wilson(sum(by[d]), len(by[d])) for d in D]
with open(f"{HERE}/table9_depth_accuracy.tsv", "w") as o:
    o.write("target_depth\tn\tcorrect\taccuracy_%\t95CI\n")
    for d, a, c in zip(D, acc, ci): o.write(f"{d}\t{len(by[d])}\t{sum(by[d])}\t{a:.1f}\t{c[0]:.1f}-{c[1]:.1f}\n")
print(open(f"{HERE}/table9_depth_accuracy.tsv").read())
INK, INK2, GRID, C = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2})
fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
lo = [a - c[0] for a, c in zip(acc, ci)]; hi = [c[1] - a for a, c in zip(acc, ci)]
ax.errorbar(D, acc, yerr=[lo, hi], color=C, lw=2, marker="o", ms=6, mfc=C, mec="white", mew=1.5, capsize=3, ecolor=INK2, elinewidth=1)
for d, a in zip(D, acc): ax.text(d * 1.08, a + 1.2, f"{a:.1f}%", ha="left", va="bottom", fontsize=7.5, color=INK)
ax.set_xscale("log"); ax.set_xticks(D, [f"{d:,}" for d in D]); ax.minorticks_off()
ax.set_ylim(0, 100); ax.set_xlabel("Target sequencing depth (reads per sample, log scale)"); ax.set_ylabel("Pooled classification accuracy (%)")
ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.savefig(f"{HERE}/figure19_depth_accuracy.png", dpi=300, facecolor="white")
