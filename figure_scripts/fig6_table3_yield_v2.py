"""Table 3 and Figure 6: read yield and mean read length per species group (species_id_verified_fastqs),
plus the infectivity set (combined_stats.tsv)."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); BASE = os.path.dirname(HERE)
SP = ["arabiensis", "gambiae", "funestus", "coluzzii", "coustani", "stephensi"]
d = {s: [] for s in SP + ["infectivity"]}
for r in csv.DictReader(open(f"{BASE}/verified_stats/read_stats.tsv"), delimiter="\t"):
    d[r["file"].split("/")[-2]].append((int(r["num_seqs"]), int(r["sum_len"])))
# infectivity set = SLCA014 BatchINF barcodes 01-50, counted from the raw MinKNOW fastq_pass chunks
for r in csv.DictReader(open(f"{BASE}/verified_stats/infectivity_BatchINF_raw_stats.tsv"), delimiter="\t"):
    d["infectivity"].append((int(r["num_seqs"]), int(r["sum_len"])))
N_INF = len(d["infectivity"])
with open(f"{HERE}/table3_read_yield_v2.tsv", "w") as o:
    o.write("group\tn_samples\ttotal_reads\tmean_reads_per_sample\tmean_read_length_bp\n")
    for g in SP + ["infectivity"]:
        n = N_INF if g == "infectivity" else len(d[g]); r = sum(x[0] for x in d[g]); b = sum(x[1] for x in d[g])
        o.write(f"{g}\t{n}\t{r}\t{r/n:.1f}\t{b/r:.1f}\n")
print(open(f"{HERE}/table3_read_yield_v2.tsv").read())

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "xtick.color": INK, "ytick.color": INK2})
order = ["arabiensis", "gambiae", "coluzzii", "funestus", "stephensi", "coustani", "infectivity"]
labs = [rf"$\it{{An.}}$" + "\n" + rf"$\it{{{g}}}$" if g != "infectivity" else "Infectivity\nset" for g in order]
rng = np.random.default_rng(1)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
for ax, f, yl, log in ((axes[0], lambda x: x[0], "Reads per sample (log scale)", True),
                       (axes[1], lambda x: x[1] / x[0], "Mean read length (bp)", False)):
    data = [np.array([f(x) for x in d[g] if x[0] > 0]) for g in order]
    ax.boxplot(data, positions=range(len(order)), widths=0.6, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="#cde2fb", edgecolor="#2a78d6"), medianprops=dict(color="#2a78d6", lw=2),
               whiskerprops=dict(color="#2a78d6"), capprops=dict(color="#2a78d6"))
    for i, v in enumerate(data):
        ax.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v, s=5, color=INK, alpha=0.3, lw=0, zorder=3)
    if log: ax.set_yscale("log")
    ax.set_xticks(range(len(order)), [f"{l}\n(n={len(v)})" for l, v in zip(labs, data)], fontsize=8)
    ax.set_ylabel(yl, color=INK); ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
for ax, t in zip(axes, "AB"):
    ax.text(0, 1.02, t, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
fig.savefig(f"{HERE}/figure6_yield_readlength_v2.png", dpi=300, facecolor="white")
