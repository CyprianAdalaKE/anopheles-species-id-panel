"""Confusion matrices (top Bracken call vs species label), from gamcol_validation/accuracy6/per_sample_calls.tsv.
A: six species. B: An. gambiae s.s. and An. coluzzii merged."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = "/home/genomics/Desktop/species_id/gamcol_validation/accuracy6/per_sample_calls.tsv"
rows = list(csv.DictReader(open(SRC), delimiter="\t"))

SP6 = ["arabiensis", "coluzzii", "coustani", "funestus", "gambiae", "stephensi"]
GC = "gambiae/coluzzii"
SP5 = ["arabiensis", GC, "coustani", "funestus", "stephensi"]
merge = lambda x: GC if x in ("gambiae", "coluzzii") else x


def matrix(species, f=lambda x: x):
    M = np.zeros((len(species),) * 2, int)
    for r in rows:
        M[species.index(f(r["true"])), species.index(f(r["pred"]))] += 1
    return M


M6, M5 = matrix(SP6), matrix(SP5, merge)
with open(f"{HERE}/FigX_confusion6_source_data.tsv", "w") as o:
    for name, sp, M in (("A_six_species", SP6, M6), ("B_gambiae_coluzzii_merged", SP5, M5)):
        o.write(f"# {name}\ntrue\\pred\t" + "\t".join(sp) + "\tn\n")
        for i, s in enumerate(sp):
            o.write(s + "\t" + "\t".join(map(str, M[i])) + f"\t{M[i].sum()}\n")

INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
cmap = LinearSegmentedColormap.from_list("blue", [SURF, "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})


def label(s):
    if s == GC:
        return r"$\it{An.\ gambiae}$ s.s./" + "\n" + r"$\it{An.\ coluzzii}$"
    return rf"$\it{{An.\ {s}}}$"


def panel(ax, M, species, tag):
    k = len(species); n = M.sum(1); P = 100 * M / n[:, None]
    im = ax.imshow(P, cmap=cmap, vmin=0, vmax=100, aspect="auto", interpolation="nearest")
    for i in range(k):
        for j in range(k):
            if M[i, j] == 0:
                ax.text(j, i, "0", ha="center", va="center", color="#a8a7a2", fontsize=8); continue
            dark = P[i, j] >= 45
            ax.text(j, i - 0.12, str(M[i, j]), ha="center", va="center", fontsize=10, fontweight="bold",
                    color="white" if dark else INK)
            ax.text(j, i + 0.2, f"{P[i, j]:.0f}%", ha="center", va="center", fontsize=7.5,
                    color="#e8eef7" if dark else INK2)
    ax.set_xticks(np.arange(-.5, k + .5), minor=True); ax.set_yticks(np.arange(-.5, k + .5), minor=True)
    ax.grid(which="minor", color="white", linewidth=2.5); ax.tick_params(which="minor", length=0)
    ax.set_xlim(-0.5, k - 0.5); ax.set_ylim(k - 0.5, -0.5)
    for v in (-0.5, k - 0.5):  # outer edges get the same white gap (covers imshow edge antialiasing)
        ax.axvline(v, color="white", lw=2.5, clip_on=False, zorder=5)
        ax.axhline(v, color="white", lw=2.5, clip_on=False, zorder=5)
    ax.set_xticks(range(k), [label(s) for s in species], rotation=35, ha="right", fontsize=8.5)
    ax.set_yticks(range(k), [f"{label(s)}\n(n={c})" if s == GC else f"{label(s)}  (n={c})"
                             for s, c in zip(species, n)], fontsize=8.5)
    ax.tick_params(length=0, colors=INK)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xlabel("Top Kraken2/Bracken call", color=INK)
    ax.set_ylabel("Species label (field/morphological ID)", color=INK, labelpad=6)
    ok = np.trace(M)
    ax.text(0, 1.025, tag, transform=ax.transAxes, fontsize=12, fontweight="bold", color=INK, ha="left", va="bottom")
    ax.text(0.06, 1.025, f"Overall accuracy {100*ok/M.sum():.1f}% ({ok}/{M.sum()})", transform=ax.transAxes,
            fontsize=9, color=INK2, ha="left", va="bottom")
    return im


fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4), constrained_layout=True,
                         gridspec_kw={"width_ratios": [6, 5]})
panel(axes[0], M6, SP6, "A")
im = panel(axes[1], M5, SP5, "B")
cb = fig.colorbar(im, ax=axes, shrink=0.7, pad=0.015); cb.set_label("% of samples in row", color=INK2)
cb.outline.set_visible(False); cb.ax.tick_params(colors=INK2, length=0)
for ext in ("png", "pdf"):
    fig.savefig(f"{HERE}/FigX_confusion6.{ext}", dpi=300, facecolor="white")
print(M6, M5, sep="\n\n")
