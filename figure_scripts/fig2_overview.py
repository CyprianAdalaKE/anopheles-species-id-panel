"""Fig. 2: whole-study workflow overview (vector schematic)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

INK, INK2 = "#1f1f1d", "#4a4945"
LAB_F, LAB_E = "#e3eefb", "#2a78d6"
SEQ_F, SEQ_E = "#e6f5ef", "#1baf7a"
OUT_F, OUT_E = "#fdeee6", "#eb6834"
plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.default": "regular"})
AN, PL, PF = r"$\it{Anopheles}$", r"$\it{Plasmodium}$", r"$\it{P.\ falciparum}$"
fig, ax = plt.subplots(figsize=(11, 4.3), dpi=300)
ax.set_xlim(0, 110); ax.set_ylim(0, 42); ax.axis("off")

def box(x, y, w, h, title, lines, fc, ec, num):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.35,rounding_size=1.4", fc=fc, ec=ec, lw=1.5))
    ax.text(x + 1.0, y + h - 0.9, num, ha="left", va="top", fontsize=8, fontweight="bold", color=ec)
    ax.text(x + w / 2, y + h - 1.0, title, ha="center", va="top", fontsize=9.5, fontweight="bold", color=INK)
    ax.text(x + w / 2, y + h - 5.2, "\n".join(lines), ha="center", va="top", fontsize=7.7, color=INK2, linespacing=1.45)

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=13, lw=1.4, color=INK2, shrinkA=0, shrinkB=0))

# top row
W, H, Y = 23.5, 13.5, 27.5
xs = [1.0, 28.3, 55.6, 82.9]
box(xs[0], Y, W, H, "Specimen", ["Field-collected mosquito", "(374 specimens, six species)", "Genomic DNA extraction", "(salt–alcohol precipitation)"], LAB_F, LAB_E, "1")
box(xs[1], Y, W, H, "One multiplex PCR", ["11 targets in a single reaction:", f"10 {AN} nuclear loci (811–1,051 bp)", f"+ {PL} 18S rDNA (846 bp)"], LAB_F, LAB_E, "2")
box(xs[2], Y, W, H, "Barcoding and library", ["Bead clean-up and quantification", "Native Barcoding Kit 96 V14", "(SQK-NBD114.96),", "up to 96 samples per pool"], LAB_F, LAB_E, "3")
box(xs[3], Y, W, H, "MinION sequencing", ["R10.4.1 flow cell, MinION Mk1B", "GPU laptop, super-accurate", "basecalling", "≈US$21 per sample consumables"], SEQ_F, SEQ_E, "4")
for i in range(3):
    arrow(xs[i] + W + 0.5, Y + H / 2, xs[i + 1] - 0.5, Y + H / 2)

# reads bar
BY, BH = 19.5, 4.4
ax.add_patch(FancyBboxPatch((1.0, BY), 105.4, BH, boxstyle="round,pad=0.35,rounding_size=1.4", fc=SEQ_F, ec=SEQ_E, lw=1.5))
ax.text(53.7, BY + BH / 2, "Demultiplexed long reads (FASTQ), one file per specimen", ha="center", va="center", fontsize=8.8, fontweight="bold", color=INK)
arrow(xs[3] + W / 2, Y - 0.4, xs[3] + W / 2, BY + BH + 0.4)

# bottom row
Wb, Hb, Yb = 33.0, 11.5, 1.0
xb = [1.0, 37.2, 73.4]
box(xb[0], Yb, Wb, Hb, "Species identification", ["Kraken2 + Bracken against six", f"{AN} and four {PL} genomes", f"→ top {AN} species per specimen"], OUT_F, OUT_E, "5")
box(xb[1], Yb, Wb, Hb, f"{PL} detection", [f"Kraken2 + Bracken → % reads {PF}", "Competitive alignment to host and", "parasite references"], OUT_F, OUT_E, "6")
box(xb[2], Yb, Wb, Hb, "Coverage and genotype checks", ["minimap2 → per-locus depth", "SNP genotyping (bcftools)", "→ PCA, Fst, discriminant analysis, tree"], OUT_F, OUT_E, "7")
for x in xb:
    arrow(x + Wb / 2, BY - 0.4, x + Wb / 2, Yb + Hb + 0.4)

fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)
fig.savefig("figure2_overview.png", facecolor="white"); fig.savefig("figure2_overview.pdf", facecolor="white")
