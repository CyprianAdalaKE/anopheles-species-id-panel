"""Figure 13: ML tree (IQ-TREE consensus, UFBoot) with tips coloured by reference species, plus per-species
clade purity / support summary for the text."""
import os, sys
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from style import SP, COL, INK, INK2, it
T = os.path.join(os.path.dirname(HERE), "analysis_verified/phylo_tree")
sp_of = dict(l.strip().split("\t") for l in open(f"{T}/sample_species_map.tsv"))
tree = Phylo.read(f"{T}/snp_tree_nm5000.contree", "newick"); tree.root_at_midpoint(); tree.ladderize()
tips = tree.get_terminals()

# per-species: smallest clade with >=90% of that species' tips; purity and support
with open(f"{T}/clade_summary_nm5000.tsv", "w") as o:
    o.write("species\tn\tbest_clade_size\tspecies_in_clade\tpurity\tufboot\n")
    for s in SP:
        mine = {t.name for t in tips if sp_of.get(t.name) == s}; best = None
        for c in tree.find_clades(order="level"):
            names = {t.name for t in c.get_terminals()}; k = len(names & mine)
            if k >= 0.9 * len(mine) and (best is None or len(names) < best[0]): best = (len(names), k, c.confidence)
        o.write(f"{s}\t{len(mine)}\t{best[0]}\t{best[1]}\t{best[1]/best[0]:.2f}\t{best[2]}\n")
print(open(f"{T}/clade_summary_nm5000.tsv").read())

# layout
y = {}; n = [0]
def walk(c):
    if c.is_terminal(): n[0] += 1; y[c] = n[0]
    else:
        for ch in c.clades: walk(ch)
        y[c] = (y[c.clades[0]] + y[c.clades[-1]]) / 2
walk(tree.root); depth = tree.depths()
if not max(depth.values()): depth = tree.depths(unit_branch_lengths=True)
fig, ax = plt.subplots(figsize=(7.5, 10), constrained_layout=True)
for c in tree.find_clades():
    for ch in c.clades:
        ax.plot([depth[c], depth[c]], [y[ch], y[c]], color=INK2, lw=0.4)
        ax.plot([depth[c], depth[ch]], [y[ch], y[ch]], color=INK2, lw=0.4)
xmax = max(depth.values()); strip = xmax * 1.04
for t in tips:
    s = sp_of.get(t.name); ax.plot([strip, strip * 1.06], [y[t], y[t]], color=COL.get(s, INK), lw=2.2, solid_capstyle="butt")
ax.set_ylim(len(tips) + 1, 0); ax.set_yticks([]); ax.set_xlim(0, strip * 1.08)
for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
ax.set_xlabel("Substitutions per site")
ax.legend(handles=[Line2D([0], [0], color=COL[s], lw=6, label=it(s)) for s in SP], loc="lower left", frameon=False, fontsize=8.5)
fig.savefig(f"{HERE}/figure13_tree_nm5000.png", dpi=300, facecolor="white")
