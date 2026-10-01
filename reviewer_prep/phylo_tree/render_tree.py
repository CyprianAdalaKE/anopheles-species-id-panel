#!/usr/bin/env python3
"""
Render the IQ-TREE ML tree as a circular/rectangular cladogram with tips
colored by species (no tip-label text -- 471 taxa is too dense to label
individually; color alone shows whether species form monophyletic groups,
directly complementing the genotype PCA).
"""
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPECIES_COLORS = {
    "gambiae":    "#1b9e77",
    "coluzzii":   "#d95f02",
    "arabiensis": "#7570b3",
    "funestus":   "#e7298a",
    "stephensi":  "#66a61e",
    "coustani":   "#e6ab02",
}

species_map = {}
with open("sample_species_map.tsv") as f:
    for line in f:
        label, sp = line.strip().split("\t")
        species_map[label] = sp

tree = Phylo.read("snp_tree.contree", "newick")
tree.root_at_midpoint()
tree.ladderize()

terminals = tree.get_terminals()
print(f"{len(terminals)} tips")

fig, ax = plt.subplots(figsize=(10, max(14, len(terminals) * 0.035)))

def get_color(clade):
    if clade.name and clade.name in species_map:
        return SPECIES_COLORS.get(species_map[clade.name], "black")
    return "black"

Phylo.draw(
    tree, axes=ax, do_show=False,
    label_func=lambda c: "",  # no tip labels -- too dense at n=471
    branch_labels=None,
)

# color tip points manually via the drawn axes
for term in terminals:
    sp = species_map.get(term.name)
    color = SPECIES_COLORS.get(sp, "black")
    y = term.__dict__.get("_y", None)

# Bio.Phylo.draw doesn't expose coordinates easily; use a simpler manual approach instead
plt.close()

# ---- manual rectangular layout with colored tip markers ----
fig, ax = plt.subplots(figsize=(9, max(14, len(terminals) * 0.032)))
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: "", show_confidence=False)

# find tip y-positions by re-walking the same order Bio.Phylo uses internally
y_positions = {}
y = 0
def assign_y(clade):
    global y
    if clade.is_terminal():
        y += 1
        y_positions[clade] = y
    else:
        for child in clade.clades:
            assign_y(child)
assign_y(tree.root)

depths = tree.depths()
for term in terminals:
    x = depths[term]
    yy = y_positions[term]
    sp = species_map.get(term.name)
    color = SPECIES_COLORS.get(sp, "black")
    ax.plot(x, yy, "o", color=color, markersize=3, zorder=5)

ax.set_yticks([])
ax.set_ylabel("")
ax.set_xlabel("Substitutions/site")
ax.set_title("Maximum-likelihood tree from panel SNP genotypes (471 samples, 1,022 SNPs)")

from matplotlib.lines import Line2D
handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=8, label=sp)
           for sp, c in SPECIES_COLORS.items()]
ax.legend(handles=handles, title="Species", loc="lower right", fontsize=9)

plt.tight_layout()
plt.savefig("phylo_tree_species_colored.png", dpi=250)
print("Saved phylo_tree_species_colored.png")
