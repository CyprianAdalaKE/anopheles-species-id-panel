#!/usr/bin/env python3
"""Largest species-pure clade per species (size, UFBoot) in a midpoint-rooted IQ-TREE consensus tree.
Tip labels end in _<species>. Usage: tree_pure_clades.py <contree>"""
import sys
from Bio import Phylo
t = Phylo.read(sys.argv[1], 'newick'); t.root_at_midpoint()
sp = lambda n: n.name.rsplit('_', 1)[1]
tot = {}
for tip in t.get_terminals(): tot[sp(tip)] = tot.get(sp(tip), 0) + 1
best = {}
for c in t.find_clades():
    tips = c.get_terminals(); s = {sp(x) for x in tips}
    if len(s) == 1 and len(tips) > 1:
        k = s.pop(); sup = c.confidence if c.confidence is not None else ''
        if len(tips) > best.get(k, (0, ''))[0]: best[k] = (len(tips), sup)
for k in sorted(tot): print(f"{k}\t{best.get(k, (1, ''))[0]}/{tot[k]}\tUFBoot {best.get(k, (0, ''))[1]}")
