#!/usr/bin/env python3
"""Per (run, barcode) species composition from per-read Kraken2 output.
Species = rank-S ancestor of the assigned taxon (reads assigned above species are 'unresolved').
Plasmodium = any taxon under genus Plasmodium (5820). Output composition.tsv (fractions of classified reads)."""
import collections
parent = {}; rank = {}; name = {}; stack = []
for l in open('subsample.kreport'):
    f = l.rstrip('\n').split('\t'); depth = (len(f[5]) - len(f[5].lstrip(' '))) // 2
    tid = int(f[4]); rank[tid] = f[3]; name[tid] = f[5].strip()
    stack = stack[:depth]; parent[tid] = stack[-1] if stack else 0; stack.append(tid)
def lineage(t):
    while t: yield t; t = parent.get(t, 0)
spmap = {}
def species_of(t):
    if t not in spmap:
        lin = list(lineage(t))
        s = next((x for x in lin if rank.get(x) == 'S'), None)
        spmap[t] = ('Plasmodium' if 5820 in lin else (name[s].replace('Anopheles ', '') if s else 'unresolved'))
    return spmap[t]
c = collections.defaultdict(collections.Counter); tot = collections.Counter()
for l in open('subsample.k2out'):
    st, rid, tid = l.split('\t')[:3]
    run, bc, _ = rid.split('|'); key = (run, bc); tot[key] += 1
    if st == 'C': c[key][species_of(int(tid))] += 1
cols = ['arabiensis', 'gambiae', 'coluzzii', 'funestus', 'stephensi', 'coustani', 'Plasmodium', 'unresolved', 'other']
with open('composition.tsv', 'w') as o:
    o.write('run\tbarcode\tsampled\tclassified\t' + '\t'.join('pct_' + x for x in cols) + '\tdominant\tdominant_pct\n')
    for key in sorted(tot):
        cl = sum(c[key].values()); row = c[key]
        vals = {x: row.get(x, 0) for x in cols[:-1]}; vals['other'] = cl - sum(vals.values())
        pct = {x: (100 * vals[x] / cl if cl else 0) for x in cols}
        anoph = {x: pct[x] for x in cols[:6]}; dom = max(anoph, key=anoph.get) if cl else '-'
        o.write(f'{key[0]}\t{key[1]}\t{tot[key]}\t{cl}\t' + '\t'.join(f'{pct[x]:.2f}' for x in cols) + f'\t{dom}\t{anoph.get(dom,0):.2f}\n')
