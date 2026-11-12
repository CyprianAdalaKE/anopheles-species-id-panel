#!/usr/bin/env python3
"""Per-sample read identity for coustani-labelled samples (5,000-read PAFs vs the six-species amplicon references).
For each read: best locus = most matching bases; identity to each species = matches/block at that locus (block>=300).
Classes: near-identical (>=97%) to a reference species = that species' DNA; best <90% = divergent (no panel reference).
Usage: coustani_identity.py <paf.gz> [...]  -> TSV on stdout"""
import gzip, sys, collections, os
print('sample\treads_aligned\tpct_ge97_' + '\tpct_ge97_'.join(['arabiensis','coluzzii','coustani','funestus','gambiae','stephensi']) + '\tpct_gambiae_complex_ge97\tpct_best_lt90\tmedian_best_identity')
for pf in sys.argv[1:]:
    hits = collections.defaultdict(list)
    for l in gzip.open(pf, 'rt'):
        x = l.split('\t'); m, blk = int(x[9]), int(x[10])
        if blk < 300: continue
        loc, s = x[5].split('|'); hits[x[0]].append((m, blk, loc, s))
    c = collections.Counter(); best = []
    for r, h in hits.items():
        loc = max(h)[2]; ident = {}
        for m, blk, l2, s in h:
            if l2 == loc: ident[s] = max(ident.get(s, 0), m / blk)
        b = max(ident.values()); best.append(b); top = [s for s, v in ident.items() if v == b]
        if b >= 0.97:
            if len(top) == 1: c[top[0]] += 1
            if set(top) <= {'arabiensis', 'coluzzii', 'gambiae'}: c['complex'] += 1
        if b < 0.90: c['lt90'] += 1
    n = len(best); best.sort()
    name = os.path.basename(pf).replace('.paf.gz', '')
    print(name + f'\t{n}\t' + '\t'.join(f"{100*c[s]/n:.1f}" for s in ['arabiensis','coluzzii','coustani','funestus','gambiae','stephensi']) +
          f"\t{100*c['complex']/n:.1f}\t{100*c['lt90']/n:.1f}\t{100*best[n//2]:.1f}")
