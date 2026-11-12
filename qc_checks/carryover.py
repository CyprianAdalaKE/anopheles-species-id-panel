#!/usr/bin/env python3
"""Barcode-specific carry-over test on reused flow cells.
For a later run L, species S and barcode b: f = % of b's classified reads assigned to S in L;
prior load = estimated S reads on barcode b in all earlier runs on the same flow cell
(raw reads on b x fraction S in the subsample). If library left in the flow cell carries over,
f should rise with the prior load on the *same* barcode. Barcodes where S is the dominant call in L
are excluded (S is the sample itself). Spearman rho, one-sided permutation p (barcode labels shuffled).
Writes carryover_tests.tsv and carryover_per_sample.tsv."""
import csv, collections, random
from scipy.stats import spearmanr
runs = [r for r in csv.DictReader(open('runs.tsv'), delimiter='\t')]
order = {r['label']: (r['flow_cell'], r['started']) for r in runs}
comp = {(r['run'], r['barcode']): r for r in csv.DictReader(open('composition.tsv'), delimiter='\t')}
raw = {(r['run'], r['barcode']): int(r['raw_reads']) for r in csv.DictReader(open('raw_barcode_counts.tsv'), delimiter='\t')}
SP = ['arabiensis', 'gambiae', 'coluzzii', 'funestus', 'stephensi', 'coustani', 'Plasmodium']
def load(run, bc, s):
    fc, t = order[run]; tot = 0.0
    for p, (pfc, pt) in order.items():
        if pfc == fc and pt < t and (p, bc) in comp:
            tot += raw.get((p, bc), 0) * float(comp[(p, bc)]['pct_' + s]) / 100
    return tot
random.seed(1); tests = []; per = []
for L, (fc, t) in sorted(order.items(), key=lambda x: (x[1][0], x[1][1])):
    if not any(pfc == fc and pt < t for pfc, pt in order.values()): continue
    bcs = [b for (r, b), c in comp.items() if r == L and b != 'unclassified' and int(c['classified']) >= 500]
    for s in SP:
        xs, ys = [], []
        for b in bcs:
            c = comp[(L, b)]
            if c['dominant'] == s: continue
            xs.append(load(L, b, s)); ys.append(float(c['pct_' + s]))
        if len(xs) < 8 or max(xs) == 0: continue
        rho = spearmanr(xs, ys).statistic
        perm = sum(spearmanr(xs, random.sample(ys, len(ys))).statistic >= rho for _ in range(2000))
        hi = sorted(zip(xs, ys)); half = len(hi) // 2
        lo_m = sorted(y for _, y in hi[:half])[half // 2]; hi_m = sorted(y for _, y in hi[half:])[(len(hi) - half) // 2]
        tests.append([fc, L, s, len(xs), f'{rho:.3f}', f'{(perm + 1) / 2001:.4f}', f'{lo_m:.2f}', f'{hi_m:.2f}'])
with open('carryover_tests.tsv', 'w') as o:
    o.write('flow_cell\trun\tspecies\tn_barcodes\tspearman_rho\tperm_p_one_sided\tmedian_pct_low_prior_half\tmedian_pct_high_prior_half\n')
    for r in tests: o.write('\t'.join(map(str, r)) + '\n')
# per-sample: prior species on the same barcode in earlier runs of the same flow cell
with open('carryover_per_sample.tsv', 'w') as o:
    o.write('run\tbarcode\tdominant\tdominant_pct\t' + '\t'.join('pct_' + s for s in SP) + '\tprior_runs_same_barcode\tprior_species_reads_est\n')
    for (L, b), c in sorted(comp.items()):
        fc, t = order[L]
        prior = [(p, comp[(p, b)]['dominant'], raw.get((p, b), 0)) for p, (pfc, pt) in sorted(order.items(), key=lambda x: x[1][1])
                 if pfc == fc and pt < t and (p, b) in comp and raw.get((p, b), 0) >= 100]
        est = {s: round(load(L, b, s)) for s in SP}
        o.write(f"{L}\t{b}\t{c['dominant']}\t{c['dominant_pct']}\t" + '\t'.join(c['pct_' + s] for s in SP) + '\t' +
                ';'.join(f'{p}:{d}({n})' for p, d, n in prior) + '\t' + ';'.join(f'{s}={v}' for s, v in est.items() if v) + '\n')
