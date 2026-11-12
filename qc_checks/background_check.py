#!/usr/bin/env python3
"""Cross-barcode background check. For each run, pooled composition = read-weighted mean over all barcodes.
A sample whose composition matches the pool (L1 distance over the six Anopheles species < 15 points) and whose
read count is within the run's background level (<= 3,500 reads; SLCA020 unused barcodes 79/80 carry 1,365-2,143
reads with pool composition) cannot be distinguished from an empty barcode: flag BACKGROUND.
Writes background_flags.tsv (all study samples) and background_runs.tsv."""
import csv, collections
SP = ['arabiensis', 'gambiae', 'coluzzii', 'funestus', 'stephensi', 'coustani']
comp, raw = {}, {}
for d in ('.', 'slca020_021'):
    for r in csv.DictReader(open(f'{d}/composition.tsv'), delimiter='\t'): comp[(r['run'], r['barcode'])] = r
    for r in csv.DictReader(open(f'{d}/raw_barcode_counts.tsv'), delimiter='\t'): raw[(r['run'], r['barcode'])] = int(r['raw_reads'])
runs = sorted({r for r, b in comp})
pool = {}
with open('background_runs.tsv', 'w') as o:
    o.write('run\ttotal_reads\t' + '\t'.join('pool_pct_' + s for s in SP) + '\tn_barcodes_pool_like_le3500\n')
    for run in runs:
        bcs = [b for (r, b) in comp if r == run and b.startswith('barcode')]
        tot = sum(raw.get((run, b), 0) for b in bcs)
        if not tot: continue
        pool[run] = {s: sum(raw[(run, b)] * float(comp[(run, b)]['pct_' + s]) for b in bcs) / tot for s in SP}
        n = sum(1 for b in bcs if 300 <= raw[(run, b)] <= 3500 and sum(abs(float(comp[(run, b)]['pct_' + s]) - pool[run][s]) for s in SP) < 15)
        o.write(f'{run}\t{tot}\t' + '\t'.join(f'{pool[run][s]:.1f}' for s in SP) + f'\t{n}\n')
with open('background_flags.tsv', 'w') as o:
    o.write('folder\tfile\trun\tbarcode\treads\tL1_to_pool\tflag\n')
    for s in csv.DictReader(open('study_samples.tsv'), delimiter='\t'):
        k = (s['run'], s['barcode'])
        if k not in comp or s['run'] not in pool or 'NTC' in s['file']: continue
        l1 = sum(abs(float(comp[k]['pct_' + x]) - pool[s['run']][x]) for x in SP); n = int(s['reads'])
        flag = 'BACKGROUND' if n <= 3500 and l1 < 15 else ('pool-like composition' if l1 < 15 else '')
        o.write(f"{s['folder']}\t{s['file']}\t{s['run']}\t{s['barcode']}\t{n}\t{l1:.1f}\t{flag}\n")
