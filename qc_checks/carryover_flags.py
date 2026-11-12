#!/usr/bin/env python3
"""Flag study samples by expected carry-over from earlier runs on the same flow cell.
Rate estimated from SLCA014_BatchINF -> SLCA015_BatchFUN (unused FUN barcodes 01-48 received median 0.99%,
range 0.76-1.14%, of INF's reads on the same barcode). Expected carry-over reads = 1% x reads on the same
barcode in the immediately preceding run on that flow cell (+0.1% x older runs, conservative).
expected_carryover_pct = expected reads / sample reads. Writes carryover_flags.tsv."""
import csv
runs = [r for r in csv.DictReader(open('runs.tsv'), delimiter='\t')]
order = {r['label']: (r['flow_cell'], r['started']) for r in runs}
raw = {(r['run'], r['barcode']): int(r['raw_reads']) for r in csv.DictReader(open('raw_barcode_counts.tsv'), delimiter='\t')}
comp = {(r['run'], r['barcode']): r for r in csv.DictReader(open('composition.tsv'), delimiter='\t')}
out = open('carryover_flags.tsv', 'w')
out.write('folder\tfile\trun\tflow_cell\tbarcode\treads\tflow_cell_use\tprevious_run_same_barcode(reads,dominant)\texpected_carryover_reads\texpected_carryover_pct\tsubsample_dominant\tsubsample_Plasmodium_pct\tflag\n')
for s in csv.DictReader(open('study_samples.tsv'), delimiter='\t'):
    run, bc = s['run'], s['barcode']
    if run not in order:
        out.write(f"{s['folder']}\t{s['file']}\t{run}\t{s['flow_cell']}\t{bc}\t{s['reads']}\tunknown\t\t\t\t\t\tno flow-cell history on this machine\n"); continue
    fc, t = order[run]
    prev = sorted([(pt, p) for p, (pfc, pt) in order.items() if pfc == fc and pt < t])
    use = len(prev) + 1
    exp = 0.0; desc = []
    for i, (pt, p) in enumerate(reversed(prev)):
        n = raw.get((p, bc), 0); exp += n * (0.01 if i == 0 else 0.001)
        if n >= 100: desc.append(f"{p}:{n}({comp[(p, bc)]['dominant'] if (p, bc) in comp else '-'})")
    reads = int(s['reads']); pct = 100 * exp / reads if reads else 0
    c = comp.get((run, bc), {})
    flag = ('CARRY-OVER: sample reads ~ expected carry-over' if pct >= 50 else
            'high (>=5%)' if pct >= 5 else 'moderate (1-5%)' if pct >= 1 else 'low (<1%)')
    if 'NTC' in s['file']: flag = 'NTC'
    out.write(f"{s['folder']}\t{s['file']}\t{run}\t{fc}\t{bc}\t{reads}\t{use}\t{';'.join(desc)}\t{exp:.0f}\t{pct:.2f}\t{c.get('dominant','')}\t{c.get('pct_Plasmodium','')}\t{flag}\n")
