#!/usr/bin/env python3
"""Subsample up to N reads per (run, barcode) for a carry-over composition check.
Reads are taken evenly across each barcode's time-ordered fastq chunks (ceil(N/nfiles) per file),
so early and late parts of the run are both represented. Writes one fastq with headers
@<run>|<barcode>|<i> and a counts table (total reads per run/barcode in raw)."""
import sys, os, gzip, glob, re, math, collections
N = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
runs = [l.rstrip('\n').split('\t') for l in open(sys.argv[1])][1:]
out = open('subsample.fastq', 'w'); cnt = open('raw_barcode_counts.tsv', 'w')
cnt.write('run\tflow_cell\tbarcode\traw_reads\tsampled\n')
def records(fh):
    while True:
        h = fh.readline()
        if not h: return
        s = fh.readline(); fh.readline(); q = fh.readline()
        yield h, s, q
for label, fc, started, fp in runs:
    if fp.startswith('BACKUP:'):          # SLCA011: raw gone; backup is demultiplexed per barcode
        groups = collections.defaultdict(list)
        for f in sorted(glob.glob(fp[7:] + '/combined.barcode*.fastq')):
            with open(f) as fh:
                h = fh.readline(); m = re.search(r'barcode=(\S+)', h)
            groups[m.group(1) if m else os.path.basename(f)].append(f)
    else:
        groups = {}
        for d in sorted(glob.glob(fp + '/barcode*')) + [fp + '/unclassified']:
            bc = os.path.basename(d)
            fs = sorted(glob.glob(f'{d}/{fc}_pass_{bc}_*.fastq.gz'),
                        key=lambda x: int(re.search(r'_(\d+)\.fastq\.gz$', x).group(1)))
            if fs: groups[bc] = fs
    for bc, fs in sorted(groups.items()):
        per = math.ceil(N / len(fs)); total = 0; taken = 0; i = 0
        for f in fs:
            op = gzip.open(f, 'rt') if f.endswith('.gz') else open(f)
            with op as fh:
                k = 0
                for h, s, q in records(fh):
                    total += 1
                    if k < per and taken < N and bc != 'unclassified':
                        i += 1; k += 1; taken += 1
                        out.write(f'@{label}|{bc}|{i}\n{s}+\n{q}')
        cnt.write(f'{label}\t{fc}\t{bc}\t{total}\t{taken}\n')
        print(label, bc, total, taken, file=sys.stderr, flush=True)
