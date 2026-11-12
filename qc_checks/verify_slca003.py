#!/usr/bin/env python3
"""Verify final-folder SLCA003 BatchKT files against raw fastq_pass on Sequence_Bkup (read ID + sequence, per barcode)."""
import csv, gzip, glob, hashlib, sys
RAW = '/media/genomics/Sequence_Bkup/2025-06-30_SeqLib_SLCA003_ATSB_BatchKT/no_sample_id/20250630_1645_MN30463_FBA66584_c3b9807d/fastq_pass'
F = '/home/genomics/Desktop/species_id/anopheles-panel/species_id_final_fastqs'
h = lambda s: hashlib.md5(s.strip().encode()).digest()[:8]
def reads(fh):
    for i, l in enumerate(fh):
        if i % 4 == 0: rid = l[1:].split()[0]
        elif i % 4 == 1: yield rid, h(l)
out = open('slca003_verification.tsv', 'w'); out.write('folder\tfile\traw_barcode\treads\traw_reads\tidentical\tstatus\n')
for s in csv.DictReader(open('study_samples.tsv'), delimiter='\t'):
    if s['run'] != 'SLCA003_BatchKT': continue
    rawset = {}
    for f in glob.glob(f"{RAW}/{s['barcode']}/FBA66584_pass_{s['barcode']}_*.fastq.gz"):
        with gzip.open(f, 'rt') as fh: rawset.update(reads(fh))
    with open(f"{F}/{s['folder']}/{s['file']}") as fh: fin = dict(reads(fh))
    same = sum(rawset.get(k) == v for k, v in fin.items())
    st = 'PASS' if same == len(fin) == len(rawset) else ('PASS (subset of raw)' if same == len(fin) else 'FAIL')
    out.write(f"{s['folder']}\t{s['file']}\t{s['barcode']}\t{len(fin)}\t{len(rawset)}\t{same}\t{st}\n"); print(s['folder'], s['file'], st, file=sys.stderr)
