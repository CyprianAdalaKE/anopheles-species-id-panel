#!/usr/bin/env python3
"""Trace final-folder files whose reads come from Anoph_Speciation (or carry no header metadata)
back to raw fastq_pass on Sequence_Bkup/SequenceData (FBA66226 runs). Index = read_id -> (run, barcode, md5(seq)).
Output: anoph_speciation_trace.tsv (one row per final file)."""
import gzip, glob, os, re, hashlib, csv, collections, sys
RAW = {'Anoph_Speciation': '/media/genomics/Sequence_Bkup/SequenceData/Anoph_Speciation/no_sample/20241108_1144_MN30463_FBA66226_0f4e7057/fastq_pass',
       'An_funestus_08_nov_2024': '/media/genomics/Sequence_Bkup/SequenceData/An_funestus_08_nov_2024/no_sample/20241108_1702_MN30463_FBA66226_ced627b6/fastq_pass'}
FINAL = '/home/genomics/Desktop/species_id/anopheles-panel/species_id_final_fastqs'
AUDIT = '/home/genomics/Desktop/species_id/anopheles-panel/species_id_final_fastqs_audit_2026-10-26.tsv'
h = lambda s: hashlib.md5(s.encode()).digest()[:8]
idx = {}
for run, fp in RAW.items():
    for f in sorted(glob.glob(fp + '/*/*.fastq.gz')):
        bc = os.path.basename(os.path.dirname(f))
        with gzip.open(f, 'rt') as fh:
            for i, l in enumerate(fh):
                if i % 4 == 0: rid = l[1:].split()[0]
                elif i % 4 == 1: idx.setdefault(rid, (run, bc, h(l.strip())))
    print(run, 'indexed', len(idx), file=sys.stderr, flush=True)
# candidate final files: audit run == Anoph_Speciation or NA (no header); current folders only
cand = [r for r in csv.DictReader(open(AUDIT), delimiter='\t') if r['run'] in ('Anoph_Speciation', 'NA', '')]
folders = {'arabiensis','funestus','gambiae','stephensi','coluzzii_rebuilt','coustani_rebuilt'}
seen = set(); rows = []
for folder in sorted(folders):
    for f in sorted(os.listdir(f'{FINAL}/{folder}')):
        p = f'{FINAL}/{folder}/{f}'
        if f.endswith('.tsv') or 'NTC' in f: continue
        with open(p) as fh: first = fh.readline()
        if 'protocol_group_id=' in first and 'Anoph_Speciation' not in first: continue
        n = hit = seqok = 0; src = collections.Counter()
        with open(p) as fh:
            for i, l in enumerate(fh):
                if i % 4 == 0: rid = l[1:].split()[0]
                elif i % 4 == 1:
                    n += 1; e = idx.get(rid)
                    if e:
                        hit += 1; src[(e[0], e[1])] += 1
                        if e[2] == h(l.strip()): seqok += 1
        top = src.most_common(3)
        rows.append([folder, f, n, hit, seqok, ';'.join(f'{r}/{b}({c})' for (r, b), c in top), len(src),
                     'PASS' if n and seqok == n and len(src) == 1 else ('PARTIAL' if hit else 'NO_MATCH')])
        print(*rows[-1], file=sys.stderr, flush=True)
with open('anoph_speciation_trace.tsv', 'w') as o:
    o.write('folder\tfile\treads\tread_ids_in_raw\tseq_identical\traw_source(top3)\tn_sources\tstatus\n')
    for r in rows: o.write('\t'.join(map(str, r)) + '\n')
