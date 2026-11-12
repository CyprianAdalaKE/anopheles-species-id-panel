#!/usr/bin/env python3
"""Per-run accuracy of the six-species k-mer/LCA call on the verified set (accuracy6/per_sample_calls.tsv).
Run per verified file: arabiensis/funestus/gambiae/stephensi verified files = final-folder files of the same name
(read counts checked); coluzzii/coustani from the verified MANIFEST source runs ('species_id' = SLCA011).
Also flags verified files that the carry-over check marks as carry-over. Writes per_run_accuracy.tsv, verified_runs.tsv."""
import csv, re, collections
V = '/media/genomics/Sequence_Bkup/species_id_verified_fastqs/MANIFEST.tsv'
calls = {(r['true'], re.sub(r'\.fastq$', '', r['file'])): r for r in csv.DictReader(open('/home/genomics/Desktop/species_id/gamcol_validation/accuracy6/per_sample_calls.tsv'), delimiter='\t')}
study = {}
for r in csv.DictReader(open('study_samples.tsv'), delimiter='\t'):
    study[(r['folder'], re.sub(r'\.(fastq|fatsq)$', '', r['file']))] = r
flags = {(r['folder'], re.sub(r'\.(fastq|fatsq)$', '', r['file'])): r for r in csv.DictReader(open('carryover_flags.tsv'), delimiter='\t')}
SHORT = {'2025-06-30_SeqLib_SLCA003_ATSB_BatchKT': 'SLCA003_BatchKT', '2025-10-07_SeqLib_SLCA013_NMCP_BatchCOL2': 'SLCA013_BatchCOL2',
         '2026-10-19_seqLib_SLCA020': 'SLCA020', '2026-10-20_SeqLib_SLCA021': 'SLCA021', 'species_id': 'SLCA011',
         'anoph_speciation_fastq_files': 'Anoph_Speciation', 'Anoph_Speciation': 'Anoph_Speciation'}
rows = []
for m in csv.DictReader(open(V), delimiter='\t'):
    sp = m['species']; f = m['file'].replace('.fastq.gz', '')
    # stephensi barcode16 was named combined.barcode15 when the calls were made (renamed 2026-11-11)
    c = calls.get((sp, 'combined.barcode15' if (sp, f) == ('stephensi', 'combined.barcode16') else f))
    if not c: continue
    flag = ''
    if sp in ('coluzzii', 'coustani'):
        srcs = [s.split('(')[0] for s in m['raw_source_run(reads)'].split(';')]
        run = SHORT.get(srcs[0], srcs[0]); bc = ''
        mm = re.search(r'raw \S+ (barcode\d+)', m['note'])
        if mm: bc = mm.group(1)
        model = ''
    else:
        s = study[(sp, f)]; assert int(s['reads']) == int(m['reads']), (sp, f)
        run, bc = s['run'], s['barcode']; flag = flags[(sp, f)]['flag']
    rows.append(dict(species=sp, verified_file=m['file'], reads=m['reads'], run=run, raw_barcode=bc,
                     pred=c['pred'], top_pct=c['top_pct'], second=c['second'], second_pct=c['second_pct'],
                     correct=int(c['pred'] == sp), carryover_flag=flag))
with open('verified_runs.tsv', 'w') as o:
    w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
g = collections.defaultdict(list)
for r in rows: g[(r['species'], r['run'])].append(r)
with open('per_run_accuracy.tsv', 'w') as o:
    o.write('species\trun\tn\tcorrect\taccuracy_pct\tmedian_top_pct\tn_excl_carryover\tcorrect_excl_carryover\n')
    for (sp, run), v in sorted(g.items()):
        ok = sum(r['correct'] for r in v); tp = sorted(float(r['top_pct']) for r in v)
        ve = [r for r in v if not r['carryover_flag'].startswith('CARRY')]
        o.write(f"{sp}\t{run}\t{len(v)}\t{ok}\t{100*ok/len(v):.1f}\t{tp[len(tp)//2]:.1f}\t{len(ve)}\t{sum(r['correct'] for r in ve)}\n")
