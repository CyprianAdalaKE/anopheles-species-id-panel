#!/usr/bin/env python3
"""Per-sample check table for one species folder: run / flow cell / raw barcode, verified-set membership,
raw verification, k-mer call, same-barcode carry-over flag, cross-barcode background flag, and a recommendation.
Usage: species_checks.py <folder> [<verified species name>]  -> <folder>_checks.tsv
Inputs: study_samples.tsv, carryover_flags.tsv, background_flags.tsv, verified MANIFEST, accuracy6 per_sample_calls.tsv,
slca003_verification.tsv, anoph_speciation_trace.tsv, composition tables (run_checks/ and slca020_021/)."""
import csv, re, sys, collections
folder = sys.argv[1]; sp = sys.argv[2] if len(sys.argv) > 2 else folder
t = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
st = lambda f: re.sub(r'\.(fastq|fatsq)$', '', f)
fl = {st(r['file']): r for r in t('carryover_flags.tsv') if r['folder'] == folder}
bg = {st(r['file']): r for r in t('background_flags.tsv') if r['folder'] == folder}
ver = {r['file'].replace('.fastq.gz', ''): r for r in t('/media/genomics/Sequence_Bkup/species_id_verified_fastqs/MANIFEST.tsv') if r['species'] == sp}
calls = {st(r['file']): r for r in t('/home/genomics/Desktop/species_id/gamcol_validation/accuracy6/per_sample_calls.tsv') if r['true'] == sp}
v3 = {st(r['file']): r['status'] for r in t('slca003_verification.tsv') if r['folder'] == folder}
tr = {st(r['file']): r['status'] for r in t('anoph_speciation_trace.tsv') if r['folder'] == folder}
comp = {}
for d in ('.', 'slca020_021'):
    for r in t(f'{d}/composition.tsv'): comp[(r['run'], r['barcode'])] = r
o = open(f'{folder}_checks.tsv', 'w')
cols = ['file', 'in_verified', 'run', 'flow_cell', 'raw_barcode', 'reads', 'basecall_model', 'raw_verification', 'kmer_call', 'top_pct',
        'second', 'second_pct', 'subsample_dominant', 'subsample_dominant_pct', 'carryover_flag', 'expected_carryover_pct',
        'previous_runs_same_barcode', 'background_flag', 'L1_to_pool', 'recommendation']
o.write('\t'.join(cols) + '\n'); summary = collections.Counter()
for s in t('study_samples.tsv'):
    if s['folder'] != folder: continue
    f = st(s['file']); v = ver.get(f); c = calls.get(f, {}); run = s['run']
    if 'NTC' in f: rv = 'PASS (byte-identical to raw)'
    elif run == 'SLCA011': rv = 'backup only (raw fastq deleted; fast5 exists)'
    elif tr.get(f) == 'PASS': rv = 'PASS (raw read-ID trace, Anoph_Speciation)'
    elif v3.get(f) == 'PASS': rv = 'PASS (raw on Sequence_Bkup)'
    elif v and v['check_vs_raw'] == 'PASS': rv = 'PASS (verified build)'
    else: rv = 'not traced to raw'
    b = bg.get(f, {}); co = fl.get(f, {}); cp = comp.get((run, s['barcode']), {})
    if 'NTC' in f: rec = 'NTC'
    elif not v: rec = 'not in verified set'
    elif b.get('flag') == 'BACKGROUND': rec = 'EXCLUDE - indistinguishable from empty-barcode background'
    elif co.get('flag', '').startswith('CARRY'): rec = 'EXCLUDE - reads explained by same-barcode carry-over'
    elif co.get('flag', '').startswith('high'): rec = 'KEEP - flag carry-over (>=5% expected)'
    elif run == 'SLCA011': rec = 'KEEP - flag: SLCA011 raw not verifiable, run-wide mixed background'
    elif c and c['pred'] != sp: rec = f"KEEP - misclassified ({c['pred']}); check specimen ID"
    else: rec = 'KEEP'
    summary[(run, rec.split(' - ')[0])] += 1
    o.write('\t'.join(map(str, [f, 'yes' if v else 'no', run, s['flow_cell'], s['barcode'], s['reads'], s['basecall_model'], rv,
        c.get('pred', ''), c.get('top_pct', ''), c.get('second', ''), c.get('second_pct', ''), cp.get('dominant', ''), cp.get('dominant_pct', ''),
        co.get('flag', ''), co.get('expected_carryover_pct', ''), co.get('previous_run_same_barcode(reads,dominant)', ''),
        b.get('flag', ''), b.get('L1_to_pool', ''), rec])) + '\n')
for k, n in sorted(summary.items()): print(f'{k[0]:22s} {k[1]:22s} {n}')
