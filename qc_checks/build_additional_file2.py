#!/usr/bin/env python3
"""Additional file 2: per-sample sequencing metadata and QC for all 374 sequenced species-identification specimens
(349 in the evaluation set + 25 excluded by QC), with run / flow-cell histories, carry-over rates, per-run accuracy.
Output: anopheles-panel/Additional_file_2_sample_metadata_QC.xlsx (yellow columns = specimen data to be filled in)."""
import csv, re, collections, statistics as st, os, math
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

P = '/home/genomics/Desktop/species_id/anopheles-panel/'; R = P + 'run_checks/'
VB = '/media/genomics/Sequence_Bkup/'; V = VB + 'species_id_verified_fastqs/'
OUT = P + 'Additional_file_2_sample_metadata_QC.xlsx'
tsv = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
strip = lambda f: re.sub(r'\.(fastq|fatsq)(\.gz)?$', '', f)

RUNS = {  # label: (protocol_group_id, flow cell, use no., started, run id, MinKNOW, kit, basecall model, raw data, NTC barcode, NTC reads)
 'Anoph_Speciation': ('Anoph_Speciation', 'FBA66226', 1, '2024-11-08 11:46', '0f4e7057-a418-4212-a6b7-66ffb9f09f34', '23.11.7', 'SQK-NBD114-96', 'hac@v4.2.0', 'Sequence_Bkup/SequenceData/Anoph_Speciation (fastq + pod5)', '', ''),
 'SLCA003_BatchKT': ('2025-06-30_SeqLib_SLCA003_ATSB_BatchKT', 'FBA66584', 1, '2025-06-30 16:45', 'c3b9807d-9dc8-418b-8de3-6257b3bd4bb7', '', '', 'hac@v4.3.0', 'Sequence_Bkup (fastq_pass; bc36/56 via nomadic copy)', 'barcode89', 1),
 'SLCA011': ('2025-09-10-SeqLib_SLCA011', 'FBA66584', 5, '2025-09-10 17:24', '8f4942d6-ffd3-4931-9cdb-4db2b04055e5', '25.03.9', '', 'hac@v5.0.0', 'fastq_pass deleted; demultiplexed backup + fast5 (48 GB) on Sequence_Bkup', 'barcode95', 1),
 'SLCA010_BatchST': ('2025-10-06_SeqLib_SLCA010_NMCP_BatchST', 'FBD61062', 1, '2025-10-13 09:18', '4744ed18-4739-430a-82c8-6423c70d17ab', '25.05.14', 'SQK-NBD114-96', 'hac@v5.0.0', '/var/lib/minknow/data + Sequence_Bkup', 'barcode95', 1),
 'SLCA013_BatchCOL2': ('2025-10-07_SeqLib_SLCA013_NMCP_BatchCOL2', 'FBD61062', 2, '2025-10-13 16:51', '036b812a-ad88-4742-9039-be5b44e3eb8f', '25.05.14', 'SQK-NBD114-96', 'hac@v5.0.0', '/var/lib/minknow/data + Sequence_Bkup', 'barcode94', 4),
 'SLCA015_BatchFUN': ('2025-10-14_SeqLib_SLCA015_ATSB_BatchFUN', 'FBD61062', 4, '2025-10-14 16:49', '30e151b3-a8dd-44b0-a493-e8a75ec79bb3', '25.05.14', 'SQK-NBD114-96', 'hac@v5.0.0', '/var/lib/minknow/data + Sequence_Bkup', 'barcode95', 1),
 'SLCA020': ('2026-10-19_seqLib_SLCA020', 'FBA58185', 1, '2026-10-19 15:11', '9dc3edd6-07db-4e20-b057-5b3c7d7e0886', '25.05.14', 'SQK-NBD114-96', 'hac@v5.0.0', '/var/lib/minknow/data', '', ''),
 'SLCA021': ('2026-10-20_SeqLib_SLCA021', 'FBA58542', 1, '2026-10-20 12:32', '4228b2cd-fe21-480c-8ca6-bf3ae15a7a79', '25.05.14', 'SQK-NBD114-96', 'hac@v5.0.0', '/var/lib/minknow/data', '', ''),
}
CODE = {'Anoph_Speciation': 'ASPEC', 'SLCA003_BatchKT': 'SLCA003', 'SLCA011': 'SLCA011', 'SLCA010_BatchST': 'SLCA010',
        'SLCA013_BatchCOL2': 'SLCA013', 'SLCA015_BatchFUN': 'SLCA015', 'SLCA020': 'SLCA020', 'SLCA021': 'SLCA021'}
ABBR = {'arabiensis': 'ARA', 'gambiae': 'GAM', 'funestus': 'FUN', 'coluzzii': 'COL', 'coustani': 'COU', 'stephensi': 'STE'}
FLOWCELLS = [
 ('FBA66226', 1, 'Anoph_Speciation', '2024-11-08 11:46', 'yes'), ('FBA66226', 2, 'An_funestus_08_nov_2024', '2024-11-08 17:07', 'no'),
 ('FBA66226', 3, '2025-06-24_SeqLib_SLMG001', '2025-06-24 18:11', 'no'), ('FBA66226', 4, 'dss_kmf_ksda', '2025-07-02 12:09', 'no'),
 ('FBA66584', 1, '2025-06-30_SeqLib_SLCA003_ATSB_BatchKT', '2025-06-30 16:45', 'yes'), ('FBA66584', 2, '2025-09-10_SeqLib_SLCA010 (stephensi trial)', '2025-09-10 12:32', 'no'),
 ('FBA66584', 3, '2025-09-10_SeqLib_SLC010_2 (trial, 0 reads)', '2025-09-10 13:57', 'no'), ('FBA66584', 4, '2025-09-10_SeqLib_SLCA010_3 (stephensi trial)', '2025-09-10 14:09', 'no'),
 ('FBA66584', 5, '2025-09-10-SeqLib_SLCA011', '2025-09-10 17:24', 'yes'), ('FBA66584', 6, '2025-09-11-SeqLib_SCLA012', '2025-09-11 11:15', 'no'),
 ('FBA66584', 7, '2025-10-07_SeqLib_SLCA013_NMCP_BatchCOL (same library as BatchCOL2; not used)', '2025-10-07 16:50', 'no'),
 ('FBD61062', 1, '2025-10-06_SeqLib_SLCA010_NMCP_BatchST', '2025-10-13 09:18', 'yes'), ('FBD61062', 2, '2025-10-07_SeqLib_SLCA013_NMCP_BatchCOL2', '2025-10-13 16:51', 'yes'),
 ('FBD61062', 3, '2025-10-13_SeqLib_SLCA014_ATSB_BatchINF (Plasmodium experiment)', '2025-10-14 10:10', 'Plasmodium only'), ('FBD61062', 4, '2025-10-14_SeqLib_SLCA015_ATSB_BatchFUN', '2025-10-14 16:49', 'yes'),
 ('FBD61062', 5, '2025-10-29_SeqLib_SLBM001_ATSB_BatchT', '2025-10-29 15:33', 'no'),
 ('FBA58185', 1, '2026-10-19_seqLib_SLCA020 (preceded only by flow-cell checks; logs start 2026-10-18)', '2026-10-19 15:11', 'yes'),
 ('FBA58542', 1, '2026-10-20_SeqLib_SLCA021 (preceded only by a flow-cell check)', '2026-10-20 12:32', 'yes')]

# ---- map every verified / excluded file to (run, raw barcode) ----
study = {(r['folder'].replace('_rebuilt', ''), strip(r['file'])): (r['run'], r['barcode']) for r in tsv(R + 'study_samples.tsv')}
FIX = {('arabiensis', 'combined.barcode08'): ('Anoph_Speciation', 'barcode01'), ('gambiae', 'combined.barcode53'): ('SLCA003_BatchKT', 'barcode42')}
COLSRC = {'02': ('Anoph_Speciation', 'barcode19'), '03': ('SLCA013_BatchCOL2', 'barcode10'), '04': ('SLCA013_BatchCOL2', 'barcode30'),
          '05': ('SLCA013_BatchCOL2', 'barcode52'), '06': ('SLCA013_BatchCOL2', 'barcode60'), '07': ('SLCA013_BatchCOL2', 'barcode77'),
          '08': ('SLCA013_BatchCOL2', 'barcode79'), '09': ('SLCA003_BatchKT', 'barcode27'), '10': ('SLCA003_BatchKT', 'barcode38'),
          '11': ('SLCA003_BatchKT', 'barcode39'), '13': ('SLCA011', 'barcode56'), '14': ('SLCA011', 'barcode58')}
COUSRC = {'01': ('SLCA011', 'barcode04'), '02': ('SLCA011', 'barcode13'), '03': ('SLCA011', 'barcode21')}
COLEX = {'18': 'barcode09', '20': 'barcode12', '23': 'barcode30', '24': 'barcode33', '25': 'barcode34', '27': 'barcode37', '28': 'barcode38', '29': 'barcode39', '31': 'barcode41'}
COUEX = {'10': 'barcode08', '17': 'barcode15'}
man_old = {(r['species'], r['file']): r for r in tsv(VB + 'species_id_verified_fastqs_pre-QC-final_2026-11-11/MANIFEST.tsv')}
def locate(sp, f, note=''):
    n = re.search(r'(\d+)$', f).group(1)
    if (sp, f) in FIX: return FIX[(sp, f)]
    if sp == 'arabiensis' and n in ('59', '67', '68'): return ('SLCA013_BatchCOL2', 'barcode' + n)
    if sp == 'coluzzii':
        if n in COLSRC: return COLSRC[n]
        if n in COLEX: return ('SLCA020', COLEX[n])
        return ('SLCA020', re.search(r'raw \S+ (barcode\d+)', note).group(1))
    if sp == 'coustani':
        if n in COUSRC: return COUSRC[n]
        if n in COUEX: return ('SLCA021', COUEX[n])
        return ('SLCA021', re.search(r'raw \S+ (barcode\d+)', note).group(1))
    return study[(sp, f)]

samples = []
for m in tsv(V + 'MANIFEST.tsv'):
    f = strip(m['file']); samples.append(dict(species=m['species'], file=f, reads=int(m['reads']), status='included', reason='', note=m['note'], src=m['raw_source_run(reads)']))
for r in tsv(VB + 'species_id_excluded_background_2026-11-11/EXCLUDED.tsv'):
    f = strip(r['file']); o = man_old.get((r['species'], r['file']), {})
    samples.append(dict(species=r['species'], file=f, reads=int(o.get('reads', 0)), status='excluded', reason=r['reason'], note=o.get('note', ''), src=o.get('raw_source_run(reads)', '')))
for f, n in (('combined.barcode01', 1212), ('combined.barcode02', 1021)):
    samples.append(dict(species='funestus', file=f, reads=n, status='excluded', reason='same-barcode flow-cell carry-over (reads = expected carry-over from SLCA014 BatchINF)', note='', src=''))
for s in samples: s['run'], s['bc'] = locate(s['species'], s['file'], s['note'])
assert len(samples) == 374, len(samples)

# ---- per-sample analysis results ----
raw = {(r['run'], r['barcode']): int(r['raw_reads']) for d in ('', 'slca020_021/') for r in tsv(R + d + 'raw_barcode_counts.tsv')}
comp = {(r['run'], r['barcode']): r for d in ('', 'slca020_021/') for r in tsv(R + d + 'composition.tsv')}
runs_t = {r['label']: (r['flow_cell'], r['started']) for r in tsv(R + 'runs.tsv')}
pool = {r['run']: r for r in tsv(R + 'background_runs.tsv')}
calls = {(r['true'], strip(r['file'])): r for r in tsv('/home/genomics/Desktop/species_id/gamcol_validation/accuracy6/per_sample_calls.tsv')}
nov = {(r['species'], r['file']): r for r in tsv(R + 'novelty_per_sample.tsv')}
SP6 = ['arabiensis', 'gambiae', 'coluzzii', 'funestus', 'stephensi', 'coustani']
def exp_co(run, bc, reads):
    if run not in runs_t: return ''
    fc, t = runs_t[run]; prev = sorted([(pt, p) for p, (pfc, pt) in runs_t.items() if pfc == fc and pt < t]); e = 0
    for i, (pt, p) in enumerate(reversed(prev)): e += raw.get((p, bc), 0) * (0.01 if i == 0 else 0.001)
    return round(100 * e / reads, 2) if reads else ''
def l1(run, bc):
    c = comp.get((run, bc)); p = pool.get(run)
    return round(sum(abs(float(c['pct_' + x]) - float(p['pool_pct_' + x])) for x in SP6), 1) if c and p else ''
FILL = ['collection_site', 'county', 'latitude', 'longitude', 'collection_date', 'collection_method', 'life_stage_sex', 'reference_ID_method',
        'reference_ID_by', 'dna_extraction_batch', 'biosample_accession', 'sra_run_accession', 'notes']
rows = []
for s in sorted(samples, key=lambda s: (s['species'], s['status'] != 'included', int(re.search(r'(\d+)$', s['file']).group(1)))):
    rr = RUNS[s['run']]; c = calls.get((s['species'], s['file']), {}) if s['status'] == 'included' else {}
    nv = nov.get((s['species'], s['file']), {}) if s['status'] == 'included' else {}
    rawv = ('verified against archived demultiplexed copy only (raw fastq deleted)' if s['run'] == 'SLCA011' else 'identical to raw basecalled output')
    unres = '' if not nv else ('unresolved (<95%)' if float(nv['median_best_identity']) < 95 else 'resolved')
    rows.append([f"{ABBR[s['species']]}_{CODE[s['run']]}_{s['bc'].replace('barcode', 'bc')}", f"An. {s['species']}", s['status'], s['reason'],
                 (f"species_id_verified_fastqs/{s['species']}/{s['file']}.fastq.gz" if s['status'] == 'included' else f"excluded: {s['file']}"),
                 rr[0], rr[3], rr[1], rr[2], rr[4], s['bc'], rr[7], s['reads'], rawv, exp_co(s['run'], s['bc'], s['reads']), l1(s['run'], s['bc']),
                 rr[9] or 'none identified', rr[10], c.get('pred', ''), float(c['top_pct']) if c else '', c.get('second', ''),
                 float(c['second_pct']) if c else '', ('' if not c else int(c['pred'] == s['species'])),
                 float(nv['median_best_identity']) if nv else '', unres] + [''] * len(FILL))
HDR = ['sample_id', 'reference_species', 'qc_status', 'exclusion_reason', 'file', 'run', 'run_started', 'flow_cell', 'flow_cell_use_no', 'run_id',
       'original_barcode', 'basecall_model', 'reads', 'raw_data_check', 'expected_carryover_pct', 'L1_to_run_pool', 'run_NTC_barcode', 'run_NTC_reads',
       'kmer_call', 'kmer_top_pct', 'kmer_second', 'kmer_second_pct', 'correct', 'median_best_read_identity_pct', 'read_identity_rule'] + FILL

# ---- workbook ----
wb = Workbook(); AR = 'Arial'
H = Font(name=AR, bold=True, color='FFFFFF'); HF = PatternFill('solid', fgColor='1F4E78'); Y = PatternFill('solid', fgColor='FFFF00')
N = Font(name=AR, size=10); B = Font(name=AR, size=10, bold=True); RED = PatternFill('solid', fgColor='F8CBAD')
thin = Side(style='thin', color='BFBFBF'); BOX = Border(left=thin, right=thin, top=thin, bottom=thin)
def sheet(ws, hdr, data, widths=None, fill=()):
    ws.append(hdr)
    for i, h in enumerate(hdr, 1):
        c = ws.cell(1, i); c.font = Font(name=AR, bold=True, color='000000') if h in fill else H; c.fill = Y if h in fill else HF
        c.alignment = Alignment(wrap_text=True, vertical='top'); c.border = BOX
    for r in data: ws.append(r)
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
        for c in row:
            c.font = N; c.border = BOX
            if hdr[c.column - 1] in fill: c.fill = Y
    for i, h in enumerate(hdr, 1): ws.column_dimensions[get_column_letter(i)].width = (widths or {}).get(h, max(10, min(36, len(h) + 2)))
    ws.freeze_panes = 'B2'; ws.auto_filter.ref = ws.dimensions; ws.row_dimensions[1].height = 42

ws = wb.active; ws.title = 'README'
for t, f in [
 ('Additional file 2 - per-sample sequencing metadata and quality control', B), ('', N),
 ('Samples: all 374 specimens sequenced for the species-identification evaluation; 349 passed QC (qc_status = included) and form the evaluation set.', N),
 ('QC rules (applied to all species): (1) same-barcode flow-cell carry-over - reads not exceeding the expected carry-over (1% of reads on the same barcode in the preceding run on the flow cell + 0.1% of earlier runs) -> excluded; '
  '(2) cross-barcode background - in runs with barcodes carrying no sample, reads <= 1.5x the highest such barcode AND species composition within 15 points (L1) of the read-weighted run pool -> excluded; (3) fewer than 1,000 reads -> excluded.', N),
 ('expected_carryover_pct: expected carried-over reads as % of the sample reads. L1_to_run_pool: sum of absolute differences (percentage points) between the barcode\'s species composition (Kraken2, up to 5,000 reads) and the run pool.', N),
 ('kmer_*: top two Bracken species and their read fractions (evaluation set only). median_best_read_identity_pct: median over reads of the best identity to the 60 species-specific amplicon references; read_identity_rule: unresolved if < 95%.', N),
 ('Yellow columns: specimen and accession data to be completed by the study team. Example of the expected format (not real data): collection_site = Ahero; county = Kisumu; latitude = -0.1745; longitude = 34.9180; '
  'collection_date = 2024-03-14; collection_method = CDC light trap (indoor); life_stage_sex = adult female; reference_ID_method = SINE200/IMP-PCR; reference_ID_by = initials; dna_extraction_batch = E03; biosample_accession = SAMN...; sra_run_accession = SRR...', N),
 ('sample_id = <species>_<run code>_<original barcode>; file names in the verified folder are renumbered and should not be used as identifiers in public repositories.', N), ('', N),
 ('Sheets: Samples; Runs; Flow_cell_history; Carryover_rates; Per_run_accuracy; QC_summary.', N)]:
    ws.append([t]); ws.cell(ws.max_row, 1).font = Font(name=AR, size=12 if t.startswith('Additional') else 10, bold=f is B); ws.cell(ws.max_row, 1).alignment = Alignment(wrap_text=True, vertical='top')
ws.column_dimensions['A'].width = 170

ws = wb.create_sheet('Samples')
sheet(ws, HDR, rows, widths={'sample_id': 20, 'exclusion_reason': 40, 'file': 50, 'run': 40, 'run_id': 38, 'raw_data_check': 40, 'notes': 30}, fill=FILL)
nS = ws.max_row
for i in range(2, nS + 1):
    if ws.cell(i, 3).value == 'excluded':
        for c in ws[i][:len(HDR) - len(FILL)]: c.fill = RED
    for col, fmt in (('reads', '#,##0'), ('kmer_top_pct', '0.0'), ('kmer_second_pct', '0.0'), ('expected_carryover_pct', '0.00'), ('median_best_read_identity_pct', '0.00')):
        ws.cell(i, HDR.index(col) + 1).number_format = fmt

ws = wb.create_sheet('Runs')
hdr = ['run', 'run_started', 'flow_cell', 'flow_cell_use_no', 'run_id', 'minknow_version', 'kit', 'basecall_model', 'raw_data', 'NTC_barcode', 'NTC_reads',
       'samples_sequenced', 'samples_included', 'reads_included', 'flow_cell_washed_before_run', 'library_prep_date', 'operator']
data = []
for i, (k, r) in enumerate(RUNS.items(), 2):
    data.append([r[0], r[3], r[1], r[2], r[4], r[5], r[6], r[7], r[8], r[9] or 'none identified', r[10],
                 f'=COUNTIFS(Samples!$F$2:$F${nS},A{i})', f'=COUNTIFS(Samples!$F$2:$F${nS},A{i},Samples!$C$2:$C${nS},"included")',
                 f'=SUMIFS(Samples!$M$2:$M${nS},Samples!$F$2:$F${nS},A{i},Samples!$C$2:$C${nS},"included")', '', '', ''])
n = len(RUNS) + 1
data.append(['TOTAL'] + [''] * 10 + [f'=SUM(L2:L{n})', f'=SUM(M2:M{n})', f'=SUM(N2:N{n})', '', '', ''])
sheet(ws, hdr, data, widths={'run': 40, 'run_id': 38, 'raw_data': 55}, fill=['flow_cell_washed_before_run', 'library_prep_date', 'operator', 'kit', 'minknow_version'])
for c in ws[ws.max_row]: c.font = B
for i in range(2, ws.max_row + 1): ws.cell(i, 14).number_format = '#,##0'

ws = wb.create_sheet('Flow_cell_history')
sheet(ws, ['flow_cell', 'use_no', 'run', 'started', 'species_identification_samples', 'washed_before_this_run'], [list(r) + [''] for r in FLOWCELLS],
      widths={'run': 70, 'started': 18}, fill=['washed_before_this_run'])

ws = wb.create_sheet('Carryover_rates')
sheet(ws, ['flow_cell', 'transition', 'barcodes', 'n', 'measure', 'median_pct', 'range_or_p90', 'note'], [
 ['FBD61062', 'SLCA014 BatchINF -> SLCA015 BatchFUN', 'BatchFUN barcodes 01-48 (used in INF, not FUN); INF > 5,000 reads', 30, 'FUN reads / INF reads, same barcode', 0.99, '0.76-1.14', 'Plasmodium % INF vs FUN, same barcode: Spearman rho 0.94 (n = 48)'],
 ['FBD61062', 'SLCA010 BatchST -> SLCA013 BatchCOL2', 'COL2 barcodes not stephensi-dominated; ST > 5,000 reads', 73, 'stephensi reads in COL2 / ST reads, same barcode', 0.60, 'p90 1.10', 'barcode-specific Spearman rho 0.49, permutation p = 0.0005'],
 ['FBD61062', 'SLCA013 BatchCOL2 -> SLCA014 BatchINF', 'INF barcodes 49-80 (unused in INF); COL2 > 5,000 reads', 31, 'INF reads / COL2 reads, same barcode', 0.43, '', ''],
 ['FBA58185', 'cross-barcode background (SLCA020)', 'barcodes 79, 80 (no sample assigned)', 2, 'reads per empty barcode', '', '1,365-2,143', 'composition equals the read-weighted run pool'],
 ['FBA58542', 'cross-barcode background (SLCA021)', 'barcode 80 (no sample assigned)', 1, 'reads per empty barcode', '', '1,998', 'composition equals the run pool']],
 widths={'transition': 38, 'barcodes': 55, 'measure': 42, 'note': 60})

ws = wb.create_sheet('Per_run_accuracy')
g = collections.defaultdict(list)
for r in rows:
    if r[2] == 'included': g[(r[1], r[5])].append(r)
data = []
for i, ((sp, run), v) in enumerate(sorted(g.items()), 2):
    data.append([sp, run, len(v), sum(x[22] for x in v), f'=D{i}/C{i}', round(st.median(x[19] for x in v), 1), sum(x[24] == 'unresolved (<95%)' for x in v)])
sheet(ws, ['reference_species', 'run', 'n', 'correct', 'accuracy', 'median_top_pct', 'unresolved_by_identity_rule'], data, widths={'run': 42})
for i in range(2, ws.max_row + 1): ws.cell(i, 5).number_format = '0.0%'

ws = wb.create_sheet('QC_summary')
q = collections.Counter((r[1], r[2], r[3].split(' (')[0]) for r in rows)
sheet(ws, ['reference_species', 'qc_status', 'reason', 'n'], [[a, b, c, n] for (a, b, c), n in sorted(q.items())], widths={'reason': 60})
wb.save(OUT); print(OUT, len(rows), 'rows')
