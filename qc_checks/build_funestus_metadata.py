#!/usr/bin/env python3
"""Build the funestus sequencing-metadata workbook (draft) from the run checks of 2026-11-11.
Inputs (all in run_checks/ or anopheles-panel/): run_summary_funestus.tsv, study_samples.tsv, carryover_flags.tsv,
verified_runs.tsv, composition.tsv, raw_barcode_counts.tsv, carryover_tests.tsv, per_run_accuracy.tsv,
species_id_final_fastqs_audit_2026-10-26.tsv. Output: ../funestus_metadata_draft_2026-11-11.xlsx"""
import csv, re, collections
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.comments import Comment

P = '/home/genomics/Desktop/species_id/anopheles-panel/'
R = P + 'run_checks/'
OUT = P + 'funestus_metadata_draft_2026-11-11.xlsx'
tsv = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
strip = lambda f: re.sub(r'\.(fastq|fatsq)$', '', f)

# ---------------------------------------------------------------- run-level facts (from MinKNOW final_summary/report files)
RUNS = [  # label, protocol_group_id, flow cell, use#, started, run_id, minknow, kit, model, raw data location, NTC
    ('Anoph_Speciation', 'Anoph_Speciation', 'FBA66226', 1, '2024-11-08 11:46', '0f4e7057-a418-4212-a6b7-66ffb9f09f34',
     '23.11.7', 'SQK-NBD114-96', 'dna_r10.4.1_e8.2_400bps_hac@v4.2.0',
     'Sequence_Bkup/SequenceData/Anoph_Speciation (fastq_pass + pod5_pass)', ''),
    ('SLCA003_BatchKT', '2025-06-30_SeqLib_SLCA003_ATSB_BatchKT', 'FBA66584', 1, '2025-06-30 16:45', 'c3b9807d-9dc8-418b-8de3-6257b3bd4bb7',
     '', '', 'dna_r10.4.1_e8.2_400bps_hac@v4.3.0',
     'Sequence_Bkup/2025-06-30_SeqLib_SLCA003_ATSB_BatchKT (fastq_pass only)', 'barcode89'),
    ('SLCA011', '2025-09-10-SeqLib_SLCA011', 'FBA66584', 5, '2025-09-10 17:24', '8f4942d6-ffd3-4931-9cdb-4db2b04055e5',
     '25.03.9', '', 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0',
     'fastq_pass deleted; demultiplexed backup Sequence_Bkup/species_id/coustani_combined; raw signal (fast5, 48 GB) in Sequence_Bkup/2025-09-10-SeqLib_SLCA011',
     'barcode95'),
    ('SLCA013_BatchCOL2', '2025-10-07_SeqLib_SLCA013_NMCP_BatchCOL2', 'FBD61062', 2, '2025-10-13 16:51', '036b812a-ad88-4742-9039-be5b44e3eb8f',
     '25.05.14', 'SQK-NBD114-96', 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0',
     '/var/lib/minknow/data + Sequence_Bkup (fastq_pass)', 'barcode94'),
    ('SLCA015_BatchFUN', '2025-10-14_SeqLib_SLCA015_ATSB_BatchFUN', 'FBD61062', 4, '2025-10-14 16:49', '30e151b3-a8dd-44b0-a493-e8a75ec79bb3',
     '25.05.14', 'SQK-NBD114-96', 'dna_r10.4.1_e8.2_400bps_hac@v5.0.0',
     '/var/lib/minknow/data + Sequence_Bkup (fastq_pass)', 'barcode95'),
]
RUN = {r[0]: r for r in RUNS}
CODE = {'Anoph_Speciation': 'ASPEC', 'SLCA003_BatchKT': 'SLCA003', 'SLCA011': 'SLCA011',
        'SLCA013_BatchCOL2': 'SLCA013C2', 'SLCA015_BatchFUN': 'SLCA015'}
FLOWCELL_HISTORY = [  # flow cell, use#, run, started, reads/notes, study run?
    ('FBA66226', 1, 'Anoph_Speciation', '2024-11-08 11:46', 'funestus 17, arabiensis 45, coluzzii 1 in study', 'yes'),
    ('FBA66226', 2, 'An_funestus_08_nov_2024', '2024-11-08 17:07', 'not in study', 'no'),
    ('FBA66226', 3, '2025-06-24_SeqLib_SLMG001', '2025-06-24 18:11', 'not in study', 'no'),
    ('FBA66226', 4, 'dss_kmf_ksda', '2025-07-02 12:09', '2 fastq files; not in study', 'no'),
    ('FBA66584', 1, 'SLCA003_BatchKT', '2025-06-30 16:45', 'no earlier run on this flow cell found on this machine', 'yes'),
    ('FBA66584', 2, '2025-09-10_SeqLib_SLCA010 (trial)', '2025-09-10 12:32', 'stephensi library; ~100-1,400 reads/barcode', 'no'),
    ('FBA66584', 3, '2025-09-10_SeqLib_SLC010_2 (trial)', '2025-09-10 13:57', '0 reads', 'no'),
    ('FBA66584', 4, '2025-09-10_SeqLib_SLCA010_3 (trial)', '2025-09-10 14:09', 'stephensi library; ~100-1,400 reads/barcode', 'no'),
    ('FBA66584', 5, 'SLCA011', '2025-09-10 17:24', 'same day as the stephensi trials', 'yes'),
    ('FBA66584', 6, '2025-09-11-SeqLib_SCLA012', '2025-09-11 11:15', 'not in study', 'no'),
    ('FBA66584', 7, 'SLCA013_BatchCOL', '2025-10-07 16:50', 'same library as BatchCOL2', 'yes'),
    ('FBD61062', 1, 'SLCA010_BatchST', '2025-10-13 09:18', 'stephensi bc01-80', 'yes'),
    ('FBD61062', 2, 'SLCA013_BatchCOL2', '2025-10-13 16:51', 'gambiae/coluzzii/funestus bc01-80', 'yes'),
    ('FBD61062', 3, 'SLCA014_BatchINF', '2025-10-14 10:10', 'Plasmodium-infected mosquitoes bc01-48', 'yes (Plasmodium)'),
    ('FBD61062', 4, 'SLCA015_BatchFUN', '2025-10-14 16:49', 'funestus bc49-76 (+ carry-over on bc01-48)', 'yes'),
    ('FBD61062', 5, '2025-10-29_SeqLib_SLBM001_ATSB_BatchT', '2025-10-29 15:33', 'TEP1 test; not in study', 'no'),
]

# ---------------------------------------------------------------- per-sample data
summ = {r['file']: r for r in tsv(P + 'run_summary_funestus.tsv')}
study = {(r['folder'], strip(r['file'])): r for r in tsv(R + 'study_samples.tsv')}
flags = {(r['folder'], strip(r['file'])): r for r in tsv(R + 'carryover_flags.tsv')}
ver = {strip(r['verified_file'].replace('.fastq.gz', '')): r for r in tsv(R + 'verified_runs.tsv') if r['species'] == 'funestus'}
comp = {(r['run'], r['barcode']): r for r in tsv(R + 'composition.tsv')}
raw = {(r['run'], r['barcode']): int(r['raw_reads']) for r in tsv(R + 'raw_barcode_counts.tsv')}
audit = {r['file']: r for r in tsv(P + 'species_id_final_fastqs_audit_2026-10-26.tsv') if r['folder'] == 'funestus'}
v3 = {r['file']: r['status'] for r in tsv(R + 'slca003_verification.tsv') if r['folder'] == 'funestus'}
ntc_reads = {('SLCA015_BatchFUN', 'barcode95'): 1, ('SLCA013_BatchCOL2', 'barcode94'): 4,
             ('SLCA003_BatchKT', 'barcode89'): 1, ('SLCA011', 'barcode95'): 1}
ntc_where = {'SLCA015_BatchFUN': 'funestus/SLCA015_BatchFUN_barcode95_NTC.fastq',
             'SLCA013_BatchCOL2': 'gambiae/SLCA013_BatchCOL2_barcode94_NTC.fastq',
             'SLCA003_BatchKT': 'arabiensis/SLCA003_BatchKT_barcode89_NTC.fastq',
             'SLCA011': 'coustani_rebuilt/SLCA011_barcode95_NTC.fastq'}

samples = []
for f in sorted(summ, key=lambda x: (('NTC' in x), x)):
    s = study[('funestus', strip(f))]; run = s['run']; bc = s['barcode']; rr = RUN[run]
    fl = flags[('funestus', strip(f))]; c = comp.get((run, bc), {}); is_ntc = 'NTC' in f
    v = ver.get(strip(f), {})
    if run == 'Anoph_Speciation':
        rawv = 'PASS - all read IDs and sequences identical to raw (read-ID trace, 2026-11-11)'
    elif run == 'SLCA011':
        rawv = 'not verifiable - fastq_pass deleted (raw signal fast5 still on backup)'
    elif is_ntc:
        rawv = 'PASS - byte-identical to raw fastq_pass'
    elif run == 'SLCA003_BatchKT':
        rawv = 'PASS - read IDs and sequences identical to raw on Sequence_Bkup (2026-11-11)' if v3.get(f) == 'PASS' else v3.get(f, '')
    else:
        rawv = 'PASS - exact match to raw fastq_pass (audit 2026-10-26)' if audit.get(f, {}).get('raw_check') == 'exact_match' else audit.get(f, {}).get('raw_check', '')
    if is_ntc:
        action = 'NTC - report reads/classification; not a sample'
    elif fl['flag'].startswith('CARRY'):
        action = ('EXCLUDE - reads are carry-over of SLCA014 BatchINF library on the same barcode '
                  f"(~1% of INF {bc}: {raw.get(('SLCA014_BatchINF', bc), 0):,} reads; {c.get('pct_Plasmodium')}% Plasmodium); no funestus specimen DNA")
    elif run == 'SLCA011':
        action = 'KEEP - flag: raw not verifiable; whole run shows high mixed background (median top call 28%)'
    else:
        action = 'KEEP'
    sid = f"NTC_{CODE[run]}_{bc.replace('barcode', 'bc')}" if is_ntc else f"FUN_{CODE[run]}_{bc.replace('barcode', 'bc')}"
    samples.append(dict(
        sample_id=sid, species_label='An. funestus' if not is_ntc else 'NTC (no template)',
        final_file=f'species_id_final_fastqs/funestus/{f}',
        verified_file=(f'species_id_verified_fastqs/funestus/{f}.fastq.gz' if not is_ntc else ''),
        run=rr[1], run_started=rr[4], flow_cell=rr[2], flow_cell_use=rr[3], run_id=rr[5], raw_barcode=bc,
        kit=rr[7], minknow=rr[6], model=rr[8], reads=int(s['reads']), raw_verification=rawv,
        ntc_barcode=rr[10], ntc_reads=ntc_reads.get((run, rr[10]), ''), ntc_file=ntc_where.get(run, ''),
        call=v.get('pred', ''), top_pct=float(v['top_pct']) if v else '', second=v.get('second', ''),
        second_pct=float(v['second_pct']) if v else '',
        sub_fun=float(c['pct_funestus']) if c else '', sub_pl=float(c['pct_Plasmodium']) if c else '',
        prev_same_bc=fl['previous_run_same_barcode(reads,dominant)'],
        exp_co=float(fl['expected_carryover_pct']) if fl['expected_carryover_pct'] else '', co_flag=fl['flag'], action=action))

# ---------------------------------------------------------------- workbook
wb = Workbook()
AR = 'Arial'
H = Font(name=AR, bold=True, color='FFFFFF'); HF = PatternFill('solid', fgColor='1F4E78')
FILL = PatternFill('solid', fgColor='FFFF00'); RED = PatternFill('solid', fgColor='F8CBAD'); AMB = PatternFill('solid', fgColor='FCE4D6')
GREY = PatternFill('solid', fgColor='EDEDED'); N = Font(name=AR, size=10); B = Font(name=AR, size=10, bold=True)
thin = Side(style='thin', color='BFBFBF'); BOX = Border(left=thin, right=thin, top=thin, bottom=thin)

def sheet(ws, headers, rows, widths=None, fill_cols=(), comments=None):
    ws.append(headers)
    for i, h in enumerate(headers, 1):
        c = ws.cell(1, i); c.font = H; c.fill = FILL if h in fill_cols else HF
        if h in fill_cols: c.font = Font(name=AR, bold=True, color='000000')
        c.alignment = Alignment(wrap_text=True, vertical='top'); c.border = BOX
        if comments and h in comments: c.comment = Comment(comments[h], 'run checks')
    for r in rows:
        ws.append(r)
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
        for c in row:
            c.font = N; c.border = BOX; c.alignment = Alignment(vertical='top', wrap_text=False)
            if headers[c.column - 1] in fill_cols: c.fill = FILL
    for i, h in enumerate(headers, 1):
        ws.column_dimensions[get_column_letter(i)].width = (widths or {}).get(h, max(10, min(40, len(h) + 2)))
    ws.freeze_panes = 'B2'; ws.auto_filter.ref = ws.dimensions; ws.row_dimensions[1].height = 45

# README -------------------------------------------------------------------------------
ws = wb.active; ws.title = 'README'
lines = [
    ('An. funestus sequencing metadata - DRAFT (2026-11-11)', B),
    ('One row per specimen in Funestus_samples. Sequencing columns are filled from read headers and MinKNOW run files; yellow columns are for you to fill (specimen/collection data).', N),
    ('', N),
    ('Sheets', B),
    ('Funestus_samples - 79 funestus files + 1 NTC: stable sample ID, run, flow cell, original barcode, verification, classification, carry-over check, recommended action.', N),
    ('Runs - one row per sequencing run contributing funestus samples (counts are formulas over Funestus_samples).', N),
    ('Flow_cell_history - every run found on the flow cells used, in order (reuse matters for carry-over).', N),
    ('Carryover_rates - observed carry-over between consecutive runs on the same flow cell.', N),
    ('Carryover_tests - barcode-specific carry-over test (Spearman) for every later run x species.', N),
    ('Per_run_accuracy - six-species call accuracy of the verified set broken down by run (all species).', N),
    ('', N),
    ('Legend', B),
    ('Yellow header/cells = to be filled in by the study team. Red row = recommended EXCLUDE. Orange row = keep with a flag.', N),
    ('sample_id = <species>_<run code>_<original barcode>; run codes: ASPEC = Anoph_Speciation, SLCA003, SLCA011, SLCA013C2 = SLCA013 BatchCOL2, SLCA015. Use this ID (not the renumbered combined.barcodeNN file name) in SRA/BioSample submissions.', N),
    ('Example of a filled specimen block (format only, not real data): collection_site = Ahero; county = Kisumu; latitude = -0.1745; longitude = 34.9180; collection_date = 2025-03-14; collection_method = CDC light trap (indoor); morphological_id_by = initials; molecular_id_method = Koekemoer 2002 cocktail PCR; dna_extraction = Qiagen DNeasy, batch E03; biosample = SAMN........; sra_run = SRR........', N),
    ('', N),
    ('Key findings of the checks', B),
    ('1. Flow cells were reused: FBD61062 carried 5 runs (stephensi ST -> COL2 -> BatchINF (Plasmodium) -> BatchFUN -> TEP1 test); FBA66584 carried 7 (SLCA003 -> 3 stephensi trial runs -> SLCA011 -> SLCA012 -> BatchCOL).', N),
    ('2. Library from the previous run carries over on the SAME barcode: BatchFUN barcodes 01-48 (unused in BatchFUN) received 0.76-1.14% (median 0.99%) of BatchINF reads on the same barcode, with the same composition (Plasmodium % Spearman rho = 0.94, n = 48). ST -> COL2: ~0.6% of stephensi reads per barcode.', N),
    ('3. funestus files combined.barcode01 and 02 (BatchFUN raw barcode02/17, ~1,000-1,200 reads) are pure carry-over of BatchINF libraries (92-94% Plasmodium, ~2% funestus) - recommend excluding them from the species-ID set (funestus n 79 -> 77).', N),
    ('4. The other 25 BatchFUN funestus samples used barcodes 49-76, which BatchINF did not use: expected carry-over <1% (from COL2/ST).', N),
    ('5. The run NTCs (bc89/94/95) sit on barcodes never used in the earlier runs, so they cannot detect same-barcode carry-over; unused barcodes that were used in the previous run (e.g. BatchFUN bc01-48) are the informative carry-over controls.', N),
    ('6. Anoph_Speciation raw data exists (Sequence_Bkup/SequenceData): all 17 funestus (and 45 arabiensis, 1 coluzzii) files are identical to raw - now verified. SLCA003 raw is on Sequence_Bkup: 35/37 study files identical (gambiae raw bc36/bc56 folders missing from the backup); all 4 funestus PASS.', N),
    ('7. SLCA011: fastq_pass deleted, but raw signal (fast5, 48 GB) is on the backup - re-basecalling would allow verification. Every SLCA011 barcode is a mixture (median top call ~28%); unused barcodes are stephensi-dominated.', N),
    ('8. Run effects: funestus is correct in every run (SLCA003 3/4), but top-call strength differs by run (median 28% SLCA011 vs 92% BatchFUN/COL2). coluzzii is 0/18 in SLCA020 and 12/12 elsewhere; coustani 0/27 in SLCA021 and 3/3 in SLCA011 - species and run are confounded.', N),
    ('', N),
    ('Method notes', B),
    ('Composition: up to 5,000 reads per barcode subsampled evenly across each barcode\'s time-ordered fastq chunks; Kraken2 (species_id_infectivity_db, same DB as accuracy analysis), species = rank-S ancestor; Plasmodium = any taxon under genus Plasmodium. Scripts and tables: anopheles-panel/run_checks/.', N),
    ('Expected carry-over % = (1% x reads on the same barcode in the preceding run on the flow cell + 0.1% x older runs) / sample reads. >=50% = sample reads are explained by carry-over.', N),
    ('Calls/top %: six-species Bracken calls of the verified set (gamcol_validation/accuracy6/per_sample_calls.tsv).', N),
]
for t, f in lines:
    ws.append([t]); ws.cell(ws.max_row, 1).font = Font(name=AR, size=12 if f is B and ws.max_row == 1 else 10, bold=(f is B))
    ws.cell(ws.max_row, 1).alignment = Alignment(wrap_text=True, vertical='top')
ws.column_dimensions['A'].width = 160

# Funestus_samples ----------------------------------------------------------------------
ws = wb.create_sheet('Funestus_samples')
FILLC = ['collection_site', 'county', 'latitude', 'longitude', 'collection_date', 'collection_method', 'morphological_id_by',
         'molecular_id_method', 'dna_extraction_method_batch', 'flow_cell_washed_before_run', 'biosample_accession', 'sra_run_accession', 'notes']
hdr = ['sample_id', 'species_label', 'final_file', 'verified_file', 'run', 'run_started', 'flow_cell', 'flow_cell_use_no',
       'run_id', 'original_barcode', 'kit', 'minknow_version', 'basecall_model', 'reads', 'raw_verification',
       'run_NTC_barcode', 'run_NTC_reads', 'run_NTC_file', 'kmer_call', 'kmer_top_pct', 'kmer_second', 'kmer_second_pct',
       'subsample_funestus_pct', 'subsample_Plasmodium_pct', 'previous_runs_same_barcode(reads,dominant)',
       'expected_carryover_pct', 'carryover_flag', 'recommended_action'] + FILLC
rows = []
for s in samples:
    rows.append([s['sample_id'], s['species_label'], s['final_file'], s['verified_file'], s['run'], s['run_started'], s['flow_cell'],
                 s['flow_cell_use'], s['run_id'], s['raw_barcode'], s['kit'], s['minknow'], s['model'], s['reads'], s['raw_verification'],
                 s['ntc_barcode'], s['ntc_reads'], s['ntc_file'], s['call'], s['top_pct'], s['second'], s['second_pct'],
                 s['sub_fun'], s['sub_pl'], s['prev_same_bc'], s['exp_co'], s['co_flag'], s['action']] + [''] * len(FILLC))
sheet(ws, hdr, rows, widths={'sample_id': 22, 'final_file': 44, 'verified_file': 50, 'run': 38, 'run_id': 38, 'basecall_model': 34,
                             'raw_verification': 50, 'run_NTC_file': 44, 'previous_runs_same_barcode(reads,dominant)': 60,
                             'carryover_flag': 28, 'recommended_action': 70, 'notes': 30},
      fill_cols=FILLC + ['kit', 'minknow_version', 'run_NTC_barcode'],
      comments={'kit': 'From MinKNOW report; blank = report not found, please fill (SQK-NBD114.96 used for the panel).',
                'run_NTC_barcode': 'Blank = no NTC identified for this run; fill if one existed.',
                'expected_carryover_pct': '(1% x reads on same barcode in preceding run + 0.1% x older runs) / sample reads.',
                'subsample_Plasmodium_pct': 'From 5,000-read Kraken2 subsample of the raw barcode.',
                'flow_cell_use_no': 'n-th run on this flow cell (see Flow_cell_history).'})
for i, s in enumerate(samples, 2):
    if s['action'].startswith('EXCLUDE'):
        for c in ws[i][:len(hdr) - len(FILLC)]: c.fill = RED
    elif s['action'].startswith('KEEP - flag'):
        for c in ws[i][:len(hdr) - len(FILLC)]: c.fill = AMB
    elif s['action'].startswith('NTC'):
        for c in ws[i][:len(hdr) - len(FILLC)]: c.fill = GREY
    else:  # plain KEEP row: pre-filled kit/version/NTC cells are not "to fill"
        for col in ('kit', 'minknow_version', 'run_NTC_barcode'):
            c = ws.cell(i, hdr.index(col) + 1)
            if c.value not in ('', None): c.fill = PatternFill(fill_type=None)
for col in ('kmer_top_pct', 'kmer_second_pct', 'subsample_funestus_pct', 'subsample_Plasmodium_pct', 'expected_carryover_pct'):
    j = hdr.index(col) + 1
    for i in range(2, ws.max_row + 1): ws.cell(i, j).number_format = '0.0'
for i in range(2, ws.max_row + 1): ws.cell(i, hdr.index('reads') + 1).number_format = '#,##0'
nS = ws.max_row

# Runs ---------------------------------------------------------------------------------
ws = wb.create_sheet('Runs')
hdr = ['run', 'run_code', 'run_started', 'flow_cell', 'flow_cell_use_no', 'run_id', 'minknow_version', 'kit', 'basecall_model',
       'raw_data_location', 'funestus_samples', 'funestus_reads', 'recommended_exclude', 'funestus_original_barcodes',
       'NTC_barcode', 'NTC_reads', 'flow_cell_washed_before_run', 'library_prep_date', 'operator', 'notes']
rows = []
for i, r in enumerate(RUNS, 2):
    bcs = sorted(int(s['raw_barcode'][7:]) for s in samples if s['run'] == r[1] and not s['sample_id'].startswith('NTC'))
    rng = [];
    for b in bcs:
        if rng and b == rng[-1][1] + 1: rng[-1][1] = b
        else: rng.append([b, b])
    txt = ', '.join(f'{a:02d}' if a == b else f'{a:02d}-{b:02d}' for a, b in rng)
    rows.append([r[1], CODE[r[0]], r[4], r[2], r[3], r[5], r[6], r[7], r[8], r[9],
                 f'=COUNTIFS(Funestus_samples!$E$2:$E${nS},A{i},Funestus_samples!$B$2:$B${nS},"An. funestus")',
                 f'=SUMIFS(Funestus_samples!$N$2:$N${nS},Funestus_samples!$E$2:$E${nS},A{i},Funestus_samples!$B$2:$B${nS},"An. funestus")',
                 f'=COUNTIFS(Funestus_samples!$E$2:$E${nS},A{i},Funestus_samples!$AB$2:$AB${nS},"EXCLUDE*")',
                 txt, r[10], ntc_reads.get((r[0], r[10]), '') if r[10] else 'none identified', '', '', '', ''])
rows.append(['TOTAL', '', '', '', '', '', '', '', '', '', f'=SUM(K2:K{len(RUNS)+1})', f'=SUM(L2:L{len(RUNS)+1})', f'=SUM(M2:M{len(RUNS)+1})', '', '', '', '', '', '', ''])
sheet(ws, hdr, rows, widths={'run': 40, 'run_id': 38, 'basecall_model': 34, 'raw_data_location': 60, 'funestus_original_barcodes': 50, 'notes': 30},
      fill_cols=['flow_cell_washed_before_run', 'library_prep_date', 'operator', 'notes', 'kit', 'minknow_version'])
for i in range(2, ws.max_row + 1): ws.cell(i, 12).number_format = '#,##0'
for c in ws[ws.max_row]: c.font = B

# Flow_cell_history -------------------------------------------------------------------
ws = wb.create_sheet('Flow_cell_history')
sheet(ws, ['flow_cell', 'use_no', 'run', 'started', 'notes', 'in_study', 'washed_before_this_run'],
      [list(r) + [''] for r in FLOWCELL_HISTORY], widths={'run': 42, 'notes': 60, 'started': 18},
      fill_cols=['washed_before_this_run'])

# Carryover_rates ---------------------------------------------------------------------
ws = wb.create_sheet('Carryover_rates')
rate_rows = [
    ['FBD61062', 'SLCA014_BatchINF -> SLCA015_BatchFUN', 'BatchFUN bc01-48 (unused in FUN, used in INF), INF >5,000 reads', 30, 'total reads FUN / INF same barcode', 0.99, '0.76-1.14', 'Plasmodium % INF vs FUN same barcode: Spearman rho 0.94 (n=48)'],
    ['FBD61062', 'SLCA010_BatchST -> SLCA013_BatchCOL2', 'all COL2 barcodes not stephensi-dominated', 73, 'stephensi reads in COL2 / ST reads same barcode', 0.60, '90th pct 1.10', 'barcode-specific test rho 0.49, p=0.0005'],
    ['FBD61062', 'SLCA013_BatchCOL2 -> SLCA014_BatchINF', 'INF bc49-80 (unused in INF), COL2 >5,000 reads', 31, 'total reads INF / COL2 same barcode', 0.43, '', ''],
    ['FBA66584', 'stephensi trial runs -> SLCA011', 'not estimable (trial runs too small); unused SLCA011 barcodes have 1,400-8,000 reads, 25-60% stephensi', '', '', '', '', 'no barcode-specific signal; high background in every barcode'],
]
sheet(ws, ['flow_cell', 'transition', 'barcodes used', 'n', 'measure', 'median_rate_pct', 'range', 'note'], rate_rows,
      widths={'transition': 40, 'barcodes used': 60, 'measure': 44, 'note': 60})

# Carryover_tests ----------------------------------------------------------------------
ws = wb.create_sheet('Carryover_tests')
t = tsv(R + 'carryover_tests.tsv')
sheet(ws, list(t[0]), [[float(v) if k not in ('flow_cell', 'run', 'species') and v not in ('nan', '') else v for k, v in r.items()] for r in t],
      widths={'run': 26})

# Per_run_accuracy ---------------------------------------------------------------------
ws = wb.create_sheet('Per_run_accuracy')
a = tsv(R + 'per_run_accuracy.tsv')
rows = []
for i, r in enumerate(a, 2):
    rows.append([r['species'], r['run'], int(r['n']), int(r['correct']), f'=D{i}/C{i}', float(r['median_top_pct']),
                 int(r['n_excl_carryover']), int(r['correct_excl_carryover'])])
sheet(ws, ['species', 'run', 'n', 'correct', 'accuracy', 'median_top_pct', 'n_excl_carryover', 'correct_excl_carryover'], rows,
      widths={'run': 24})
for i in range(2, ws.max_row + 1): ws.cell(i, 5).number_format = '0.0%'

wb.save(OUT); print(OUT)
