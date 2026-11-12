#!/usr/bin/env python3
"""NCBI submission metadata for the reads analysed in the manuscript (374 sequenced specimens; 349 in the evaluation set).
BioSample: Invertebrate.1.0 package (as used for the existing samples in PRJNA1517867). SRA: one single-end ONT amplicon run per specimen.
Yellow cells must be completed before upload (per-specimen collection data are not in the sequencing records).
Output: SRA_submission_metadata_PRJNA1517867.xlsx"""
import csv
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
M = list(csv.DictReader(open('files_manifest.tsv'), delimiter='\t'))
REGION = {'arabiensis': 'western Kenya (Siaya County)', 'gambiae': 'western Kenya (Siaya County)', 'funestus': 'western Kenya (Siaya County)',
          'coustani': 'western Kenya (Siaya County)', 'stephensi': 'northern Kenya (Turkana County)', 'coluzzii': 'northern Kenya (Turkana County)'}
BS = ['*sample_name', 'sample_title', 'bioproject_accession', '*organism', '*isolate', '*collection_date', '*geo_loc_name', '*tissue',
      'dev_stage', 'sex', 'isolation_source', 'lat_lon', 'identified_by', 'description']
SRA = ['sample_name', 'library_ID', 'title', 'library_strategy', 'library_source', 'library_selection', 'library_layout', 'platform',
       'instrument_model', 'design_description', 'filetype', 'filename']
TOFILL = {'*collection_date', '*geo_loc_name', '*tissue', 'lat_lon', 'identified_by', 'sex', 'dev_stage'}
bs, sra = [], []
for m in M:
    sp = m['species']; sid = m['sample_id']
    desc = (f"Reference ID: An. {sp}. Sequencing run {m['run']}, original barcode {m['original_barcode']}, {int(m['reads']):,} reads. "
            + ('Included in the species-identification evaluation set (n = 349).' if m['qc_status'] == 'included'
               else f"Sequenced but excluded from the evaluation set by quality control ({m['exclusion_reason']})."))
    bs.append([sid, f'An. {sp} specimen {sid}, ten-locus species-identification amplicon panel', 'PRJNA1517867', f'Anopheles {sp}', sid,
               '', '', '', 'adult', 'female', f'Field-collected mosquito, {REGION[sp]}', '', '', desc])
    sra.append([sid, sid, f'Nanopore amplicon sequencing of Anopheles {sp}: ten nuclear species-identification loci and Plasmodium 18S',
                'AMPLICON', 'GENOMIC', 'PCR', 'single', 'OXFORD_NANOPORE', 'MinION',
                'Single multiplex PCR of ten Anopheles nuclear loci (811-1,051 bp) and Plasmodium 18S rDNA (846 bp); Native Barcoding Kit 96 V14 '
                '(SQK-NBD114.96); R10.4.1 flow cell; MinION Mk1B; HAC basecalling in MinKNOW; reads demultiplexed by barcode, unfiltered.',
                'fastq', m['upload_filename']])
wb = Workbook(); Y = PatternFill('solid', fgColor='FFFF00'); HF = PatternFill('solid', fgColor='1F4E78')
def sheet(ws, hdr, rows):
    ws.append(hdr)
    for i, h in enumerate(hdr, 1):
        c = ws.cell(1, i); c.font = Font(name='Arial', bold=True, color='000000' if h in TOFILL else 'FFFFFF'); c.fill = Y if h in TOFILL else HF
        ws.column_dimensions[get_column_letter(i)].width = 22 if h not in ('description', 'design_description', 'title', 'sample_title') else 60
    for r in rows: ws.append(r)
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = Font(name='Arial', size=10)
            if hdr[c.column - 1] in TOFILL: c.fill = Y
    ws.freeze_panes = 'B2'
ws = wb.active; ws.title = 'README'
for t in ['NCBI submission metadata for the reads analysed in the manuscript (all 374 sequenced specimens; qc_status in the BioSample description).',
          'Sheet BioSample_Invertebrate: NCBI BioSample package Invertebrate.1.0 (same package as the existing samples in PRJNA1517867). Columns marked * are mandatory.',
          'Yellow columns must be completed from field records before upload: collection_date (YYYY-MM-DD, or YYYY / YYYY-MM if exact date unknown), geo_loc_name (e.g. "Kenya: Siaya, Rarieda"), '
          'tissue (e.g. "whole body" or "legs and wings", as extracted), lat_lon (e.g. "0.3524 S 34.3418 E"), identified_by. dev_stage = adult and sex = female are pre-filled: confirm or correct.',
          'NCBI does not accept blank mandatory fields; if a value is genuinely unknown use "missing", "not collected" or "not provided".',
          'Sheet SRA_metadata: one single-end ONT run per BioSample; filename = file to upload (see files_manifest.tsv for source paths and md5_source_files.txt for checksums).',
          'Save each sheet as tab-delimited text (or upload the sheet in the Submission Portal) after completing the yellow columns.']:
    ws.append([t]); ws.cell(ws.max_row, 1).alignment = Alignment(wrap_text=True); ws.cell(ws.max_row, 1).font = Font(name='Arial', size=10)
ws.column_dimensions['A'].width = 160
sheet(wb.create_sheet('BioSample_Invertebrate'), BS, bs)
sheet(wb.create_sheet('SRA_metadata'), SRA, sra)
wb.save('SRA_submission_metadata_PRJNA1517867.xlsx'); print(len(bs), 'samples')
