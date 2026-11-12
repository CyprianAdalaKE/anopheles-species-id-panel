#!/usr/bin/env python3
"""Map every final-folder study file to (run label, raw barcode) using read headers (run_summary_*.tsv)
or, for header-less/Anoph_Speciation files, the raw read-ID trace. Output study_samples.tsv."""
import csv, re
P = '/home/genomics/Desktop/species_id/anopheles-panel/'
RUN = {'2025-06-30_SeqLib_SLCA003_ATSB_BatchKT': 'SLCA003_BatchKT', '2025-09-10-SeqLib_SLCA011': 'SLCA011',
       '2025-10-07_SeqLib_SLCA013_NMCP_BatchCOL': 'SLCA013_BatchCOL', '2025-10-07_SeqLib_SLCA013_NMCP_BatchCOL2': 'SLCA013_BatchCOL2',
       '2025-10-06_SeqLib_SLCA010_NMCP_BatchST': 'SLCA010_BatchST', '2025-10-14_SeqLib_SLCA015_ATSB_BatchFUN': 'SLCA015_BatchFUN',
       '2026-10-19_seqLib_SLCA020': 'SLCA020', '2026-10-20_SeqLib_SLCA021': 'SLCA021', 'Anoph_Speciation': 'Anoph_Speciation'}
trace = {(r['folder'], r['file']): r for r in csv.DictReader(open(P + 'run_checks/anoph_speciation_trace.tsv'), delimiter='\t')}
out = open(P + 'run_checks/study_samples.tsv', 'w')
out.write('folder\tfile\treads\trun\tflow_cell\tbarcode\tbasecall_model\tprovenance\n')
for folder in ['arabiensis', 'coluzzii_rebuilt', 'coustani_rebuilt', 'funestus', 'gambiae', 'stephensi']:
    rows = list(csv.DictReader(open(P + f'run_summary_{folder}.tsv'), delimiter='\t'))
    by = {}
    for r in rows: by.setdefault(r['file'], []).append(r)
    for f, rs in by.items():
        r = max(rs, key=lambda x: int(x['reads'])); n = sum(int(x['reads']) for x in rs)
        t = trace.get((folder, f))
        if t and t['status'] == 'PASS':
            run, bc = re.match(r'([^/]+)/(\w+)\(', t['raw_source(top3)']).groups()
            out.write(f"{folder}\t{f}\t{n}\t{run}\tFBA66226\t{bc}\t{r['basecall_model']}\traw read-ID trace (Anoph_Speciation)\n")
        else:
            prov = 'read headers' if r['run'] != 'NA' else 'no header; source run unknown'
            if len(rs) > 1: prov += f'; {len(rs)} header groups'
            out.write(f"{folder}\t{f}\t{n}\t{RUN.get(r['run'], r['run'])}\t{r['flow_cell_id']}\t{r['barcode']}\t{r['basecall_model']}\t{prov}\n")
