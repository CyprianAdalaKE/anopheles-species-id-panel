#!/usr/bin/env python3
"""Rows for Additional file 1 Tables S4 (read identity) and S5 (per-locus accuracy) from foreign_reads_per_sample.tsv,
identity_profile.txt and per_locus_calls.tsv / per_sample_alignment_kraken_calls.tsv. Usage: make_tables_S4_S5.py <dir> <out.json>"""
import csv, json, re, sys, statistics as st, collections
D = sys.argv[1]; SP = ['arabiensis', 'coluzzii', 'coustani', 'funestus', 'gambiae', 'stephensi']
fr = list(csv.DictReader(open(f'{D}/foreign_reads_per_sample.tsv'), delimiter='\t'))
idp = {m.group(1): (m.group(2), m.group(3)) for m in re.finditer(r'^(\w+)\s+reads\s+(\d+)\s+best-identity median ([\d.]+)%', open(f'{D}/identity_profile.txt').read(), re.M)}
S4 = [["Labelled species", "Samples", "Reads analysed", "Median best identity (%)", "Reads ≥95% to own species (%)",
       "Reads ≥95% to another panel species (%)", "Range across samples (%)", "Reads <95% to all references (%)"]]
for s in SP:
    r = [x for x in fr if x['true'] == s]
    oth = [sum(float(x[k]) for k in x if k.startswith('hi_')) for x in r]
    S4.append([f'*An. {s}*', str(len(r)), f'{int(idp[s][0]):,}', idp[s][1], f"{st.median(float(x['own']) for x in r):.1f}",
               f'{st.median(oth):.1f}', f'{min(oth):.1f}–{max(oth):.1f}', f"{st.median(float(x['90_95']) + float(x['lt90']) for x in r):.1f}"])
L = list(csv.DictReader(open(f'{D}/per_locus_calls.tsv'), delimiter='\t'))
R = list(csv.DictReader(open(f'{D}/per_sample_alignment_kraken_calls.tsv'), delimiter='\t'))
base = sum(r['true'] == r['kraken_vote_call'] for r in R) / len(R)
S5 = [["Locus", "*An. arab.*", "*An. col.*", "*An. cou.*", "*An. fun.*", "*An. gam.*", "*An. ste.*", "All samples", "Accuracy without locus (change)"]]
for loc in sorted({r['locus'] for r in L}):
    sub = [r for r in L if r['locus'] == loc and int(r.get('n_reads', r.get('n', 0)) or 0) >= 10] if 'n_reads' in L[0] or 'n' in L[0] else [r for r in L if r['locus'] == loc]
    allL = [r for r in L if r['locus'] == loc]
    by = [f"{100*sum(r['true']==r['kraken_locus_call'] for r in sub if r['true']==s)/max(1,sum(r['true']==s for r in sub)):.0f}" for s in SP]
    a = sum(r['true'] == r['kraken_locus_call'] for r in sub) / len(sub); lo = sum(r['true'] == r['kraken_lolo_call'] for r in allL) / len(allL)
    d = 100 * (lo - base); S5.append([loc] + by + [f'{100*a:.1f}', f"{100*lo:.1f} ({'−' if d < -0.05 else ''}{abs(d):.1f})"])
S5.append(["All loci", "", "", "", "", "", "", f'{100*base:.1f}', ""])
json.dump({'S4': S4, 'S5': S5}, open(sys.argv[2], 'w'), indent=1, ensure_ascii=False)
for t in (S4, S5):
    for r in t: print(' | '.join(r))
