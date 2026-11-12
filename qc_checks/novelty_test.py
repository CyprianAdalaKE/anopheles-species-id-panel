#!/usr/bin/env python3
"""Novelty ("not in panel") rule: per sample, median best read identity to the 60 species-specific amplicon
references (5,000-read PAFs, alignment blocks >= 300 bp). Samples whose median best identity is below a cut-off
are reported as 'unresolved' instead of being forced into one of the six species.
Writes novelty_per_sample.tsv and prints the distribution by species/run and the effect of cut-offs on accuracy."""
import gzip, glob, os, csv, collections, statistics as st, re
RA = '/home/genomics/Desktop/species_id/anopheles-panel/analysis_verified/reviewer_analyses/per_read'
V = '/media/genomics/Sequence_Bkup/species_id_verified_fastqs'
calls = {(r['true'], re.sub(r'\.fastq$', '', r['file'])): r for r in
         csv.DictReader(open('/home/genomics/Desktop/species_id/gamcol_validation/accuracy6/per_sample_calls.tsv'), delimiter='\t')}
runs = {}
for m in csv.DictReader(open(f'{V}/MANIFEST.tsv'), delimiter='\t'):
    runs[(m['species'], m['file'].replace('.fastq.gz', ''))] = m['raw_source_run(reads)'].split('(')[0]
rows = []
for sp_f in sorted(runs):
    sp, f = sp_f; pf = f'{RA}/{sp}__{f}.paf.gz'
    best = collections.defaultdict(float); bestm = collections.defaultdict(int)
    for l in gzip.open(pf, 'rt'):
        x = l.split('\t'); m, blk = int(x[9]), int(x[10])
        if blk < 300: continue
        if m > bestm[x[0]]: bestm[x[0]] = m; best[x[0]] = m / blk   # best locus = most matching bases
        elif m == bestm[x[0]]: best[x[0]] = max(best[x[0]], m / blk)
    ids = sorted(best.values()); c = calls[sp_f]
    rows.append(dict(species=sp, file=f, run=runs[sp_f], reads_aligned=len(ids),
                     median_best_identity=round(100 * st.median(ids), 2) if ids else '',
                     pct_ge95=round(100 * sum(i >= 0.95 for i in ids) / len(ids), 1) if ids else '',
                     pred=c['pred'], top_pct=c['top_pct'], correct=int(c['pred'] == sp)))
with open('novelty_per_sample.tsv', 'w') as o:
    w = csv.DictWriter(o, list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
g = collections.defaultdict(list)
for r in rows: g[(r['species'], r['run'])].append(r['median_best_identity'])
print('median best identity by species / run: n, min, median, max')
for k, v in sorted(g.items()): print(f'  {k[0]:11s} {k[1]:42s} {len(v):3d} {min(v):6.2f} {st.median(v):6.2f} {max(v):6.2f}')
print('\ncut-off  unresolved(by species)  accuracy among resolved  correct calls lost')
for cut in (93, 94, 94.5, 95, 95.5, 96):
    un = [r for r in rows if r['median_best_identity'] < cut]; res = [r for r in rows if r['median_best_identity'] >= cut]
    print(f'  {cut:5}  {len(un):3d} {dict(collections.Counter(r["species"] for r in un))}  '
          f'{100*sum(r["correct"] for r in res)/len(res):.1f}% ({sum(r["correct"] for r in res)}/{len(res)})  {sum(r["correct"] for r in un)}')
