"""Summarise per-read results (5,000 reads/sample):
 (1) alignment-based species assignment: each read's best locus = hit with most matching bases; read assigned to the species whose
     amplicon reference has the highest identity (matches / alignment block length) at that locus; ties = ambiguous.
 (2) Kraken2 per-read species assignment (same reads), majority vote over species-level reads.
 Per-sample, per-locus and leave-one-locus-out calls for both methods; identity of each read to every species' reference."""
import gzip, glob, os, csv, collections, statistics as st
SP = ['arabiensis','coluzzii','coustani','funestus','gambiae','stephensi']
TAX = {'7173':'arabiensis','1518534':'coluzzii','139045':'coustani','62324':'funestus','7165':'gambiae','30069':'stephensi'}
LOCI = sorted({l[1:].split('|')[0] for l in open('amplicons_ref.fa') if l.startswith('>')})
MINBLOCK = 300
vote = lambda c: max(c, key=c.get) if c and max(c.values()) > 0 else 'none'
rows, locrows, idrows = [], [], []
for pf in sorted(glob.glob('per_read/*.paf.gz')):
    base = os.path.basename(pf)[:-7]; true, f = base.split('__')
    k2 = dict(l.rstrip('\n').split('\t') for l in gzip.open(base.join(['per_read/', '.k2.gz']), 'rt'))
    hits = collections.defaultdict(list)
    for l in gzip.open(pf, 'rt'):
        x = l.split('\t'); m, blk = int(x[9]), int(x[10])
        if blk < MINBLOCK: continue
        loc, s = x[5].split('|'); hits[x[0]].append((m, blk, loc, s))
    reads = []  # (locus, align_species or 'ambig', kraken_species or None, ident dict)
    for r, h in hits.items():
        loc = max(h)[2]; ident = {}
        for m, blk, l2, s in h:
            if l2 == loc: ident[s] = max(ident.get(s, 0), m / blk)
        srt = sorted(ident.items(), key=lambda kv: -kv[1])
        a = srt[0][0] if len(srt) == 1 or srt[0][1] > srt[1][1] else 'ambig'
        reads.append((loc, a, TAX.get(k2.get(r)), ident))
    def calls(sub):
        ca = collections.Counter(a for _, a, _, _ in sub if a in SP)
        ck = collections.Counter(k for _, _, k, _ in sub if k)
        return vote(ca), vote(ck), ca, ck
    a_call, k_call, ca, ck = calls(reads)
    n = len(reads)
    rows.append(dict(true=true, file=f, aligned_reads=n, align_call=a_call, kraken_vote_call=k_call,
                     **{f'align_{s}': round(100 * ca[s] / n, 1) if n else 0 for s in SP},
                     ambiguous_pct=round(100 * sum(a == 'ambig' for _, a, _, _ in reads) / n, 1) if n else 0,
                     **{f'meanid_{s}': round(100 * st.mean(i[s] for *_, i in reads if s in i), 2) if any(s in i for *_, i in reads) else '' for s in SP}))
    for L in LOCI:
        sub = [r for r in reads if r[0] == L]; rest = [r for r in reads if r[0] != L]
        a1, k1, _, _ = calls(sub) if len(sub) >= 10 else ('none', 'none', 0, 0)
        a2, k2c, _, _ = calls(rest)
        locrows.append(dict(true=true, file=f, locus=L, n_reads=len(sub), align_locus_call=a1, kraken_locus_call=k1,
                            align_lolo_call=a2, kraken_lolo_call=k2c))
def w(name, rr):
    with open(name, 'w') as o:
        d = csv.DictWriter(o, rr[0].keys(), delimiter='\t'); d.writeheader(); d.writerows(rr)
w('per_sample_alignment_kraken_calls.tsv', rows); w('per_locus_calls.tsv', locrows)

out = []
p = lambda s='': (print(s), out.append(s))
br = {(r['true'], r['file'].removesuffix('.fastq')): r['pred'] for r in csv.DictReader(open('/home/genomics/Desktop/species_id/gamcol_validation/accuracy6/per_sample_calls.tsv'), delimiter='\t')}
br9 = {(r['true'], r['file'].removesuffix('.fastq')): r['pred'] for r in csv.DictReader(open('bracken900/per_sample_calls.tsv'), delimiter='\t')}
for r in rows:
    k = (r['true'], r['file'].removesuffix('.fastq.gz')); r['bracken100'] = br.get(k, '?'); r['bracken900'] = br9.get(k, '?')
merge = lambda s: 'gam/col' if s in ('gambiae', 'coluzzii') else s
for meth in ('bracken100', 'bracken900', 'kraken_vote_call', 'align_call'):
    ok = sum(r['true'] == r[meth] for r in rows); okm = sum(merge(r['true']) == merge(r[meth]) for r in rows)
    p(f'== {meth}: six-species {100*ok/len(rows):.1f}% ({ok}/{len(rows)}); gam+col merged {100*okm/len(rows):.1f}% ({okm})')
    p('true\\pred'.ljust(12) + ''.join(s[:10].rjust(11) for s in SP + ['none']))
    for s in SP:
        c = collections.Counter(r[meth] for r in rows if r['true'] == s)
        p(s.ljust(12) + ''.join(str(c[x]).rjust(11) for x in SP + ['none']))
p()
p('== Mean read identity (%) to each species reference, by true species (median of per-sample means)')
p('true'.ljust(12) + ''.join(s[:10].rjust(11) for s in SP))
for s in SP:
    p(s.ljust(12) + ''.join(f"{st.median(float(r['meanid_'+t]) for r in rows if r['true']==s and r['meanid_'+t]!=''):11.2f}" for t in SP))
p()
p('== Per-locus sample-level accuracy (vote of reads at that locus; samples with >=10 reads at the locus)')
p('locus'.ljust(12) + 'n'.rjust(5) + 'align%'.rjust(8) + 'kraken%'.rjust(9) + '  align% by species ' + ' '.join(s[:4] for s in SP) + '   LOLO align% / kraken%')
base_a = sum(r['true'] == r['align_call'] for r in rows) / len(rows); base_k = sum(r['true'] == r['kraken_vote_call'] for r in rows) / len(rows)
for L in LOCI:
    sub = [r for r in locrows if r['locus'] == L and r['align_locus_call'] != 'none']
    allL = [r for r in locrows if r['locus'] == L]
    a = sum(r['true'] == r['align_locus_call'] for r in sub) / len(sub); k = sum(r['true'] == r['kraken_locus_call'] for r in sub) / len(sub)
    bysp = ' '.join(f"{100*sum(r['true']==r['align_locus_call'] for r in sub if r['true']==s)/max(1,sum(r['true']==s for r in sub)):4.0f}" for s in SP)
    la = sum(r['true'] == r['align_lolo_call'] for r in allL) / len(allL); lk = sum(r['true'] == r['kraken_lolo_call'] for r in allL) / len(allL)
    p(f'{L:12s}{len(sub):5d}{100*a:8.1f}{100*k:9.1f}  {bysp}   {100*la:.1f} ({100*(la-base_a):+.1f}) / {100*lk:.1f} ({100*(lk-base_k):+.1f})')
p(f'all loci: align {100*base_a:.1f}%, kraken vote {100*base_k:.1f}%')
open('summary.txt', 'w').write('\n'.join(out) + '\n')
