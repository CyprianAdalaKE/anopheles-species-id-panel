"""Per sample: % of aligned reads that match (>=95% identity) a reference outside the labelled species' group
(group: gambiae complex = gambiae/coluzzii/arabiensis; funestus; stephensi; coustani), and % with best identity <90% to all references."""
import gzip, glob, collections, statistics as st, csv
GRP={'gambiae':'gamcx','coluzzii':'gamcx','arabiensis':'gamcx','funestus':'funestus','stephensi':'stephensi','coustani':'coustani'}
rows=[]
for pf in sorted(glob.glob('per_read/*.paf.gz')):
    true=pf.split('/')[-1].split('__')[0]
    hits=collections.defaultdict(list)
    for x in (l.split('\t') for l in gzip.open(pf,'rt')):
        if int(x[10])>=300: loc,s=x[5].split('|'); hits[x[0]].append((int(x[9]),int(x[10]),loc,s))
    c=collections.Counter(); n=0
    for r,h in hits.items():
        loc=max(h)[2]; best={}
        for m,b,l2,s in h:
            if l2==loc: best[s]=max(best.get(s,0),m/b)
        sp,v=max(best.items(),key=lambda kv:kv[1]); n+=1
        if v>=0.95: c['own' if GRP[sp]==GRP[true] else 'hi_'+GRP[sp]]+=1
        elif v<0.90: c['lt90']+=1
        else: c['90_95']+=1
    rows.append(dict(true=true,file=pf.split('__')[1][:-7],n=n,**{k:round(100*c[k]/n,1) for k in ('own','hi_gamcx','hi_funestus','hi_stephensi','hi_coustani','90_95','lt90')}))
with open('foreign_reads_per_sample.tsv','w') as o:
    w=csv.DictWriter(o,rows[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(rows)
K=('own','hi_gamcx','hi_funestus','hi_stephensi','90_95','lt90')
print('median % of reads per sample'); print('true'.ljust(12)+''.join(k.rjust(13) for k in K)+'   samples with >5% foreign')
for t in ['arabiensis','coluzzii','gambiae','funestus','stephensi','coustani']:
    sub=[r for r in rows if r['true']==t]
    fr=[sum(r[k] for k in ('hi_gamcx','hi_funestus','hi_stephensi','hi_coustani') if not k.endswith(GRP[t])) for r in sub]
    print(t.ljust(12)+''.join(f"{st.median(r[k] for r in sub):13.1f}" for k in K)+f"   {sum(f>5 for f in fr)}/{len(sub)} (max {max(fr):.1f})")
