"""For each true species: share of reads per locus, and per-read best identity (max over the six references at the read's best locus),
with the reference species giving that best identity."""
import gzip, glob, collections, statistics as st
SP=['arabiensis','coluzzii','coustani','funestus','gambiae','stephensi']
for true in SP:
    loc_n=collections.Counter(); bestid=[]; hi_sp=collections.Counter(); n=0
    for pf in glob.glob(f'per_read/{true}__*.paf.gz'):
        hits=collections.defaultdict(list)
        for x in (l.split('\t') for l in gzip.open(pf,'rt')):
            if int(x[10])>=300: loc,s=x[5].split('|'); hits[x[0]].append((int(x[9]),int(x[10]),loc,s))
        for r,h in hits.items():
            loc=max(h)[2]; best={}
            for m,b,l2,s in h:
                if l2==loc: best[s]=max(best.get(s,0),m/b)
            sp,idv=max(best.items(),key=lambda kv:kv[1]); n+=1; loc_n[loc]+=1; bestid.append(idv)
            if idv>=0.95: hi_sp[sp]+=1
    q=sorted(bestid)
    print(f"{true:11s} reads {n:7d}  best-identity median {100*st.median(q):.1f}%  >=95%: {100*sum(v>=0.95 for v in q)/n:5.1f}%  (those best match: {dict(hi_sp.most_common(3))})")
    print('            reads by locus %: '+' '.join(f"{l[4:]}:{100*c/n:.0f}" for l,c in sorted(loc_n.items())))
