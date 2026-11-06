"""Per-locus Kraken2 vote accuracy by species; per-locus read identity of An. coustani-labelled samples to each reference;
per-sample spread of identity for coustani samples."""
import csv, gzip, glob, os, collections, statistics as st
SP=['arabiensis','coluzzii','coustani','funestus','gambiae','stephensi']
L=list(csv.DictReader(open('per_locus_calls.tsv'),delimiter='\t'))
loci=sorted({r['locus'] for r in L})
print('Kraken per-locus vote accuracy by species (%), samples with >=10 reads at locus')
print('locus'.ljust(12)+''.join(s[:5].rjust(7) for s in SP)+'  all')
for l in loci:
    sub=[r for r in L if r['locus']==l and r['kraken_locus_call']!='none']
    print(l.ljust(12)+''.join(f"{100*sum(r['true']==r['kraken_locus_call'] for r in sub if r['true']==s)/max(1,sum(r['true']==s for r in sub)):7.0f}" for s in SP)+f"{100*sum(r['true']==r['kraken_locus_call'] for r in sub)/len(sub):6.1f}")
# coustani & funestus samples: per-locus median read identity to each reference
for true in ('coustani','funestus'):
    ids=collections.defaultdict(list); persample=[]
    for pf in glob.glob(f'per_read/{true}__*.paf.gz'):
        hits=collections.defaultdict(list)
        for x in (l.split('\t') for l in gzip.open(pf,'rt')):
            if int(x[10])>=300: loc,s=x[5].split('|'); hits[x[0]].append((int(x[9]),int(x[10]),loc,s))
        smp=collections.defaultdict(list)
        for r,h in hits.items():
            loc=max(h)[2]; best={}
            for m,b,l2,s in h:
                if l2==loc: best[s]=max(best.get(s,0),m/b)
            for s,v in best.items(): ids[(loc,s)].append(v); smp[s].append(v)
        persample.append({s:st.mean(v) for s,v in smp.items()})
    print(f'\n{true}-labelled samples: median read identity (%) to each reference, by locus')
    print('locus'.ljust(12)+''.join(s[:5].rjust(7) for s in SP))
    for l in loci: print(l.ljust(12)+''.join(f"{100*st.median(ids[(l,s)]):7.1f}" if ids[(l,s)] else '      -' for s in SP))
    for s in ('coustani','funestus','gambiae'):
        v=sorted(100*d.get(s,0) for d in persample); print(f'  per-sample mean identity to {s}: min {v[0]:.1f} median {st.median(v):.1f} max {v[-1]:.1f}')
    best_own=sum(1 for d in persample if max(d,key=d.get)==true); print(f'  samples whose highest mean identity is to the {true} reference: {best_own}/{len(persample)}; to funestus: {sum(1 for d in persample if max(d,key=d.get)=="funestus")}')
