import csv,glob,re,itertools,json
D='/home/genomics/Desktop/manuscript_docs_CAG/target_selection'
def read(f):
    s={};k=None
    for l in open(f):
        l=l.strip()
        if l.startswith('>'): k=l[1:]; s[k]=''
        elif l: s[k]+=l.upper()
    return s
rc=lambda s: s[::-1].translate(str.maketrans('ACGTN-','TGCAN-'))
SP={'gambiae':'gam','coluzzii':'col','coustani':'cou','arabiensis':'ara','funestus':'fun','stephensi':'ste'}
prim={}
for r in csv.DictReader(open(f'{D}/multiply_outputs/multiply_outputs/2024-03-19_ANOPH-SP-INF/select/multiplex00/table.multiplex00_overview.ag-SID.csv')):
    prim.setdefault(r['target_id'],{})[r['direction']]=r['seq']
og={r['target_id']:r['ortho_id'] for r in csv.DictReader(open(f'{D}/targets/targets/anoph-speciation-targets.csv',encoding='utf-8-sig'))}
out={}
for tid,p in prim.items():
    if tid not in og: continue
    msa=read(f'{D}/targets/targets/msas/{og[tid]}.mafft.fasta')
    seqs={SP[[k for k in SP if k in h][0]]:s for h,s in msa.items()}
    g=seqs['gam']; ung=g.replace('-',''); cols=[i for i,c in enumerate(g) if c!='-']
    orient=None
    for o,(F,R) in (('+',(p['F'],rc(p['R']))),('-',(p['R'],rc(p['F'])))):
        a,b=ung.find(F[:15]),ung.find(R[:15])
        if a>=0 and b>=0: orient=o; f0,r0=a,b; fl=min(len(F),len(ung)-a); rl=min(len(R),len(ung)-b); break
    if orient is None: print(tid,'primer not found'); continue
    fcols=cols[f0:f0+fl]; rcols=cols[r0:r0+rl]; amp=range(fcols[0],rcols[-1]+1)
    res={'orient':orient,'amp_cols':len(amp)}
    # primer mismatches per species; 3' end position
    for name,pc,three in (('F' if orient=='+' else 'R',fcols,'end'),('R' if orient=='+' else 'F',rcols,'start')):
        mm={}
        for sp,s in seqs.items():
            pos=[i for i,c in enumerate(pc) if s[c]!=g[c]]
            # distance from 3' end: for left primer 3' is last col; right primer (rc) 3' is first col
            d3=[(len(pc)-1-i) if three=='end' else i for i in pos]
            mm[sp]={'n':len(pos),'gap':sum(s[c]=='-' for c in pc),'min_from3':min(d3) if d3 else None}
        res[name]=mm
    # interspecies differences inside amplicon (excluding primers), ignoring gap-vs-base? count both
    inner=[c for c in amp if c not in set(fcols)|set(rcols)]
    pw={}
    for a,b in itertools.combinations(['gam','col','ara','fun','ste','cou'],2):
        sub=sum(1 for c in inner if seqs[a][c] in 'ACGT' and seqs[b][c] in 'ACGT' and seqs[a][c]!=seqs[b][c])
        gap=sum(1 for c in inner if (seqs[a][c]=='-')!=(seqs[b][c]=='-'))
        pw[f'{a}-{b}']=(sub,gap)
    res['pw']=pw
    out[tid]=res
json.dump(out,open('primer_msa.json','w'),indent=1)
for t,r in out.items():
    fm={k:(v['n'],v['min_from3']) for k,v in r['F'].items() if v['n']}; rm={k:(v['n'],v['min_from3']) for k,v in r['R'].items() if v['n']}
    print(t,r['orient'],'F mm',fm,'R mm',rm,'| gam-col',r['pw']['gam-col'],'gam-ara',r['pw']['gam-ara'],'gam-fun',r['pw']['gam-fun'],'gam-ste',r['pw']['gam-ste'],'gam-cou',r['pw']['gam-cou'])
