"""Extract the 10 nuclear amplicons (primer to primer, inclusive) from each species' sequence in the MAFFT ortholog MSAs."""
import csv
D='/home/genomics/Desktop/manuscript_docs_CAG/target_selection'
def read(f):
    s={};k=None
    for l in open(f):
        l=l.strip()
        if l.startswith('>'): k=l[1:]; s[k]=''
        elif l: s[k]+=l.upper()
    return s
rc=lambda s: s[::-1].translate(str.maketrans('ACGTN-','TGCAN-'))
SP={'gambiae':'gambiae','coluzzii':'coluzzii','coustani':'coustani','arabiensis':'arabiensis','funestus':'funestus','stephensi':'stephensi'}
prim={}
for r in csv.DictReader(open(f'{D}/multiply_outputs/multiply_outputs/2024-03-19_ANOPH-SP-INF/select/multiplex00/table.multiplex00_overview.ag-SID.csv')):
    prim.setdefault(r['target_id'],{})[r['direction']]=r['seq']
og={r['target_id']:r['ortho_id'] for r in csv.DictReader(open(f'{D}/targets/targets/anoph-speciation-targets.csv',encoding='utf-8-sig'))}
with open('amplicons_6sp.fa','w') as o:
    for tid,p in prim.items():
        if tid not in og: continue
        msa=read(f'{D}/targets/targets/msas/{og[tid]}.mafft.fasta')
        seqs={SP[[k for k in SP if k in h][0]]:s for h,s in msa.items()}
        g=seqs['gambiae']; ung=g.replace('-',''); cols=[i for i,c in enumerate(g) if c!='-']
        for F,R in ((p['F'],rc(p['R'])),(p['R'],rc(p['F']))):
            a,b=ung.find(F[:15]),ung.find(R[:15])
            if a>=0 and b>=0: break
        c0=cols[a]; c1=cols[min(b+len(R),len(ung))-1]
        for sp,s in seqs.items():
            amp=s[c0:c1+1].replace('-','')
            o.write(f'>{tid}|{sp}\n{amp}\n')
        print(tid, {sp:len(s[c0:c1+1].replace('-','')) for sp,s in seqs.items()})
