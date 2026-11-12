"""In-silico PCR: for each species/locus, find both panel primers (best Hamming match, <=6 mm) within +/-700 bp of the
MSA-derived amplicon's own-genome hit, and extract the primer-to-primer genomic sequence."""
import csv, subprocess, re
G='/home/genomics/species_id/species_id_infectivity_db/genomes'
GF=dict(gambiae='AgambiaePEST',coluzzii='AcoluzziiAcolN3',coustani='AcoustaniAcouGA1',arabiensis='AarabiensisDONGOLA2021',funestus='AfunestusAfunGA1',stephensi='AstephensiUCISS2018')
D='/home/genomics/Desktop/manuscript_docs_CAG/target_selection'
rc=lambda s: s[::-1].translate(str.maketrans('ACGTN','TGCAN'))
prim={}
for r in csv.DictReader(open(f'{D}/multiply_outputs/multiply_outputs/2024-03-19_ANOPH-SP-INF/select/multiplex00/table.multiplex00_overview.ag-SID.csv')):
    prim.setdefault(r['target_id'],{})[r['direction']]=r['seq']
def best(seq,p):
    b=(99,None)
    for i in range(len(seq)-len(p)+1):
        d=sum(1 for x,y in zip(seq[i:i+len(p)],p) if x!=y)
        if d<b[0]: b=(d,i)
    return b
out=open('amplicons_6sp_insilico.fa','w'); log=open('insilico_pcr.tsv','w'); log.write('locus\tspecies\tlen\tF_mm\tR_mm\n')
for sp,gf in GF.items():
    fa=f'{G}/VectorBase-67_{gf}_Genome.k2.fasta'
    for l in open(f'hits_{sp}.paf'):
        q,ql,qs,qe,ctg,ts,te,*_=l.split('\t'); loc=q.split('|')[0]
        a,b=max(0,int(ts)-700),int(te)+700
        reg=subprocess.run(['samtools','faidx',fa,f'{ctg}:{a+1}-{b}'],capture_output=True,text=True).stdout
        s=''.join(reg.split('\n')[1:]).upper()
        F,R=prim[loc]['F'],prim[loc]['R']
        res=[]
        for strand,(P1,P2) in (('+',(F,rc(R))),('-',(R,rc(F)))):
            d1,i1=best(s,P1); d2,i2=best(s,P2)
            if i1 is not None and i2 is not None and i2>i1: res.append((d1+d2,d1,d2,s[i1:i2+len(P2)],strand))
        d,d1,d2,amp,strand=min(res)
        if strand=='-': amp=rc(amp); d1,d2=d2,d1
        out.write(f'>{loc}|{sp}\n{amp}\n'); log.write(f'{loc}\t{sp}\t{len(amp)}\t{d1}\t{d2}\n')
out.close(); log.close()
