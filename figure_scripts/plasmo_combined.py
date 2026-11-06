"""Combined Plasmodium-specificity figure (replaces former Figs 12 and 13): A every sample by group,
B share of samples called positive across all thresholds, C classified reads vs % P. falciparum.
Data: plasmo_rerun_taxfix/per_sample.tsv (corrected database)."""
import csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
rows=[r for r in csv.DictReader(open('/home/genomics/Desktop/species_id/plasmo_rerun_taxfix/per_sample.tsv'),delimiter='\t') if r['corrected_pf_pct']!='']
for r in rows: r['pf']=float(r['corrected_pf_pct']); r['n']=int(r['corrected_classified'])
G=[('2 negative ctrl','Negative controls\n(mosquito DNA only)','#2a78d6'),
   ('1 dilution rep1','3D7 dilution series,\nreplicate 1','#eb6834'),
   ('3 dilution rep2','3D7 dilution series,\nreplicate 2','#1baf7a'),
   ('4 wild-infected','Wild-caught\nmosquitoes','#eda100'),
   ('5 SID primers+3D7','3D7 + species-ID\nprimers only','#e87ba4')]
by={g:[r for r in rows if r['group']==g] for g,_,_ in G}
INK,INK2,GRID='#1f1f1d','#5c5b56','#e6e5e0'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.edgecolor':INK2,'xtick.color':INK2,'ytick.color':INK2,'axes.labelcolor':INK})
fig=plt.figure(figsize=(10.5,7.2),dpi=300); gs=fig.add_gridspec(2,2,height_ratios=[1,1],hspace=.42,wspace=.25)
axA=fig.add_subplot(gs[0,:]); axB=fig.add_subplot(gs[1,0]); axC=fig.add_subplot(gs[1,1])
rng=np.random.default_rng(3)
for i,(g,lab,c) in enumerate(G):
    v=np.array([r['pf'] for r in by[g]])
    axA.scatter(i+rng.uniform(-.17,.17,len(v)),v,s=34,color=c,edgecolor='white',lw=1.2,zorder=3)
    axA.hlines(np.median(v),i-.28,i+.28,color=INK,lw=1.6,zorder=4)
axA.set_xticks(range(len(G)),[f"{lab}\n(n = {len(by[g])})" for g,lab,_ in G])
axA.set_ylabel(r'Reads classified as $\it{P.\ falciparum}$ (%)'); axA.set_ylim(-4,102); axA.set_xlim(-.6,len(G)-.4)
# B threshold sweep
t=np.linspace(0,100,1001)
for g,lab,c in G:
    v=np.array([r['pf'] for r in by[g]])
    axB.plot(t,[100*np.mean(v>x) for x in t],color=c,lw=2,label=lab.replace('\n',' '))
axB.set_xlabel(r'Positivity threshold (% reads classified as $\it{P.\ falciparum}$)'); axB.set_ylabel('Samples called positive (%)')
axB.set_xlim(0,100); axB.set_ylim(-3,103)
# C depth vs pct
for g,lab,c in G:
    d=by[g]; axC.scatter([r['n'] for r in d],[r['pf'] for r in d],s=34,color=c,edgecolor='white',lw=1.2,zorder=3,label=lab.replace('\n',' '))
axC.set_xscale('log'); axC.set_xlabel('Reads with a taxonomic classification (log scale)'); axC.set_ylabel(r'Reads classified as $\it{P.\ falciparum}$ (%)'); axC.set_ylim(-4,102)
for ax,L in ((axA,'A'),(axB,'B'),(axC,'C')):
    ax.grid(color=GRID,lw=.6); ax.set_axisbelow(True)
    for s in ('top','right'): ax.spines[s].set_visible(False)
    ax.text(-0.06 if ax is not axA else -0.035,1.04,L,transform=ax.transAxes,fontsize=12,fontweight='bold',color=INK)
h,l=axC.get_legend_handles_labels()
fig.legend(h,l,loc='lower center',ncol=5,frameon=False,fontsize=7.5,bbox_to_anchor=(.5,-.005),markerscale=1.2)
fig.subplots_adjust(bottom=.12,top=.95,left=.07,right=.98)
fig.savefig('figure_plasmo_combined.png'); fig.savefig('figure_plasmo_combined.pdf')
# numbers for text
neg=np.array([r['pf'] for r in by['2 negative ctrl']]); wild=np.array([r['pf'] for r in by['4 wild-infected']])
dil=np.array([r['pf'] for g in ('1 dilution rep1','3 dilution rep2') for r in by[g]])
for x in (1,5,10,20,30,40,50):
    print(f"thr {x}%: neg {np.sum(neg>x)}/10 wild {np.sum(wild>x)}/10 dil {np.sum(dil>x)}/20")
print('max wild',wild.max(),'neg above max wild',np.sum(neg>wild.max()))
