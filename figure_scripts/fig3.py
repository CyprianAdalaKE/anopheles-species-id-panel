import pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
D='/home/genomics/Desktop/manuscript_docs_CAG/target_selection'
a=pd.read_csv(f'{D}/summary_tables/summary_tables/table.all_ortholog.scores.csv')
t=pd.read_csv(f'{D}/targets/targets/anoph-speciation-targets.csv',encoding='utf-8-sig')
SEL=['AGAP005134','AGAP001804','AGAP001805','AGAP002141','AGAP002151','AGAP003240','AGAP003262','AGAP004422','AGAP004522','AGAP007827']
passed=(a.max_bp_conserved>=40)&(a.n_20bp_conserved>=10)
short=a.ortho_id.isin(t.ortho_id); sel=a.agam_id.isin(SEL)
GREY,BLUE,ORANGE,INK='#b9b8b3','#2a78d6','#eb6834','#3a3a38'
fig,ax=plt.subplots(figsize=(7.2,5.4),dpi=300)
r=a[~passed]; ax.scatter(r.max_bp_conserved,r.n_20bp_conserved,s=10,c=GREY,alpha=.5,lw=0,label=f'Not retained (n = {len(r):,})')
r=a[passed&~short]; ax.scatter(r.max_bp_conserved,r.n_20bp_conserved,s=22,c=BLUE,alpha=.75,lw=0,label=f'Passed conservation thresholds (n = {passed.sum()})')
r=a[short&~sel]; ax.scatter(r.max_bp_conserved,r.n_20bp_conserved,s=40,facecolors='white',edgecolors=BLUE,lw=1.4,label=f'Shortlisted, not taken forward (n = {len(r)})')
r=a[sel]; ax.scatter(r.max_bp_conserved,r.n_20bp_conserved,s=48,c=ORANGE,edgecolors='white',lw=1,zorder=5,label=f'Selected panel targets (n = {len(r)})')
off={'AGAP005134':(-6,4,'right'),'AGAP001804':(-6,4,'right'),'AGAP004422':(6,-1,'left'),'AGAP007827':(6,2,'left'),'AGAP002151':(6,0,'left'),
     'AGAP002141':(-6,4,'right'),'AGAP001805':(6,-2,'left'),'AGAP004522':(6,-8,'left'),'AGAP003240':(-6,-10,'right'),'AGAP003262':(-6,-8,'right')}
for _,p in r.iterrows():
    if p.agam_id=='AGAP001805': continue
    dx,dy,ha=off[p.agam_id]
    lab='AGAP002141, AGAP001805' if p.agam_id=='AGAP002141' else p.agam_id
    LEAD={'AGAP004522':(80,23),'AGAP002141':(92,6.5),'AGAP003240':(78,3.5),'AGAP003262':(22,25)}
    if p.agam_id in LEAD:
        ax.annotate(lab,(p.max_bp_conserved,p.n_20bp_conserved),xytext=LEAD[p.agam_id],textcoords='data',ha='left',va='center',fontsize=6.5,color=INK,
                    arrowprops=dict(arrowstyle='-',color=INK,lw=.5,shrinkA=1,shrinkB=3))
    else:
        ax.annotate(lab,(p.max_bp_conserved,p.n_20bp_conserved),xytext=(dx,dy),textcoords='offset points',ha=ha,fontsize=6.5,color=INK)
ax.axvline(40,color=INK,lw=.8,ls='--'); ax.axhline(10,color=INK,lw=.8,ls='--')
ax.set_xlabel('Longest conserved region (bp)',fontsize=9,color=INK); ax.set_ylabel('Conserved segments longer than 20 bp',fontsize=9,color=INK)
ax.tick_params(labelsize=8,colors=INK); [s.set_visible(False) for k,s in ax.spines.items() if k in('top','right')]
ax.grid(color='#e6e5e0',lw=.5); ax.set_axisbelow(True); ax.set_xlim(0,128); ax.set_ylim(-1.5,44)
ax.legend(fontsize=7,frameon=False,loc='upper left')
fig.tight_layout(); fig.savefig('figure3_ortholog_screen.png'); fig.savefig('figure3_ortholog_screen.pdf')
print(passed.sum(), short.sum(), sel.sum())
print(a[sel][['agam_id','max_bp_conserved','n_20bp_conserved']].to_string())
