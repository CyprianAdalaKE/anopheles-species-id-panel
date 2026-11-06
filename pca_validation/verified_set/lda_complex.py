"""Within-complex (gambiae/coluzzii/arabiensis) discrimination from panel SNP dosages:
per-locus Hudson Fst (from locus_fst_table.tsv) and 5-fold stratified CV LDA accuracy, per locus and all 10 loci
(PCA inside each fold, 10 PCs, then LDA). Majority-class baseline reported."""
import csv, numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score
man = {r["sample_id"]: r["species"] for r in csv.DictReader(open("manifest.csv"))}
order = [l.strip() for l in open("results/sample_order.txt") if l.strip()]
ids = [f"{p.split('/')[-2]}_{p.split('/')[-1].split('.')[0]}" for p in order]
sp = np.array([man[i] for i in ids])
def dos(g):
    if g.startswith("."): return np.nan
    a, b = g.replace("|", "/").split("/"); return np.nan if "." in (a, b) else int(a) + int(b)
pos, rows = [], []
for l in open("results/genotype_matrix.tsv"):
    c = l.rstrip("\n").split("\t"); pos.append((c[0], int(c[1]))); rows.append([dos(g) for g in c[4:]])
G = np.array(rows, float).T
bed = [l.split() for l in open("/home/genomics/Desktop/species_id/anopheles-panel/nomadsSID.amplicons.bed")]
keep = np.isin(sp, ["gambiae", "coluzzii", "arabiensis"]); G, y = G[keep], sp[keep]
G = np.where(np.isnan(G), np.nanmean(G, 0), G)
cv = StratifiedKFold(5, shuffle=True, random_state=1)
base = max(np.mean(y == s) for s in set(y))
print(f"n={len(y)} ({', '.join(f'{s} {sum(y==s)}' for s in sorted(set(y)))}); majority baseline {100*base:.1f}%")
out = []
for ch, s, e, loc in bed:
    m = np.array([(c == ch and int(s) <= p <= int(e)) for c, p in pos])
    X = G[:, m]; X = X[:, X.std(0) > 0]
    acc = cross_val_score(make_pipeline(StandardScaler(), LDA(solver="lsqr", shrinkage="auto")), X, y, cv=cv)
    out.append((loc, X.shape[1], acc.mean(), acc.std()))
X = G[:, G.std(0) > 0]
acc = cross_val_score(make_pipeline(StandardScaler(), PCA(10, random_state=1), LDA()), X, y, cv=cv)
with open("results/locus_lda_5fold.tsv", "w") as o:
    o.write("locus\tn_snps\tLDA_5fold_accuracy\tsd\n")
    for r in sorted(out, key=lambda r: -r[2]): o.write(f"{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]:.4f}\n"); print(f"{r[0]}  {r[1]:4d} SNPs  {100*r[2]:.1f}% ±{100*r[3]:.1f}")
    o.write(f"ALL_10_LOCI_PCA10_LDA\t{X.shape[1]}\t{acc.mean():.4f}\t{acc.std():.4f}\nmajority_baseline\t\t{base:.4f}\t\n")
print(f"All loci (PCA10 + LDA): {100*acc.mean():.1f}% ±{100*acc.std():.1f}")

# breakdown: CV-predicted confusion (all loci) and two-way gambiae vs coluzzii
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix, roc_auc_score
pipe = make_pipeline(StandardScaler(), PCA(10, random_state=1), LDA())
pred = cross_val_predict(pipe, X, y, cv=cv); labs = ["arabiensis", "coluzzii", "gambiae"]
print("confusion (rows true):", labs); print(confusion_matrix(y, pred, labels=labs))
m2 = np.isin(y, ["gambiae", "coluzzii"]); X2, y2 = X[m2], y[m2]
cv2 = StratifiedKFold(5, shuffle=True, random_state=1)
a2 = cross_val_score(pipe, X2, y2, cv=cv2)
p2 = cross_val_predict(pipe, X2, y2, cv=cv2, method="predict_proba")[:, 0]
print(f"gambiae vs coluzzii: 5-fold accuracy {100*a2.mean():.1f}% (baseline {100*max(np.mean(y2=='gambiae'),np.mean(y2=='coluzzii')):.1f}%), AUC {roc_auc_score(y2=='coluzzii', p2):.2f}")
print("confusion:", confusion_matrix(y2, cross_val_predict(pipe, X2, y2, cv=cv2), labels=["coluzzii","gambiae"]).tolist())
