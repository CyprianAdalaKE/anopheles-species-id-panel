"""Accuracy of top-Bracken-species calls on the verified six-species set (manuscript method, no refinement).
coluzzii and coustani: n=30 each (species_id_verified_fastqs, 2026-11-04). True label = folder (field/morphological ID)."""
import csv, math, sys
from collections import Counter
from pathlib import Path

W = Path(sys.argv[1])
SP = ["arabiensis", "coluzzii", "coustani", "funestus", "gambiae", "stephensi"]


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)


rows = []
for sp, f in csv.reader(open(W / "files.tsv"), delimiter="\t"):
    stem = f.removesuffix(".fastq")
    b, k = W / "k2" / f"{sp}__{stem}.bracken", W / "k2" / f"{sp}__{stem}.kreport"
    recs = sorted(csv.DictReader(open(b), delimiter="\t"), key=lambda r: -float(r["fraction_total_reads"])) if b.exists() else []
    recs = [r for r in recs if r["name"].startswith("Anopheles")]
    classified = sum(int(l.split("\t")[1]) for l in open(k) if l.split("\t")[3] == "R") if k.exists() else 0
    top = recs[0]["name"].split()[-1] if recs else "none"
    sec = recs[1] if len(recs) > 1 else None
    rows.append({"true": sp, "file": f, "classified_reads": classified, "pred": top,
                 "top_pct": round(100 * float(recs[0]["fraction_total_reads"]), 1) if recs else 0,
                 "second": sec["name"].split()[-1] if sec else "-",
                 "second_pct": round(100 * float(sec["fraction_total_reads"]), 1) if sec else 0})

with open(W / "per_sample_calls.tsv", "w") as o:
    wr = csv.DictWriter(o, rows[0].keys(), delimiter="\t"); wr.writeheader(); wr.writerows(rows)

out = []
p = lambda s="": (print(s), out.append(s))


def report(rows, species, title):
    n = len(rows); ok = sum(r["true"] == r["pred"] for r in rows); lo, hi = wilson(ok, n)
    p(f"== {title}")
    p(f"Samples: {n}  ({', '.join(f'{s} {sum(r['true']==s for r in rows)}' for s in species)})")
    p(f"Overall accuracy: {100*ok/n:.1f}% ({ok}/{n}), 95% CI {lo:.1f}-{hi:.1f}")
    cols = species + (["none"] if any(r["pred"] == "none" for r in rows) else [])
    p("\nConfusion matrix (rows = true label, columns = top Bracken call):")
    p("true\\pred".ljust(16) + "".join(c[:14].rjust(15) for c in cols))
    for s in species:
        c = Counter(r["pred"] for r in rows if r["true"] == s)
        p(s.ljust(16) + "".join(str(c[x]).rjust(15) for x in cols))
    p(f"\n{'species':16s}{'n':>4s}{'TP':>5s}{'sens%':>8s}{'95% CI':>14s}{'FP':>5s}{'spec%':>8s}{'95% CI':>14s}{'PPV%':>8s}")
    for s in species:
        tp = sum(r["true"] == s and r["pred"] == s for r in rows); pos = sum(r["true"] == s for r in rows)
        fp = sum(r["true"] != s and r["pred"] == s for r in rows); neg = n - pos; tn = neg - fp
        a, b = wilson(tp, pos); c, d = wilson(tn, neg)
        ppv = f"{100*tp/(tp+fp):.1f}" if tp + fp else "n/a"
        p(f"{s:16s}{pos:4d}{tp:5d}{100*tp/pos:8.1f}{f'{a:.1f}-{b:.1f}':>14s}{fp:5d}{100*tn/neg:8.1f}{f'{c:.1f}-{d:.1f}':>14s}{ppv:>8s}")
    p()


report(rows, SP, "Six species (true label = species folder)")
merge = lambda x: "gambiae/coluzzii" if x in ("gambiae", "coluzzii") else x
mrows = [dict(r, true=merge(r["true"]), pred=merge(r["pred"])) for r in rows]
report(mrows, ["arabiensis", "coustani", "funestus", "gambiae/coluzzii", "stephensi"],
       "gambiae and coluzzii merged (An. gambiae s.s./An. coluzzii)")
four = [r for r in rows if r["true"] in ("arabiensis", "funestus", "gambiae", "stephensi")]
report(four, ["arabiensis", "funestus", "gambiae", "stephensi"], "Four species only (previous table, for reference)")
wrong = [r for r in rows if r["true"] != r["pred"]]
p(f"Misclassified samples ({len(wrong)}):")
for r in wrong:
    p(f"  {r['true']:11s} {r['file']:26s} -> {r['pred']:10s} {r['top_pct']:5.1f}%  (2nd {r['second']} {r['second_pct']}%)  classified reads {r['classified_reads']}")
(W / "accuracy_summary.txt").write_text("\n".join(out) + "\n")
