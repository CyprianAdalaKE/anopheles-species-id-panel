import os, glob, csv

BASE = os.path.join(os.environ.get("PANEL_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "aligned")
SPECIES = ["arabiensis", "coluzzii", "coustani", "funestus", "gambiae", "stephensi"]

EXCLUDE = {
    ("coluzzii", "barcode80"),   # negative control (no-template)
    ("coustani", "barcode84"),   # negative control (no-template)
    ("coustani", "barcode95"),   # negative control (no-template)
    ("arabiensis", "barcode06"), # confirmed amplification/sequencing failure
}

rows = []
for sp in SPECIES:
    bams = sorted(glob.glob(os.path.join(BASE, sp, "*.sorted.bam")))
    for bam in bams:
        bc = os.path.basename(bam).split(".")[0]
        if (sp, bc) in EXCLUDE:
            continue
        rows.append({"species": sp, "barcode": bc, "bam": bam,
                     "sample_id": f"{sp}_{bc}"})

with open("manifest.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sample_id", "species", "barcode", "bam"])
    w.writeheader()
    w.writerows(rows)

from collections import Counter
c = Counter(r["species"] for r in rows)
print("Samples per species (post-exclusion):")
for sp in SPECIES:
    print(f"  {sp}: {c[sp]}")
print(f"Total: {len(rows)}")
