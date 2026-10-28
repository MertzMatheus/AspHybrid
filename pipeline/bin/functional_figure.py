#!/usr/bin/env python3
"""functional_figure.py --dir 10_functional_annotation --fig-dir figures --species A B C
Combined table + figure of the eggNOG annotation of every species done so far:
 A  % of genes with each annotation type
 B  COG functional categories (% of genes)
outputs: functional_annotation_summary.tsv, COG_summary.tsv, Fig_functional_annotation.{png,pdf}
"""
import argparse, csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "axes.titleweight": "bold", "axes.titlesize": 10})
COL = {"A_fumigatus": "#0072B2", "A_spinulosporus": "#E69F00", "A_latus": "#009E73"}
ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True); ap.add_argument("--fig-dir", required=True)
ap.add_argument("--species", nargs="+", required=True)
a = ap.parse_args()
sps = [s for s in a.species if os.path.exists(os.path.join(a.dir, s, f"{s}.annotation_coverage.tsv"))]
if not sps:
    raise SystemExit("no species annotated yet")
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
cov = {s: rd(os.path.join(a.dir, s, f"{s}.annotation_coverage.tsv")) for s in sps}
cog = {s: rd(os.path.join(a.dir, s, f"{s}.COG_counts.tsv")) for s in sps}

with open(os.path.join(a.dir, "functional_annotation_summary.tsv"), "w") as fo:
    keys = [r["annotation"] for r in cov[sps[0]]]
    fo.write("annotation\t" + "\t".join(f"{s}_genes\t{s}_percent" for s in sps) + "\n")
    for i, k in enumerate(keys):
        fo.write(k + "\t" + "\t".join(f"{cov[s][i]['genes']}\t{cov[s][i]['percent_of_genes']}" for s in sps) + "\n")
with open(os.path.join(a.dir, "COG_summary.tsv"), "w") as fo:
    fo.write("COG\tdescription\t" + "\t".join(f"{s}_genes\t{s}_percent" for s in sps) + "\n")
    for i, r in enumerate(cog[sps[0]]):
        fo.write(f"{r['COG']}\t{r['description']}\t" + "\t".join(f"{cog[s][i]['genes']}\t{cog[s][i]['percent_of_genes']}" for s in sps) + "\n")

show = ["eggNOG hit", "named product (not hypothetical)", "COG category", "GO (any)", "GO biological process",
        "EC number", "KEGG KO", "KEGG pathway", "Pfam domain", "CAZy"]
fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1, 1.5]})
ax = axes[0]; h = 0.8 / len(sps); y = np.arange(len(show))
for j, s in enumerate(sps):
    d = {r["annotation"]: float(r["percent_of_genes"]) for r in cov[s]}
    n = int(cov[s][0]["genes"])
    ax.barh(y + j * h, [d.get(k, 0) for k in show], height=h, color=COL.get(s, "#777777"),
            label=f"{s.replace('_', '. ', 1)} (n = {n:,})")
ax.set_yticks(y + h * (len(sps) - 1) / 2, show); ax.invert_yaxis(); ax.set_xlim(0, 100)
ax.set_xlabel("% of protein-coding genes"); ax.set_title("eggNOG-mapper annotation")
ax.legend(frameon=False, loc="lower right", prop={"style": "italic", "size": 7.5})
ax = axes[1]
pct = {s: {r["COG"]: float(r["percent_of_genes"]) for r in cog[s]} for s in sps}
cats = [r["COG"] for r in cog[sps[0]] if any(pct[s][r["COG"]] > 0 for s in sps)]
desc = {r["COG"]: r["description"] for r in cog[sps[0]]}
x = np.arange(len(cats)); w = 0.8 / len(sps)
for j, s in enumerate(sps):
    ax.barh(x + j * w, [pct[s][c] for c in cats], height=w, color=COL.get(s, "#777777"))
ax.set_yticks(x + w * (len(sps) - 1) / 2, [f"{c}  {desc[c][:55]}" for c in cats], fontsize=7)
ax.invert_yaxis(); ax.set_xlabel("% of genes"); ax.set_title("COG functional categories")
for i, axx in enumerate(axes):
    axx.text(-0.02 if i == 0 else -0.62, 1.04, "AB"[i], transform=axx.transAxes, fontsize=13, weight="bold")
plt.tight_layout()
os.makedirs(a.fig_dir, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(a.fig_dir, f"Fig_functional_annotation.{ext}"), dpi=300, bbox_inches="tight")
