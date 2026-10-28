#!/usr/bin/env python3
"""
ortholog_summary.py - A. fumigatus-anchored ortholog tables + statistics + figure.

For every A. fumigatus gene (official AFUA_... or de novo AFUM_...): its orthologue(s) in A. spinulosporus
and in A. latus (one proteome), with the relation type. A. latus is a hybrid,
so the expected relation for most genes is 1:2 (one A. fumigatus gene, two
A. latus homeologs, one per subgenome). If the A. latus id_map carries a
subgenome label (step 05/06), each A. latus orthologue is tagged S or Q.

Relation (A. fumigatus genes : orthologues in the other species):
  1:1  1:2  1:>2  many:1  many:many  none

usage: ortholog_summary.py --tables 07_orthology/tables --out 07_orthology --fig-dir figures [--id-map A_latus.id_map.tsv]
outputs: Table_orthologs_Afumigatus.tsv, Table_pairs_Afumigatus_Alatus.tsv,
         Table_pairs_Afumigatus_Aspinulosporus.tsv, ortholog_statistics.tsv, Fig_orthologs.{png,pdf}
"""
import argparse, csv, os
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

csv.field_size_limit(10**9)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "axes.titleweight": "bold", "axes.titlesize": 10})
ap = argparse.ArgumentParser()
ap.add_argument("--tables", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--fig-dir", required=True); ap.add_argument("--id-map")
a = ap.parse_args()

sub = {}
if a.id_map and os.path.exists(a.id_map):
    for r in csv.DictReader(open(a.id_map), delimiter="\t"):
        s = (r.get("subgenome") or "").strip()
        if s and s != "NA":
            sub[r["locus_tag"]] = s
has_sub = bool(sub)

rows = list(csv.DictReader(open(os.path.join(a.tables, "correspondence_A_fumigatus.tsv")), delimiter="\t"))
REL = ["1:1", "1:2", "1:>2", "many:1", "many:many", "none"]


def rel(genes, orth_rel):
    """refine OrthoFinder's 1:many into 1:2 / 1:>2"""
    if not genes:
        return "none"
    if orth_rel == "1:many":
        return "1:2" if len(genes) == 2 else "1:>2"
    if orth_rel == "many:1":
        return "many:1"
    return orth_rel


out, pairs_l, pairs_s = [], [], []
c_as, c_al, n_lat, c_sub = Counter(), Counter(), Counter(), Counter()
for r in rows:
    af = r["locus_tag"]
    as_g = [g for g in (r.get("A_spinulosporus_orthologs") or "").split(",") if g]
    al_g = [g for g in (r.get("A_latus_orthologs") or "").split(",") if g]
    ra = rel(as_g, r.get("A_spinulosporus_relation", "none"))
    rl = rel(al_g, r.get("A_latus_relation", "none"))
    d = {"A_fumigatus": af, "orthogroup": r["orthogroup"],
         "A_spinulosporus": ",".join(as_g), "A_spinulosporus_relation": ra,
         "A_latus": ",".join(al_g), "A_latus_n": len(al_g), "A_latus_relation": rl}
    if has_sub:
        labs = [sub.get(g, "NA") for g in al_g]
        d["A_latus_subgenome"] = ",".join(labs)
        d["A_latus_S"] = ",".join(g for g, l in zip(al_g, labs) if l == "S")
        d["A_latus_Q"] = ",".join(g for g, l in zip(al_g, labs) if l == "Q")
        k = ("S+Q" if "S" in labs and "Q" in labs else "S only" if "S" in labs else
             "Q only" if "Q" in labs else "unassigned only" if al_g else "absent")
        c_sub[k] += 1
    d["one_to_one_all"] = "yes" if ra == "1:1" and rl == "1:1" else "no"
    d["A_fumigatus_1to1_As_and_1to2_Al"] = "yes" if ra == "1:1" and rl == "1:2" else "no"
    out.append(d)
    c_as[ra] += 1; c_al[rl] += 1
    n_lat[min(len(al_g), 3)] += 1
    for g in al_g:
        pairs_l.append({"A_fumigatus": af, "A_latus": g, "subgenome": sub.get(g, "NA"),
                        "relation": rl, "orthogroup": r["orthogroup"]})
    for g in as_g:
        pairs_s.append({"A_fumigatus": af, "A_spinulosporus": g, "relation": ra, "orthogroup": r["orthogroup"]})

n = len(out)


def write(name, rr):
    with open(os.path.join(a.out, name), "w", newline="") as fo:
        w = csv.DictWriter(fo, list(rr[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rr)


write("Table_orthologs_Afumigatus.tsv", out)
if pairs_l:
    write("Table_pairs_Afumigatus_Alatus.tsv", pairs_l)
if pairs_s:
    write("Table_pairs_Afumigatus_Aspinulosporus.tsv", pairs_s)

with open(os.path.join(a.out, "ortholog_statistics.tsv"), "w") as fo:
    fo.write("comparison\tcategory\tA_fumigatus_genes\tpercent\n")
    for name, c in (("A_fumigatus -> A_spinulosporus", c_as), ("A_fumigatus -> A_latus", c_al)):
        for k in REL:
            fo.write(f"{name}\t{k}\t{c[k]}\t{100*c[k]/n:.2f}\n")
    for k, lab in ((0, "0"), (1, "1"), (2, "2"), (3, ">=3")):
        fo.write(f"A_fumigatus -> A_latus\tA. latus orthologues = {lab}\t{n_lat[k]}\t{100*n_lat[k]/n:.2f}\n")
    for k, v in sorted(c_sub.items(), key=lambda x: -x[1]):
        fo.write(f"A_fumigatus -> A_latus subgenomes\t{k}\t{v}\t{100*v/n:.2f}\n")
    k1 = sum(1 for d in out if d["A_fumigatus_1to1_As_and_1to2_Al"] == "yes")
    fo.write(f"A_fumigatus -> all\t1:1 with As and 1:2 with Al (expected hybrid pattern)\t{k1}\t{100*k1/n:.2f}\n")
    fo.write(f"A_fumigatus\ttotal genes\t{n}\t100\n")

# ---------- figure ----------
fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw={"width_ratios": [1.4, 1]})
ax = axes[0]
pal = {"1:1": "#2f4b7c", "1:2": "#4b9b6e", "1:>2": "#9fd3b0", "many:1": "#a5c8dd", "many:many": "#f0a35e", "none": "#d9d9d9"}
left = np.zeros(2)
for k in REL:
    v = np.array([100 * c[k] / n for c in (c_as, c_al)])
    ax.barh(range(2), v, left=left, color=pal[k], label=k, edgecolor="white", lw=0.5)
    for i, (l, x) in enumerate(zip(left, v)):
        if x >= 5:
            ax.text(l + x / 2, i, f"{x:.0f}%", ha="center", va="center", fontsize=7.5,
                    color="white" if k in ("1:1", "1:2") else "black")
    left += v
ax.set_yticks(range(2), ["A. spinulosporus", "A. latus"], style="italic"); ax.invert_yaxis(); ax.set_xlim(0, 100)
ax.set_xlabel(f"% of A. fumigatus genes (n = {n:,})")
ax.set_title("Orthology relation to A. fumigatus genes")
ax.legend(ncol=6, fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.25))
ax = axes[1]
labs = ["0", "1", "2", "≥3"]
v = [100 * n_lat[k] / n for k in range(4)]
ax.bar(range(4), v, color=["#d9d9d9", "#2f4b7c", "#4b9b6e", "#9fd3b0"])
for i, x in enumerate(v):
    ax.text(i, x, f"{x:.1f}%", ha="center", va="bottom", fontsize=7.5)
ax.set_xticks(range(4), labs); ax.set_xlabel("A. latus orthologues per A. fumigatus gene")
ax.set_ylabel("% of A. fumigatus genes"); ax.set_title("Copies in the hybrid A. latus")
for i, axx in enumerate(axes):
    axx.text(-0.12, 1.08, "AB"[i], transform=axx.transAxes, fontsize=13, weight="bold")
os.makedirs(a.fig_dir, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(a.fig_dir, f"Fig_orthologs.{ext}"), dpi=300, bbox_inches="tight")
print(f"A. fumigatus genes: {n}; with A. spinulosporus orthologue {100*(n-c_as['none'])/n:.1f}% "
      f"(1:1 {100*c_as['1:1']/n:.1f}%); with A. latus orthologue {100*(n-c_al['none'])/n:.1f}% "
      f"(1:2 {100*c_al['1:2']/n:.1f}%, 1:1 {100*c_al['1:1']/n:.1f}%)")
