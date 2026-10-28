#!/usr/bin/env python3
"""subgenome_figure.py --dir 05_latus_subgenomes --fig-dir figures
Summary table and figure of the A. latus S/Q subgenome assignment:
 A  histogram of identity(S) - identity(Q) per 5-kb window (bimodal = clean phasing)
 B  genome fraction S / Q / ambiguous / unassigned
 C  BUSCO of each subgenome (each should be ~complete and single-copy)
 D  S/Q blocks along the 20 largest contigs
"""
import argparse, csv, glob, os, re
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "axes.titleweight": "bold", "axes.titlesize": 10})
COL = {"S": "#c8553d", "Q": "#4b9b6e", "ambiguous": "#bdbdbd", "unassigned": "#ececec"}


def busco(d):
    f = glob.glob(os.path.join(d, "short_summary*.txt"))
    if not f:
        return None
    m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%", open(f[0]).read())
    return dict(zip("CSDFM", map(float, m.groups()))) if m else None


ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True); ap.add_argument("--fig-dir", required=True)
a = ap.parse_args()
win = list(csv.DictReader(open(os.path.join(a.dir, "A_latus.windows.tsv")), delimiter="\t"))
bed = [l.split("\t") for l in open(os.path.join(a.dir, "A_latus.blocks.bed"))]
diff = [100 * (float(r["ident_S"]) - float(r["ident_Q"])) for r in win if float(r["ident_S"]) or float(r["ident_Q"])]
bp = defaultdict(int)
for r in win:
    bp[r["label"]] += int(r["end"]) - int(r["start"])
tot = sum(bp.values())
bS, bQ = busco(os.path.join(a.dir, "busco_S")), busco(os.path.join(a.dir, "busco_Q"))

rows = [("windows_5kb", len(win))] + [(f"{k}_Mb", round(bp[k] / 1e6, 2)) for k in COL] + \
       [(f"{k}_percent", round(100 * bp[k] / tot, 2)) for k in COL] + \
       [("median_identity_gap_S_windows_percent", round(float(np.median([d for d in diff if d >= 2])), 2) if diff else "NA"),
        ("median_identity_gap_Q_windows_percent", round(float(np.median([-d for d in diff if d <= -2])), 2) if diff else "NA")]
for lab, b in (("S", bS), ("Q", bQ)):
    if b:
        rows += [(f"BUSCO_{lab}_{k}", b[k]) for k in "CSDFM"]
with open(os.path.join(a.dir, "subgenome_summary.tsv"), "w") as fo:
    fo.write("metric\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in rows))

fig = plt.figure(figsize=(11, 7))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.1], hspace=0.55, wspace=0.4)
ax = fig.add_subplot(gs[0, 0])
ax.hist(diff, bins=np.arange(-15, 15.25, 0.25), color="#555555")
ax.axvspan(2, 15, color=COL["S"], alpha=0.12); ax.axvspan(-15, -2, color=COL["Q"], alpha=0.12)
for x in (-2, 2):
    ax.axvline(x, ls="--", c="k", lw=0.7)
ax.set_xlim(-15, 15)
ax.set_xlabel("identity to A. spinulosporus − identity to\nA. quadrilineatus (%)"); ax.set_ylabel("5-kb windows")
ax.set_title("Parental assignment")
ax.text(0.97, 0.95, "S", transform=ax.transAxes, ha="right", va="top", color=COL["S"], fontsize=12, weight="bold")
ax.text(0.03, 0.95, "Q", transform=ax.transAxes, ha="left", va="top", color=COL["Q"], fontsize=12, weight="bold")

ax = fig.add_subplot(gs[0, 1])
ks = [k for k in COL if bp[k]]
ax.bar(range(len(ks)), [bp[k] / 1e6 for k in ks], color=[COL[k] for k in ks], edgecolor="#999999", lw=0.4)
for i, k in enumerate(ks):
    ax.text(i, bp[k] / 1e6, f"{100 * bp[k] / tot:.1f}%", ha="center", va="bottom", fontsize=8)
ax.set_xticks(range(len(ks)), ks); ax.set_ylabel("Mb of A. latus genome"); ax.set_title("Genome fraction")

ax = fig.add_subplot(gs[0, 2])
if bS and bQ:
    left = np.zeros(2)
    for k, c, l in (("S", "#2f4b7c", "single"), ("D", "#7aa6c2", "duplicated"), ("F", "#f0a35e", "fragmented"), ("M", "#d1495b", "missing")):
        v = np.array([bS[k], bQ[k]])
        ax.barh([0, 1], v, left=left, color=c, label=l, edgecolor="white", lw=0.5)
        left += v
    ax.set_yticks([0, 1], ["subgenome S", "subgenome Q"]); ax.invert_yaxis(); ax.set_xlim(0, 100)
    ax.set_xlabel("% BUSCOs"); ax.legend(fontsize=7, frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.25))
ax.set_title("Subgenome completeness")

ax = fig.add_subplot(gs[1, :])
clen = defaultdict(int)
for f in bed:
    clen[f[0]] = max(clen[f[0]], int(f[2]))
top = sorted(clen, key=clen.get, reverse=True)[:20]
for i, c in enumerate(top):
    for f in bed:
        if f[0] == c:
            ax.barh(i, (int(f[2]) - int(f[1])) / 1e6, left=int(f[1]) / 1e6, color=COL.get(f[3].strip(), "#ececec"), height=0.75)
ax.set_yticks(range(len(top)), top, fontsize=6.5); ax.invert_yaxis()
ax.set_xlabel("Position (Mb)"); ax.set_title("Subgenome blocks along the 20 largest contigs")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=COL[k], label=k) for k in ("S", "Q", "ambiguous")], fontsize=7, frameon=False, loc="lower right")
for i, axx in enumerate(fig.axes):
    axx.text(-0.15 if i < 3 else -0.06, 1.08, "ABCD"[i], transform=axx.transAxes, fontsize=13, weight="bold")
os.makedirs(a.fig_dir, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(a.fig_dir, f"Fig_latus_subgenomes.{ext}"), dpi=300, bbox_inches="tight")
print(f"S {100*bp['S']/tot:.1f}%  Q {100*bp['Q']/tot:.1f}%  ambiguous {100*bp['ambiguous']/tot:.1f}%")
