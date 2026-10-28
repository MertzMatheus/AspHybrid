#!/usr/bin/env python3
"""
de_analysis.py - DEG calls, functional enrichment, figures and tables for ONE
species, from the DESeq2 output of deseq2_strains.R.

DEG: padj < --padj AND |shrunken log2FC| >= --lfc   (up = higher with macrophage)

Enrichment (over-representation, per strain x direction x collection):
  one-sided hypergeometric test P(X >= k), BH correction within each
  strain x direction x collection. Universe = genes tested by DESeq2
  (padj not NA) that have >= 1 annotation in the collection; only those
  annotated DEGs are the test set. Gene sets with --min-set..--max-set genes
  in the universe (max not applied to COG / KEGG, whose broad categories are
  meaningful). Collections: GO_BP, GO_MF, GO_CC (propagated), KEGG_pathway, COG.

usage: de_analysis.py --species SP --de-dir 11_differential_expression/SP/deseq2
          --annotation 10_functional_annotation/SP/SP.gene_annotation.tsv
          --genesets 10_functional_annotation/SP/genesets --out 11_differential_expression/SP
          --fig-dir figures [--padj 0.05 --lfc 1 --enrich-padj 0.05 --min-set 5 --max-set 500]
"""
import argparse, csv, glob, math, os, re
from collections import defaultdict, Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
import numpy as np

csv.field_size_limit(10**9)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "axes.titleweight": "bold", "axes.titlesize": 10})
UP, DOWN, NS = "#C0392B", "#2E6DB4", "#C8C8C8"

ap = argparse.ArgumentParser()
ap.add_argument("--species", required=True); ap.add_argument("--de-dir", required=True)
ap.add_argument("--annotation"); ap.add_argument("--genesets")
ap.add_argument("--out", required=True); ap.add_argument("--fig-dir", required=True)
ap.add_argument("--padj", type=float, default=0.05); ap.add_argument("--lfc", type=float, default=1.0)
ap.add_argument("--enrich-padj", type=float, default=0.05)
ap.add_argument("--min-set", type=int, default=5); ap.add_argument("--max-set", type=int, default=500)
a = ap.parse_args()
SP = a.species; TITLE = SP.replace("_", ". ", 1)
os.makedirs(a.out, exist_ok=True); os.makedirs(a.fig_dir, exist_ok=True)
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
num = lambda x: float(x) if x not in ("", "NA", "NaN") else float("nan")


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.fig_dir, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write(path, rows, cols=None):
    if not rows:
        cols = cols or ["no_rows"]
    cols = cols or list(rows[0].keys())
    with open(path, "w", newline="") as fo:
        w = csv.DictWriter(fo, cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)


# ---------- inputs ----------
meta = rd(os.path.join(a.de_dir, "samples_used.tsv"))
smeta = {m["sample_old"]: m for m in meta}
strains = list(dict.fromkeys(m["strain"] for m in meta))
scol = {m["strain"]: m.get("strain_color") or "#555555" for m in meta}
de = {}
for st in strains:
    p = os.path.join(a.de_dir, f"DE_{st}.tsv")
    if os.path.exists(p):
        de[st] = rd(p)
strains = [s for s in strains if s in de]
ann, ann_cols = {}, []
if a.annotation and os.path.exists(a.annotation):
    r_ = csv.DictReader(open(a.annotation), delimiter="\t")
    ann_cols = [c for c in r_.fieldnames if c != "locus_tag"]
    for r in r_: ann[r["locus_tag"]] = r
label = lambda g: (ann.get(g, {}).get("eggNOG_preferred_name") or g)

# ---------- DEG calls + annotated tables ----------
summary, calls = [], defaultdict(dict)
for st in strains:
    rows = []
    for r in de[st]:
        p, l = num(r["padj"]), num(r["log2FC"])
        call = "up" if p < a.padj and l >= a.lfc else "down" if p < a.padj and l <= -a.lfc else "ns"
        if math.isnan(p): call = "not_tested"
        calls[st][r["locus_tag"]] = (call, l, p)
        rr = {"species": SP, "strain": st, "contrast": f"{st} Macrophage vs DMEN", **r, "DE": call}
        rr.update({c: ann.get(r["locus_tag"], {}).get(c, "") for c in ann_cols})
        rows.append(rr)
    rows.sort(key=lambda r: (r["DE"] not in ("up", "down"), num(r["padj"]) if not math.isnan(num(r["padj"])) else 9))
    cols = list(rows[0].keys())
    write(os.path.join(a.out, f"DE_{st}.annotated.tsv"), rows, cols)
    write(os.path.join(a.out, f"DEGs_{st}.tsv"), [r for r in rows if r["DE"] in ("up", "down")], cols)
    c = Counter(r["DE"] for r in rows)
    tested = sum(1 for r in rows if r["DE"] != "not_tested")
    summary.append({"species": SP, "strain": st, "contrast": "Macrophage vs DMEN", "genes_in_model": len(rows),
                    "genes_tested": tested, "up": c["up"], "down": c["down"], "DEGs": c["up"] + c["down"],
                    "percent_DEGs_of_tested": f"{100 * (c['up'] + c['down']) / tested:.2f}" if tested else "0",
                    "criteria": f"padj<{a.padj} & |log2FC_shrunk|>={a.lfc}"})
write(os.path.join(a.out, "DE_summary.tsv"), summary)
for s in summary:
    print(f"  {SP} {s['strain']:>6}: {s['up']:>5} up  {s['down']:>5} down  ({s['percent_DEGs_of_tested']}% of {s['genes_tested']} tested)")

# ---------- sharing between strains ----------
genes_all = sorted({g for st in strains for g in calls[st]})
memb = []
for g in genes_all:
    r = {"locus_tag": g}
    for st in strains:
        c, l, p = calls[st].get(g, ("absent", float("nan"), float("nan")))
        r[f"{st}_DE"] = c; r[f"{st}_log2FC"] = "" if math.isnan(l) else f"{l:.3f}"; r[f"{st}_padj"] = "" if math.isnan(p) else f"{p:.3g}"
    ups = [st for st in strains if r[f"{st}_DE"] == "up"]; dns = [st for st in strains if r[f"{st}_DE"] == "down"]
    if not ups and not dns: continue
    r["up_in"] = "+".join(ups); r["down_in"] = "+".join(dns)
    r["pattern"] = ("opposite" if ups and dns else
                    ("up in all strains" if len(ups) == len(strains) else f"up: {r['up_in']}") if ups else
                    ("down in all strains" if len(dns) == len(strains) else f"down: {r['down_in']}"))
    for c in ("product", "eggNOG_preferred_name", "COG_category", "subgenome"):
        if c in ann_cols: r[c] = ann.get(g, {}).get(c, "")
    memb.append(r)
write(os.path.join(a.out, "DEG_strain_membership.tsv"), memb)
pat = Counter()
for r in memb:
    for d in ("up", "down"):
        if r[f"{d}_in"]: pat[(d, r[f"{d}_in"])] += 1
write(os.path.join(a.out, "DEG_sharing.tsv"),
      [{"species": SP, "direction": d, "strains": s, "n_strains": s.count("+") + 1, "genes": n}
       for (d, s), n in sorted(pat.items(), key=lambda x: (x[0][0] != "up", -x[0][1].count("+"), x[0][1]))])


# ---------- enrichment ----------
def logC(n, k): return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hyper_sf(k, N, K, n):          # P(X >= k), X ~ Hypergeom(N, K, n)
    if k <= 0: return 1.0
    hi = min(K, n); lo = k
    if lo > hi: return 0.0
    den = logC(N, n)
    terms = [logC(K, i) + logC(N - K, n - i) - den for i in range(lo, hi + 1)]
    m = max(terms)
    return min(1.0, math.exp(m) * sum(math.exp(t - m) for t in terms))


def bh(ps):
    n = len(ps); order = sorted(range(n), key=lambda i: ps[i]); q = [0] * n; prev = 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        prev = min(prev, ps[i] * n / rank); q[i] = prev
    return q


COLLS = ["GO_BP", "GO_MF", "GO_CC", "KEGG_pathway", "COG"]
gs = {}
if a.genesets:
    for c in COLLS:
        f = os.path.join(a.genesets, f"{SP}.{c}.TERM2GENE.tsv")
        if not os.path.exists(f): continue
        t2g = defaultdict(set)
        for r in rd(f): t2g[r["term"]].add(r["gene"])
        nm = {r["term"]: r["name"] for r in rd(os.path.join(a.genesets, f"{SP}.{c}.TERM2NAME.tsv"))}
        gs[c] = (t2g, nm)
enr = []
for st in strains:
    tested = {g for g, (c, l, p) in calls[st].items() if c != "not_tested"}
    for d in ("up", "down"):
        degs = {g for g in tested if calls[st][g][0] == d}
        for c, (t2g, nm) in gs.items():
            annotated = set().union(*t2g.values()) & tested if t2g else set()
            N = len(annotated); sel = degs & annotated; n = len(sel)
            if n == 0 or N == 0: continue
            rows = []
            for t, genes in t2g.items():
                tg = genes & annotated; K = len(tg)
                if K < a.min_set or (c.startswith("GO") and K > a.max_set): continue
                hit = tg & sel; k = len(hit)
                if k == 0: continue
                rows.append({"species": SP, "strain": st, "direction": d, "collection": c, "term": t,
                             "name": nm.get(t, t), "DEGs_in_term": k, "DEGs_annotated": n, "term_size": K,
                             "universe": N, "gene_ratio": f"{k}/{n}", "fold_enrichment": f"{(k / n) / (K / N):.3f}",
                             "pvalue": hyper_sf(k, N, K, n),
                             "genes": ",".join(sorted(hit)),
                             "gene_names": ",".join(sorted({label(g) for g in hit}))})
            for r, q in zip(rows, bh([r["pvalue"] for r in rows])):
                r["padj"] = q
            enr += rows
for r in enr:
    r["significant"] = "yes" if r["padj"] < a.enrich_padj else "no"
enr.sort(key=lambda r: (r["strain"], r["direction"], r["collection"], r["padj"]))
ecols = ["species", "strain", "direction", "collection", "term", "name", "DEGs_in_term", "DEGs_annotated",
         "term_size", "universe", "gene_ratio", "fold_enrichment", "pvalue", "padj", "significant", "genes", "gene_names"]
fmt = lambda rows: [{**r, "pvalue": f"{r['pvalue']:.3g}", "padj": f"{r['padj']:.3g}"} for r in rows]
write(os.path.join(a.out, "enrichment_all_tests.tsv"), fmt(enr), ecols)
sig = [r for r in enr if r["significant"] == "yes"]
write(os.path.join(a.out, "enrichment_significant.tsv"), fmt(sig), ecols)
es = Counter((r["strain"], r["direction"], r["collection"]) for r in sig)
write(os.path.join(a.out, "enrichment_summary.tsv"),
      [{"species": SP, "strain": st, "direction": d, "collection": c, "significant_terms": es[(st, d, c)]}
       for st in strains for d in ("up", "down") for c in gs])
print(f"  enrichment: {len(sig)} significant terms (padj < {a.enrich_padj}) in {len(enr)} tests")

# ---------- figure 1: PCA, DEG bars, volcanoes, sharing ----------
vst = rd(os.path.join(a.de_dir, "vst.tsv"))
samples = [s for s in vst[0] if s != "locus_tag"]
X = np.array([[float(r[s]) for s in samples] for r in vst])
top = np.argsort(X.var(1))[::-1][:500]
Y = X[top].T; Y = Y - Y.mean(0)
U, S, _ = np.linalg.svd(Y, full_matrices=False); pcs = U * S; ve = 100 * S ** 2 / (S ** 2).sum()
with open(os.path.join(a.out, "PCA_coordinates.tsv"), "w") as fo:
    fo.write("sample\tsample_id\tstrain\ttreatment\t" + "\t".join(f"PC{i+1}_{ve[i]:.1f}pct" for i in range(min(5, len(ve)))) + "\n")
    for j, s in enumerate(samples):
        m = smeta[s]
        fo.write(f"{s}\t{m['sample_id']}\t{m['strain']}\t{m['treatment']}\t" + "\t".join(f"{pcs[j, i]:.4f}" for i in range(min(5, len(ve)))) + "\n")

nst = len(strains)
up_sets = [{g for g, v in calls[st].items() if v[0] == "up"} for st in strains]
dn_sets = [{g for g, v in calls[st].items() if v[0] == "down"} for st in strains]
if nst > 1:     # row 1: PCA, bars, Venn up, Venn down; row 2: one volcano per strain
    fig = plt.figure(figsize=(16, 9.4))
    gsp = fig.add_gridspec(2, 12, hspace=0.62, wspace=1.2)
    axes_top = [fig.add_subplot(gsp[0, i * 3:(i + 1) * 3]) for i in range(4)]
    vw = 12 // nst
    axes_vol = [fig.add_subplot(gsp[1, i * vw:(i + 1) * vw]) for i in range(nst)]
else:           # one strain: PCA, bars, volcano in one row (no sharing between strains)
    fig = plt.figure(figsize=(15, 5.2))
    gsp = fig.add_gridspec(1, 7, wspace=1.4)
    axes_top = [fig.add_subplot(gsp[0, 0:3]), fig.add_subplot(gsp[0, 3])]
    axes_vol = [fig.add_subplot(gsp[0, 4:7])]

ax = axes_top[0]
for j, s in enumerate(samples):
    m = smeta[s]; mac = m["treatment"] == "Macrophage"
    ax.scatter(pcs[j, 0], pcs[j, 1], s=55, marker="^" if mac else "o", color=scol[m["strain"]],
               edgecolor="black" if mac else scol[m["strain"]], lw=0.6 if mac else 0, zorder=3)
ax.axhline(0, color="#e5e5e5", lw=0.6, zorder=0); ax.axvline(0, color="#e5e5e5", lw=0.6, zorder=0)
ax.set_xlabel(f"PC1 ({ve[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({ve[1]:.1f}%)")
ax.set_title("PCA (vst, top 500 genes)")
h = [Line2D([], [], ls="", marker="s", color=scol[st], label=f"strain {st}") for st in strains] + \
    [Line2D([], [], ls="", marker="o", color="#777777", label="DMEN (fungus only)"),
     Line2D([], [], ls="", marker="^", color="#777777", markeredgecolor="black", label="+ macrophage")]
ax.legend(handles=h, frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)

ax = axes_top[1]
xs = np.arange(nst); up = [int(s["up"]) for s in summary]; dn = [int(s["down"]) for s in summary]
ax.bar(xs - 0.2, up, 0.4, color=UP, label="up"); ax.bar(xs + 0.2, dn, 0.4, color=DOWN, label="down")
for x, u, d in zip(xs, up, dn):
    ax.text(x - 0.2, u, f"{u:,}", ha="center", va="bottom", fontsize=7.5)
    ax.text(x + 0.2, d, f"{d:,}", ha="center", va="bottom", fontsize=7.5)
ax.set_xticks(xs, [f"strain {s}" for s in strains]); ax.set_ylabel("DEGs"); ax.set_xlim(-0.7, nst - 0.3)
ax.set_title("DEGs, macrophage vs DMEN" if nst > 1 else "DEGs"); ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper center")
ax.set_ylim(0, max(up + dn + [1]) * 1.3)


def venn(ax, sets, names, color, title):
    ax.set_title(title); ax.axis("off")
    if len(sets) == 2:
        cen = [(-0.55, 0), (0.55, 0)]; lab = [(-1.0, 0), (1.0, 0), (0, 0)]
        A, B = sets; vals = [len(A - B), len(B - A), len(A & B)]
        npos = [(-0.75, 1.3), (0.75, 1.3)]
    else:
        cen = [(-0.6, 0.35), (0.6, 0.35), (0, -0.6)]
        A, B, C = sets[:3]
        vals = [len(A - B - C), len(B - A - C), len(C - A - B), len((A & B) - C), len((A & C) - B), len((B & C) - A), len(A & B & C)]
        lab = [(-1.15, 0.65), (1.15, 0.65), (0, -1.25), (0, 0.85), (-0.75, -0.35), (0.75, -0.35), (0, 0.05)]
        npos = [(-1.1, 1.6), (1.1, 1.6), (0, -1.95)]
    for c in cen:
        ax.add_patch(Circle(c, 1.05, color=color, alpha=0.18, lw=0)); ax.add_patch(Circle(c, 1.05, fill=False, color=color, lw=1.2))
    for (x, y), v in zip(lab, vals):
        ax.text(x, y, f"{v:,}", ha="center", va="center", fontsize=9.5)
    for (x, y), n in zip(npos, names):
        ax.text(x, y, n, ha="center", va="center", fontsize=9, weight="bold")
    ax.set_xlim(-2.2, 2.2); ax.set_ylim(-2.2, 2.0); ax.set_aspect("equal")


if nst > 1:
    venn(axes_top[2], up_sets, [f"strain {s}" for s in strains], UP, "Up-regulated DEGs shared")
    venn(axes_top[3], dn_sets, [f"strain {s}" for s in strains], DOWN, "Down-regulated DEGs shared")
for ax, st in zip(axes_vol, strains):
    items = [(g, v) for g, v in calls[st].items() if v[0] != "not_tested"]
    gl = [g for g, _ in items]; l = np.array([v[1] for _, v in items]); p = np.array([v[2] for _, v in items])
    c = np.array([v[0] for _, v in items]); yv = -np.log10(np.clip(p, 1e-300, 1))
    for k, col in (("ns", NS), ("up", UP), ("down", DOWN)):
        m = c == k
        ax.scatter(l[m], yv[m], s=4 if k == "ns" else 7, color=col, alpha=0.6 if k == "ns" else 0.85, lw=0,
                   rasterized=True, label=f"{k} ({m.sum():,})")
    for x in (-a.lfc, a.lfc): ax.axvline(x, ls="--", lw=0.6, c="#777777")
    ax.axhline(-math.log10(a.padj), ls="--", lw=0.6, c="#777777")
    for j in [j for j in np.argsort(p) if c[j] in ("up", "down")][:8]:
        ax.annotate(label(gl[j]), (l[j], yv[j]), fontsize=6, xytext=(2, 2), textcoords="offset points", color="#333333")
    ax.set_xlabel("log2 fold change (macrophage / DMEN)"); ax.set_ylabel("-log10 adjusted p")
    ax.set_title(f"Volcano, strain {st}")
    ax.legend(frameon=False, fontsize=7, markerscale=2, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
fig.suptitle(f"{TITLE}: fungus + macrophage vs fungus alone (DMEN)", style="italic", fontsize=12, y=1.04 if nst == 1 else 1.0)
for i, axx in enumerate(axes_top + axes_vol):
    axx.text(-0.1 if nst > 1 or axx is not axes_top[1] else -0.35, 1.07, "ABCDEFGHIJ"[i], transform=axx.transAxes, fontsize=13, weight="bold")
save(fig, f"Fig_DE_{SP}")

# ---------- figure 2: enrichment dot plot ----------
show = [r for r in sig if r["collection"] in ("GO_BP", "KEGG_pathway")]
if show:
    pick = []
    for st in strains:
        for d in ("up", "down"):
            for c in ("GO_BP", "KEGG_pathway"):
                pick += [r for r in show if r["strain"] == st and r["direction"] == d and r["collection"] == c][:8]
    terms = list(dict.fromkeys((r["collection"], r["term"], r["name"]) for r in pick))
    cols_x = [(st, d) for st in strains for d in ("up", "down")]
    fig, ax = plt.subplots(figsize=(2.0 + 1.0 * len(cols_x), 1.2 + 0.24 * len(terms)))
    look = {(r["strain"], r["direction"], r["term"]): r for r in sig}
    for yi, (c, t, nm) in enumerate(terms):
        for xi, (st, d) in enumerate(cols_x):
            r = look.get((st, d, t))
            if r:
                k, n = map(int, r["gene_ratio"].split("/"))
                ax.scatter(xi, yi, s=20 + 900 * k / n, c=-math.log10(max(r["padj"], 1e-30)), cmap="viridis",
                           vmin=0, vmax=10, edgecolor="black", lw=0.3)
    ax.set_yticks(range(len(terms)), [f"{'GO' if c.startswith('GO') else 'KEGG'} | {nm[:60]}" for c, t, nm in terms], fontsize=7)
    ax.set_xticks(range(len(cols_x)), [f"{st}\n{d}" for st, d in cols_x], fontsize=8)
    ax.set_xlim(-0.6, len(cols_x) - 0.4); ax.set_ylim(len(terms) - 0.5, -0.5); ax.grid(alpha=0.25, lw=0.5)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, 10))
    plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, label="-log10 adjusted p")
    for fr in (0.05, 0.2, 0.4):
        ax.scatter([], [], s=20 + 900 * fr, c="#999999", label=f"{int(fr * 100)}% of DEGs")
    ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.18, 1), title="gene ratio", title_fontsize=7)
    ax.set_title(f"{TITLE}: enriched GO BP / KEGG terms (top 8 per strain x direction)", style="italic")
    save(fig, f"Fig_enrichment_{SP}")
else:
    print("  no significant GO BP / KEGG terms: enrichment figure skipped")

# ---------- README for the supplementary workbook ----------
with open(os.path.join(a.out, "README_tables.tsv"), "w") as fo:
    fo.write("item\tdescription\n")
    for k, v in [
        ("species", TITLE), ("design", "DESeq2 ~ group (strain_treatment); per strain: Macrophage vs DMEN (fungus alone, control)"),
        ("DEG", f"padj (BH) < {a.padj} and |log2FC| >= {a.lfc}; log2FC = shrunken (apeglm) estimate; log2FC_MLE = unshrunken"),
        ("up / down", "up = higher with macrophage; down = lower with macrophage"),
        ("not_tested", "gene removed by DESeq2 independent filtering / outlier detection (padj = NA)"),
        ("enrichment", "one-sided hypergeometric test; BH within strain x direction x collection; universe = tested genes annotated in the collection"),
        ("gene sets", f"eggNOG-mapper GO (propagated to ancestors), KEGG pathways, COG categories; set size >= {a.min_set} (GO also <= {a.max_set})"),
    ]:
        fo.write(f"{k}\t{v}\n")
