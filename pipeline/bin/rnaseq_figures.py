#!/usr/bin/env python3
"""
rnaseq_figures.py - publication QC figures for one species' RNA-seq.

Inputs: outputs of make_count_matrix.py
  --prefix  Results/RNAseq/<sp>/counts/<sp>   (uses .counts.tsv, .mapping_stats.tsv,
                                              .sample_metadata.tsv, .counts.annotated.tsv)
  --species name used in titles
  --outdir  figure directory

Figures (PNG 300 dpi + PDF, vector):
  1_read_fate          % of reads per sample: filtered, unaligned, multi-mapped,
                       no feature, ambiguous, assigned (stacked)
  2_library_size       assigned reads (M) and genes detected (>=10 reads)
  3_expression_density log2(CPM+1) distribution per sample
  4_PCA                PCA of log2(CPM+1), top N most variable expressed genes,
                       % variance explained on axes, coloured by condition
  5_sample_correlation Pearson correlation heatmap, hierarchical clustering
                       (average linkage, 1 - r)
  6_subgenome_share    A. latus only: % of assigned reads from S and Q genes
  QC_panel             all of the above in one multi-panel figure
Also writes PCA coordinates and the correlation matrix as TSV.

Statistics notes:
  * Gene filter: CPM >= 1 in at least min(3, smallest group) samples
    (edgeR filterByExpr-like rule), applied before PCA/correlation.
  * log2(CPM+1) is used here for transparency and speed; DESeq2 vst is
    used for the differential-expression PCA in deseq2_analysis.R (step 13).
"""
import argparse
import csv
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none", "axes.titleweight": "bold"})
PAL = ["#3b6ea5", "#c8553d", "#4b9b6e", "#8a6fb3", "#d4a72c", "#5aa5b8", "#a35d84", "#6b6b6b"]


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def save(fig, out, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top", type=int, default=500)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    title = a.species.replace("_", ". ", 1)

    rows = read_tsv(a.prefix + ".counts.tsv")
    samples = [k for k in rows[0] if k != "locus_tag"]
    genes = [r["locus_tag"] for r in rows]
    C = np.array([[float(r[s]) for s in samples] for r in rows])
    meta = {r["sample"]: r for r in read_tsv(a.prefix + ".sample_metadata.tsv")}
    # colour = group (strain x condition); in the PCA colour = strain, marker = condition
    cond = [meta.get(s, {}).get("group") or meta.get(s, {}).get("condition", "NA") for s in samples]
    strain = [meta.get(s, {}).get("strain", "NA") for s in samples]
    treat = [meta.get(s, {}).get("condition", "NA") for s in samples]
    levels = list(dict.fromkeys(cond))
    ccol = {c: PAL[i % len(PAL)] for i, c in enumerate(levels)}
    stats = {r["sample"]: r for r in read_tsv(a.prefix + ".mapping_stats.tsv")} \
        if os.path.exists(a.prefix + ".mapping_stats.tsv") else {}

    lib = C.sum(0)
    cpm = C / np.where(lib == 0, 1, lib) * 1e6
    logc = np.log2(cpm + 1)
    min_n = max(2, min(3, min(cond.count(c) for c in levels)))
    keep = (cpm >= 1).sum(1) >= min_n
    X = logc[keep]

    panels = {}

    # 1 read fate
    def fate(s):
        st = stats.get(s, {})
        raw, clean = num(st.get("raw_reads")), num(st.get("clean_reads"))
        pe = st.get("layout") == "PE"
        unit = 2 if pe else 1
        cleanu = clean / unit
        uniq, multi = num(st.get("aligned_unique")), num(st.get("aligned_multi"))
        assigned = num(st.get("assigned_to_genes"))
        nofeat, amb = num(st.get("unassigned_no_feature")), num(st.get("unassigned_ambiguous"))
        if np.isnan(cleanu) or cleanu == 0:
            return None
        filt = (raw - clean) / unit
        unal = cleanu - uniq - multi
        other = max(uniq - assigned - nofeat - amb, 0)
        parts = [filt, unal, multi, nofeat, amb, other, assigned]
        tot = raw / unit
        return [100 * p / tot for p in parts]
    labels = ["filtered (fastp)", "unaligned", "multi-mapped", "no feature", "ambiguous", "other unassigned", "assigned to genes"]
    fcol = ["#bdbdbd", "#8c8c8c", "#f0a35e", "#c7dbe8", "#7aa6c2", "#e6e6e6", "#2f4b7c"]
    F = [fate(s) for s in samples]

    def draw_fate(ax):
        if not any(F):
            ax.text(0.5, 0.5, "mapping stats not available", ha="center"); ax.axis("off"); return
        left = np.zeros(len(samples))
        for k, lab in enumerate(labels):
            v = np.array([f[k] if f else 0 for f in F])
            ax.barh(range(len(samples)), v, left=left, color=fcol[k], label=lab, edgecolor="white", lw=0.4)
            left += v
        ax.set_yticks(range(len(samples)), samples, fontsize=7)
        ax.invert_yaxis(); ax.set_xlim(0, 100); ax.set_xlabel("% of sequenced reads (pairs)")
        ax.legend(ncol=4, fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12))
        ax.set_title("Read fate")
    panels["1_read_fate"] = draw_fate

    # 2 library size / detection
    def draw_lib(ax):
        x = np.arange(len(samples))
        ax.bar(x, lib / 1e6, color=[ccol[c] for c in cond])
        ax.set_xticks(x, samples, rotation=90, fontsize=7)
        ax.set_ylabel("Assigned reads (millions)")
        ax2 = ax.twinx()
        det = (C >= 10).sum(0)
        ax2.plot(x, det, "k.", ms=6)
        ax2.set_ylabel("Genes with >=10 reads (dots)")
        ax2.set_ylim(0, len(genes) * 1.05)
        ax2.spines["top"].set_visible(False)
        ax.set_title("Library size and detected genes")
    panels["2_library_size"] = draw_lib

    # 3 density
    def draw_density(ax):
        grid = np.linspace(0, max(1, logc.max()), 200)
        for j, s in enumerate(samples):
            v = logc[:, j][logc[:, j] > 0]
            try:
                from scipy.stats import gaussian_kde
                y = gaussian_kde(v)(grid)
            except Exception:
                h, e = np.histogram(v, bins=grid, density=True)
                y = np.append(h, 0)
            ax.plot(grid, y, color=ccol[cond[j]], lw=0.9, alpha=0.8)
        ax.set_xlabel("log2(CPM + 1), expressed genes"); ax.set_ylabel("density")
        ax.set_title("Expression distributions")
        for c in levels:
            ax.plot([], [], color=ccol[c], label=c)
        ax.legend(frameon=False, fontsize=7)
    panels["3_expression_density"] = draw_density

    # 4 PCA
    pca_out = None
    if X.shape[0] >= 3 and len(samples) >= 3:
        v = X.var(1)
        top = np.argsort(v)[::-1][:min(a.top, X.shape[0])]
        Y = X[top].T
        Y = Y - Y.mean(0)
        U, S, Vt = np.linalg.svd(Y, full_matrices=False)
        pcs = U * S
        ve = 100 * S ** 2 / (S ** 2).sum()
        pca_out = (pcs, ve, len(top))
        with open(os.path.join(a.outdir, "PCA_coordinates.tsv"), "w") as fo:
            fo.write("sample\tgroup\t" + "\t".join(f"PC{i+1} ({ve[i]:.1f}%)" for i in range(min(5, len(ve)))) + "\n")
            for j, s in enumerate(samples):
                fo.write(f"{s}\t{cond[j]}\t" + "\t".join(f"{pcs[j, i]:.4f}" for i in range(min(5, len(ve)))) + "\n")

    def draw_pca(ax):
        if not pca_out:
            ax.text(0.5, 0.5, "PCA needs >=3 samples", ha="center"); ax.axis("off"); return
        pcs, ve, n = pca_out
        slev, tlev = list(dict.fromkeys(strain)), sorted(set(treat), key=lambda x: (x != "control", x))
        scol = {c: PAL[i % len(PAL)] for i, c in enumerate(slev)}
        mk = {t: m for t, m in zip(tlev, "osD^v")}
        for j in range(len(samples)):
            ctl = treat[j] == "control"
            ax.scatter(pcs[j, 0], pcs[j, 1] if pcs.shape[1] > 1 else 0, s=48, marker=mk[treat[j]],
                       facecolor="white" if ctl else scol[strain[j]], edgecolor=scol[strain[j]], lw=1.3, zorder=3)
        from matplotlib.lines import Line2D
        h = [Line2D([], [], ls="", marker="o", color=scol[s], label=f"strain {s}") for s in slev] + \
            [Line2D([], [], ls="", marker=mk[t], markerfacecolor="white" if t == "control" else "#555555",
                    markeredgecolor="#555555", label=("control (no macrophage)" if t == "control" else t)) for t in tlev]
        for j, s in enumerate(samples):
            ax.annotate(s, (pcs[j, 0], pcs[j, 1] if pcs.shape[1] > 1 else 0), fontsize=6,
                        xytext=(3, 3), textcoords="offset points", color="#555555")
        ax.set_xlabel(f"PC1 ({ve[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({ve[1]:.1f}%)" if len(ve) > 1 else "PC2")
        ax.axhline(0, color="#dddddd", lw=0.6, zorder=0); ax.axvline(0, color="#dddddd", lw=0.6, zorder=0)
        ax.legend(handles=h, frameon=False, fontsize=7)
        ax.set_title(f"PCA, top {n} variable genes")
    panels["4_PCA"] = draw_pca

    # 5 correlation heatmap with clustering
    R = np.corrcoef(X.T) if X.shape[0] > 1 else np.eye(len(samples))
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
        D = np.clip(1 - R, 0, None)
        np.fill_diagonal(D, 0)
        order = leaves_list(linkage(squareform(D, checks=False), "average")) if len(samples) > 2 else np.arange(len(samples))
    except ImportError:
        order = np.argsort(cond, kind="stable")
    Ro = R[np.ix_(order, order)]
    so = [samples[i] for i in order]
    with open(os.path.join(a.outdir, "sample_correlation.tsv"), "w") as fo:
        fo.write("sample\t" + "\t".join(samples) + "\n")
        for i, s in enumerate(samples):
            fo.write(s + "\t" + "\t".join(f"{x:.4f}" for x in R[i]) + "\n")

    def draw_corr(ax):
        lo = max(0.0, np.floor(np.nanmin(Ro) * 100) / 100)
        im = ax.imshow(Ro, cmap="viridis", vmin=lo, vmax=1)
        ax.set_xticks(range(len(so)), so, rotation=90, fontsize=6.5)
        ax.set_yticks(range(len(so)), so, fontsize=6.5)
        for i, j in enumerate(order):     # condition strip above the columns
            ax.add_patch(plt.Rectangle((i - 0.5, -0.9), 1, 0.3, color=ccol[cond[j]], clip_on=False))
        if len(so) <= 16:
            for i in range(len(so)):
                for k in range(len(so)):
                    ax.text(k, i, f"{Ro[i, k]:.2f}", ha="center", va="center", fontsize=5,
                            color="white" if Ro[i, k] < (lo + 1) / 2 else "black")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Pearson r, log2(CPM+1)")
        ax.set_title("Sample correlation (clustered)", pad=30)
        ax.spines[:].set_visible(False)
    panels["5_sample_correlation"] = draw_corr

    # 6 subgenome share (hybrid)
    sub = {}
    ann = a.prefix + ".counts.annotated.tsv"
    if os.path.exists(ann):
        for r in read_tsv(ann):
            if "subgenome" in r:
                sub[r["locus_tag"]] = r["subgenome"]
    if sub and {"S", "Q"} <= set(sub.values()):
        labs = ["S", "Q", "ambiguous", "unassigned"]
        share = defaultdict(lambda: np.zeros(len(samples)))
        for i, g in enumerate(genes):
            share[sub.get(g, "unassigned")] += C[i]
        tot = sum(share[l] for l in labs)
        n_genes = {l: sum(1 for g in genes if sub.get(g) == l) for l in labs}
        gfrac = {l: 100 * n_genes[l] / len(genes) for l in labs}

        def draw_sub(ax):
            left = np.zeros(len(samples))
            cols = {"S": "#c8553d", "Q": "#4b9b6e", "ambiguous": "#bdbdbd", "unassigned": "#e6e6e6"}
            for l in labs:
                v = 100 * share[l] / np.where(tot == 0, 1, tot)
                ax.barh(range(len(samples)), v, left=left, color=cols[l],
                        label=f"{l} ({gfrac[l]:.1f}% of genes)", edgecolor="white", lw=0.4)
                left += v
            ax.axvline(gfrac["S"], color="k", ls="--", lw=0.8)
            ax.set_yticks(range(len(samples)), samples, fontsize=7); ax.invert_yaxis()
            ax.set_xlim(0, 100); ax.set_xlabel("% of assigned reads")
            ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.15))
            ax.set_title("Subgenome expression share (dashed = S share of genes)")
        panels["6_subgenome_share"] = draw_sub
        with open(os.path.join(a.outdir, "subgenome_read_share.tsv"), "w") as fo:
            fo.write("sample\tcondition\t" + "\t".join(f"pct_reads_{l}" for l in labs) + "\n")
            for j, s in enumerate(samples):
                fo.write(f"{s}\t{cond[j]}\t" + "\t".join(f"{100 * share[l][j] / tot[j]:.2f}" if tot[j] else "NA" for l in labs) + "\n")

    sizes = {"1_read_fate": (7, 0.28 * len(samples) + 1.6), "2_library_size": (max(5, 0.35 * len(samples) + 2), 3.6),
             "3_expression_density": (4.5, 3.4), "4_PCA": (4.8, 4.2),
             "5_sample_correlation": (max(4.5, 0.3 * len(samples) + 2.5), max(4, 0.3 * len(samples) + 2)),
             "6_subgenome_share": (7, 0.28 * len(samples) + 1.8)}
    for name, fn in panels.items():
        fig, ax = plt.subplots(figsize=sizes[name])
        fn(ax)
        fig.suptitle(title, style="italic", fontsize=9, x=0.01, ha="left")
        save(fig, a.outdir, name)

    keys = list(panels)
    n = len(keys)
    ncol = 2
    nrow = (n + 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(13, 4.6 * nrow))
    for k, ax in zip(keys, axes.flat):
        panels[k](ax)
        ax.text(-0.12, 1.06, chr(65 + keys.index(k)), transform=ax.transAxes, fontsize=13, fontweight="bold")
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    fig.suptitle(f"{title} RNA-seq quality control", style="italic", fontsize=12)
    fig.tight_layout()
    save(fig, a.outdir, "QC_panel")
    print(f"{a.species}: {len(genes)} genes, {int(keep.sum())} expressed (CPM>=1 in >={min_n} samples); figures in {a.outdir}")


if __name__ == "__main__":
    main()
