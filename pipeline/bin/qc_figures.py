#!/usr/bin/env python3
"""
qc_figures.py - cross-species QC tables and publication figures.

  qc_figures.py genome     --qc-dir 01_genome_qc --fig-dir figures --species ...
      -> 01_genome_qc/genome_qc_summary.tsv, figures/Fig_genome_qc.{png,pdf}
  qc_figures.py annotation --ref-dir 06_reference_files --fig-dir figures --species ...
      -> 06_reference_files/annotation_summary_all.tsv, figures/Fig_annotation_qc.{png,pdf}

Species that are not finished yet are skipped, so the figures can be rebuilt
after each species.
"""
import argparse
import csv
import glob
import os
import re
import statistics as st
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none", "axes.titleweight": "bold",
                     "axes.titlesize": 10})
SPCOL = {"A_fumigatus": "#3b6ea5", "A_spinulosporus": "#c8553d", "A_latus": "#4b9b6e",
         "A_quadrilineatus": "#8a6fb3"}
BUSCO_COL = {"S": "#2f4b7c", "D": "#7aa6c2", "F": "#f0a35e", "M": "#d1495b"}
BUSCO_LAB = {"S": "Complete, single-copy", "D": "Complete, duplicated", "F": "Fragmented", "M": "Missing"}


def label(sp):
    return sp.replace("_", ". ", 1)


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(p, rows):
    if not rows:
        return
    cols = list(rows[0].keys())
    for r in rows[1:]:
        cols += [c for c in r if c not in cols]
    with open(p, "w", newline="") as fo:
        w = csv.DictWriter(fo, cols, delimiter="\t", extrasaction="ignore", restval="")
        w.writeheader(); w.writerows(rows)


def busco(path_glob):
    f = sorted(glob.glob(path_glob))
    if not f:
        return None
    txt = open(f[0]).read()
    m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)", txt)
    if not m:
        return None
    C, S, D, F, M, n = m.groups()
    return {"C": float(C), "S": float(S), "D": float(D), "F": float(F), "M": float(M), "n": int(n)}


def save(fig, figdir, name):
    os.makedirs(figdir, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(figdir, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def busco_panel(ax, names, data, title):
    y = np.arange(len(names))
    left = np.zeros(len(names))
    for k in ("S", "D", "F", "M"):
        v = np.array([data[s][k] for s in names])
        ax.barh(y, v, left=left, color=BUSCO_COL[k], label=BUSCO_LAB[k], edgecolor="white", lw=0.5)
        for i, (l, x) in enumerate(zip(left, v)):
            if x >= 6:
                ax.text(l + x / 2, i, f"{x:.1f}", ha="center", va="center", fontsize=7,
                        color="white" if k in ("S", "M") else "black")
        left += v
    ax.set_yticks(y, [label(s) for s in names], style="italic")
    ax.invert_yaxis(); ax.set_xlim(0, 100); ax.set_xlabel("% of BUSCOs")
    n = data[names[0]]["n"]
    ax.set_title(f"{title} (eurotiales_odb12.2, n = {n:,})")
    ax.legend(ncol=4, fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22))


def cmd_genome(a):
    rows, bus = [], {}
    for sp in a.species:
        p = os.path.join(a.qc_dir, sp, f"{sp}.genome_stats.tsv")
        if not os.path.exists(p):
            continue
        r = read_tsv(p)[0]
        b = busco(os.path.join(a.qc_dir, sp, "busco", "short_summary*.txt"))
        row = {"species": sp, "strain": r["strain"], "contigs": int(r["contigs"]),
               "size_Mb": round(int(r["total_bp"]) / 1e6, 2), "N50_kb": round(int(r["N50"]) / 1e3, 1),
               "L50": int(r["L50"]), "N90_kb": round(int(r["N90"]) / 1e3, 1),
               "largest_kb": round(int(r["largest"]) / 1e3, 1), "GC_percent": r["GC_percent"],
               "N_percent": r["N_percent"], "contigs_lt200bp": r["contigs_lt200bp"]}
        if b:
            row.update({f"BUSCO_{k}": b[k] for k in ("C", "S", "D", "F", "M")})
            row["BUSCO_n"] = b["n"]; bus[sp] = b
        rows.append(row)
    if not rows:
        print("no species finished yet"); return
    # A. latus vs parents
    size = {r["species"]: r["size_Mb"] for r in rows}
    if all(k in size for k in ("A_latus", "A_spinulosporus", "A_quadrilineatus")):
        ratio = size["A_latus"] / ((size["A_spinulosporus"] + size["A_quadrilineatus"]) / 2)
        for r in rows:
            if r["species"] == "A_latus":
                r["size_vs_mean_parents"] = round(ratio, 2)
        print(f"A. latus genome = {ratio:.2f} x mean parental size (allodiploid expected ~2)")
    write_tsv(os.path.join(a.qc_dir, "genome_qc_summary.tsv"), rows)

    names = [r["species"] for r in rows]
    fig = plt.figure(figsize=(11, 6.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.65, wspace=0.45)
    cols = [SPCOL.get(s, "grey") for s in names]
    x = np.arange(len(names))
    for i, (key, ylab, title) in enumerate([("size_Mb", "Assembly size (Mb)", "Genome size"),
                                            ("N50_kb", "N50 (kb)", "Contiguity"),
                                            ("GC_percent", "GC (%)", "GC content")]):
        ax = fig.add_subplot(gs[0, i])
        v = [float(r[key]) for r in rows]
        ax.bar(x, v, color=cols, width=0.65)
        for xi, vi, r in zip(x, v, rows):
            txt = f"{vi:,.1f}" if key != "N50_kb" else f"{vi:,.0f}\n({r['contigs']} contigs)"
            ax.text(xi, vi, txt, ha="center", va="bottom", fontsize=7)
        ax.set_xticks(x, [label(s) for s in names], rotation=30, ha="right", style="italic")
        ax.set_ylabel(ylab); ax.set_title(title)
        if key == "GC_percent":
            ax.set_ylim(min(v) - 2, max(v) + 1.5)
        else:
            ax.set_ylim(0, max(v) * 1.25)
        ax.text(-0.25, 1.08, "ABC"[i], transform=ax.transAxes, fontsize=13, fontweight="bold")
    if bus:
        ax = fig.add_subplot(gs[1, :])
        busco_panel(ax, [s for s in names if s in bus], bus, "Genome completeness, BUSCO")
        ax.text(-0.08, 1.08, "D", transform=ax.transAxes, fontsize=13, fontweight="bold")
    save(fig, a.fig_dir, "Fig_genome_qc")


def fasta_lengths(p):
    L, n = {}, None
    for l in open(p):
        if l.startswith(">"):
            n = l[1:].split()[0]; L[n] = 0
        else:
            L[n] += len(l.strip())
    return L


def gtf_models(p):
    ex = defaultdict(list); tx2g = {}
    for l in open(p):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        tx2g[t] = g; ex[t].append((int(f[3]), int(f[4])))
    return ex, tx2g


def cmd_annotation(a):
    rows, bus, plen, nex, glen, intr = [], {}, {}, {}, {}, {}
    for sp in a.species:
        d = os.path.join(a.ref_dir, sp)
        s = os.path.join(d, f"{sp}.annotation_summary.tsv")
        if not os.path.exists(s):
            continue
        r = read_tsv(s)[0]
        rows.append(r)
        b = busco(os.path.join(d, "busco_proteins", "short_summary*.txt"))
        if b:
            bus[sp] = b
        plen[sp] = list(fasta_lengths(os.path.join(d, f"{sp}.proteins.longest.faa")).values())
        ex, tx2g = gtf_models(os.path.join(d, f"{sp}.gtf"))
        t1 = [t for t in ex if t.endswith("-T1")]
        nex[sp] = [len(ex[t]) for t in t1]
        glen[sp] = [max(e for _, e in ex[t]) - min(s for s, _ in ex[t]) + 1 for t in t1]
        intr[sp] = [b2[0] - b1[1] - 1 for t in t1 for b1, b2 in zip(sorted(ex[t]), sorted(ex[t])[1:])]
    if not rows:
        print("no species finished yet"); return
    write_tsv(os.path.join(a.ref_dir, "annotation_summary_all.tsv"), rows)
    names = [r["species"] for r in rows]
    cols = [SPCOL.get(s, "grey") for s in names]
    x = np.arange(len(names))
    fig = plt.figure(figsize=(11, 7))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.7, wspace=0.45)

    ax = fig.add_subplot(gs[0, 0])
    bottoms = np.zeros(len(names))
    subcols = {"S": "#c8553d", "Q": "#4b9b6e", "ambiguous": "#bdbdbd", "unassigned": "#e6e6e6"}
    for i, r in enumerate(rows):
        subs = {k[len("subgenome_"):]: int(v) for k, v in r.items() if k.startswith("subgenome_") and v}
        if subs:
            b = 0
            for k in ("S", "Q", "ambiguous", "unassigned"):
                if k in subs:
                    ax.bar(i, subs[k], bottom=b, color=subcols[k], width=0.65, label=f"A. latus {k}")
                    b += subs[k]
        else:
            ax.bar(i, int(r["genes"]), color=cols[i], width=0.65)
        ax.text(i, int(r["genes"]), f"{int(r['genes']):,}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x, [label(s) for s in names], rotation=30, ha="right", style="italic")
    ax.set_ylabel("Protein-coding genes"); ax.set_title("Gene models")
    ax.set_ylim(0, max(int(r["genes"]) for r in rows) * 1.2)
    h, l = ax.get_legend_handles_labels()
    if h:
        ax.legend(dict(zip(l, h)).values(), dict(zip(l, h)).keys(), fontsize=6.5, frameon=False, loc="upper left")

    for j, (data, ylab, title, logy) in enumerate([(plen, "Protein length (aa)", "Protein length", True),
                                                    (nex, "Exons per gene", "Exons per gene", False)]):
        ax = fig.add_subplot(gs[0, j + 1])
        parts = ax.violinplot([data[s] for s in names], positions=x, showmedians=True, widths=0.75)
        for pc, c in zip(parts["bodies"], cols):
            pc.set_facecolor(c); pc.set_alpha(0.6); pc.set_edgecolor("none")
        for k in ("cmedians", "cmins", "cmaxes", "cbars"):
            parts[k].set_color("#333333"); parts[k].set_linewidth(0.7)
        if logy:
            ax.set_yscale("log")
        else:
            ax.set_ylim(0, np.percentile(np.concatenate([data[s] for s in names]), 99.5) + 1)
        for xi, s in zip(x, names):
            ax.text(xi, ax.get_ylim()[1], f"median {st.median(data[s]):.0f}", ha="center", va="bottom", fontsize=6.5)
        ax.set_xticks(x, [label(s) for s in names], rotation=30, ha="right", style="italic")
        ax.set_ylabel(ylab); ax.set_title(title, pad=14)

    ax = fig.add_subplot(gs[1, 0])
    bins = np.logspace(np.log10(10), np.log10(5000), 60)
    for s, c in zip(names, cols):
        if intr[s]:
            ax.hist([v for v in intr[s] if v >= 10], bins=bins, histtype="step", color=c, lw=1.2,
                    density=True, label=label(s))
    ax.set_xscale("log"); ax.set_xlabel("Intron length (bp)"); ax.set_ylabel("Density")
    ax.set_title("Intron length"); ax.legend(fontsize=7, frameon=False, prop={"style": "italic", "size": 7})

    if bus:
        ax = fig.add_subplot(gs[1, 1:])
        busco_panel(ax, [s for s in names if s in bus], bus, "Proteome completeness, BUSCO")
    for i, axx in enumerate(fig.axes):
        axx.text(-0.22 if i < 4 else -0.1, 1.1, "ABCDE"[i], transform=axx.transAxes, fontsize=13, fontweight="bold")
    save(fig, a.fig_dir, "Fig_annotation_qc")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    g = sp.add_parser("genome"); g.add_argument("--qc-dir", required=True)
    g.add_argument("--fig-dir", required=True); g.add_argument("--species", nargs="+", required=True)
    n = sp.add_parser("annotation"); n.add_argument("--ref-dir", required=True)
    n.add_argument("--fig-dir", required=True); n.add_argument("--species", nargs="+", required=True)
    a = ap.parse_args()
    {"genome": cmd_genome, "annotation": cmd_annotation}[a.cmd](a)


if __name__ == "__main__":
    main()
