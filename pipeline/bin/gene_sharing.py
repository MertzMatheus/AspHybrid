#!/usr/bin/env python3
"""
gene_sharing.py - shared and unique genes between species (from OrthoFinder).

Species are defined as groups of OrthoFinder "species", so the hybrid is
counted once at the species level:
    --groups A_fumigatus=A_fumigatus A_spinulosporus=A_spinulosporus \
             A_latus=A_latus_S,A_latus_Q
An extra reference (e.g. A_quadrilineatus) can be given with --reference:
it is NOT part of the sharing statistics, but genes that look species-
specific are flagged if this reference has them (for A. latus, Q-subgenome
genes absent from A. fumigatus/A. spinulosporus are often shared with the
quadrilineatus-like parent, not truly novel).

Sharing unit: orthogroup (gene family). A gene is "shared with species X"
when its orthogroup contains at least one gene of X.

Outputs (--outdir):
  sharing_by_species.tsv      per species: genes in each sharing category, % of genome
  sharing_patterns.tsv        every presence/absence pattern: orthogroups and genes per species
  sharing_pairwise.tsv        % of genes of A whose family is present in B
  sharing_subgenomes.tsv      same patterns at subgenome level (Af, As, Al_S, Al_Q)
  gene_membership.tsv         EVERY gene: species, locus_tag, orthogroup, category, pattern,
                              copies per species, flag for the reference species
  unique_genes_<species>.tsv  species-specific genes (lists for the paper / enrichment)
  sharing_venn.png/.pdf       Venn (orthogroups) + % genes per category
"""
import argparse
import csv
import re
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

csv.field_size_limit(10**9)


def locus(x):
    return re.sub(r"-T\d+$", "", x.strip())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--of", required=True)
    ap.add_argument("--groups", nargs="+", required=True)
    ap.add_argument("--reference", default=None, help="OrthoFinder species used only as a flag")
    ap.add_argument("--short", nargs="*", default=[], help="labels e.g. A_fumigatus=Af")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    groups = {}
    for g in a.groups:
        k, v = g.split("=", 1)
        groups[k] = v.split(",")
    sp_order = list(groups)
    member_of = {m: k for k, v in groups.items() for m in v}
    short = dict(s.split("=", 1) for s in a.short) if a.short else {s: s for s in sp_order}

    ogs = []   # (og, {group: [genes]}, {raw species: [genes]}, ref genes)
    of = Path(a.of)
    for fn in ("Orthogroups/Orthogroups.tsv", "Orthogroups/Orthogroups_UnassignedGenes.tsv"):
        p = of / fn
        if not p.exists():
            continue
        with open(p) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                byg = defaultdict(list)
                raw = {}
                for col in rd.fieldnames[1:]:
                    genes = [locus(x) for x in (r.get(col) or "").split(",") if x.strip()]
                    raw[col] = genes
                    if col in member_of:
                        byg[member_of[col]] += genes
                ref = raw.get(a.reference, []) if a.reference else []
                if any(byg.values()):
                    ogs.append((r["Orthogroup"], byg, raw, ref))

    def pattern(present):
        return "+".join(short[s] for s in sp_order if s in present)

    # ---- per gene membership ----
    n_all = len(sp_order)
    cat_genes = defaultdict(Counter)       # species -> category -> genes
    pat_rows = defaultdict(lambda: {"orthogroups": 0, **{s: 0 for s in sp_order}})
    pair_hits = Counter()
    total = Counter()
    uniq = defaultdict(list)
    with open(out / "gene_membership.tsv", "w") as fo:
        cols = ["species", "locus_tag", "orthogroup", "category", "pattern"] + \
               [f"copies_{short[s]}" for s in sp_order]
        if a.reference:
            cols.append(f"in_{a.reference}")
        fo.write("\t".join(cols) + "\n")
        for og, byg, raw, ref in ogs:
            present = {s for s in sp_order if byg.get(s)}
            pat = pattern(present)
            pat_rows[pat]["orthogroups"] += 1
            for s in sp_order:
                for g in byg.get(s, []):
                    pat_rows[pat][s] += 1
                    total[s] += 1
                    if len(present) == n_all:
                        cat = "core (all species)"
                    elif len(present) == 1:
                        cat = "species-specific" + (" (multi-copy family)" if len(byg[s]) > 1 else " (singleton)")
                    else:
                        cat = "shared with " + " & ".join(short[x] for x in sp_order if x in present and x != s) + " only"
                    cat_genes[s][cat] += 1
                    for t in present - {s}:
                        pair_hits[(s, t)] += 1
                    row = [s, g, og, cat, pat] + [str(len(byg.get(x, []))) for x in sp_order]
                    if a.reference:
                        row.append("yes" if ref else "no")
                    fo.write("\t".join(row) + "\n")
                    if len(present) == 1:
                        uniq[s].append((g, og, "yes" if ref else "no"))

    with open(out / "sharing_by_species.tsv", "w") as fo:
        fo.write("species\tcategory\tgenes\tpercent_of_genes\n")
        for s in sp_order:
            for cat, n in sorted(cat_genes[s].items(), key=lambda x: (not x[0].startswith("core"), x[0])):
                fo.write(f"{s}\t{cat}\t{n}\t{100 * n / total[s] if total[s] else 0:.2f}\n")
            fo.write(f"{s}\tTOTAL\t{total[s]}\t100.00\n")

    with open(out / "sharing_patterns.tsv", "w") as fo:
        fo.write("pattern\tn_species\torthogroups\t" + "\t".join(f"genes_{short[s]}" for s in sp_order) + "\n")
        for pat, r in sorted(pat_rows.items(), key=lambda x: (-x[0].count("+"), x[0])):
            fo.write(f"{pat}\t{pat.count('+') + 1}\t{r['orthogroups']}\t" + "\t".join(str(r[s]) for s in sp_order) + "\n")

    with open(out / "sharing_pairwise.tsv", "w") as fo:
        fo.write("species_A\tspecies_B\tgenes_A\tgenes_A_with_family_in_B\tpercent\n")
        for s in sp_order:
            for t in sp_order:
                if s != t:
                    fo.write(f"{s}\t{t}\t{total[s]}\t{pair_hits[(s, t)]}\t{100 * pair_hits[(s, t)] / total[s] if total[s] else 0:.2f}\n")

    for s, lst in uniq.items():
        with open(out / f"unique_genes_{s}.tsv", "w") as fo:
            fo.write("locus_tag\torthogroup" + (f"\tin_{a.reference}" if a.reference else "") + "\n")
            for g, og, r in sorted(lst):
                fo.write(f"{g}\t{og}" + (f"\t{r}" if a.reference else "") + "\n")

    # ---- subgenome level ----
    raw_species = sorted({m for v in groups.values() for m in v}, key=lambda m: (sp_order.index(member_of[m]), m))
    sub_pat = Counter()
    sub_genes = defaultdict(Counter)
    for og, byg, raw, ref in ogs:
        present = tuple(m for m in raw_species if raw.get(m))
        sub_pat[present] += 1
        for m in present:
            sub_genes[present][m] += len(raw[m])
    with open(out / "sharing_subgenomes.tsv", "w") as fo:
        fo.write("pattern\torthogroups\t" + "\t".join(f"genes_{m}" for m in raw_species) + "\n")
        for p, n in sorted(sub_pat.items(), key=lambda x: (-len(x[0]), x[0])):
            fo.write("+".join(p) + f"\t{n}\t" + "\t".join(str(sub_genes[p][m]) for m in raw_species) + "\n")

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={"width_ratios": [1, 1.25]})
        cols = ["#3b6ea5", "#c8553d", "#4b9b6e", "#8a6fb3"]
        if n_all == 3:
            centers = [(-0.6, 0.35), (0.6, 0.35), (0, -0.6)]
            for (x, y), c, s in zip(centers, cols, sp_order):
                ax1.add_patch(Circle((x, y), 1.05, alpha=0.25, color=c, lw=0))
                ax1.add_patch(Circle((x, y), 1.05, fill=False, color=c, lw=1.5))
            lab = [short[s] for s in sp_order]
            pos = {(0,): (-1.15, 0.7), (1,): (1.15, 0.7), (2,): (0, -1.25), (0, 1): (0, 0.85),
                   (0, 2): (-0.75, -0.35), (1, 2): (0.75, -0.35), (0, 1, 2): (0, 0.05)}
            for idx, (x, y) in pos.items():
                pat = "+".join(lab[i] for i in idx)
                ax1.text(x, y, f"{pat_rows[pat]['orthogroups'] if pat in pat_rows else 0:,}", ha="center", va="center", fontsize=11)
            lab_pos = [(-1.05, 1.62), (1.05, 1.62), (0, -1.95)]
            for (dx, dy), s, c in zip(lab_pos, sp_order, cols):
                ax1.text(dx, dy, s.replace("_", ". ", 1), ha="center", va="center", fontsize=11, style="italic", color=c)
            ax1.set_xlim(-2.3, 2.3); ax1.set_ylim(-2.2, 2.0); ax1.set_aspect("equal"); ax1.axis("off")
            ax1.set_title("Orthogroups (gene families)", fontsize=12)
        else:
            ax1.axis("off")
        cats = ["core (all species)"]
        others = sorted({c for s in sp_order for c in cat_genes[s] if c.startswith("shared")})
        cats += others + ["species-specific (multi-copy family)", "species-specific (singleton)"]
        pal = ["#2f4b7c", "#7aa6c2", "#a5c8dd", "#c7dbe8", "#f0a35e", "#d1495b"]
        left = [0] * n_all
        for k, cat in enumerate(cats):
            vals = [100 * cat_genes[s][cat] / total[s] if total[s] else 0 for s in sp_order]
            ax2.barh(range(n_all), vals, left=left, color=pal[k % len(pal)], label=cat, edgecolor="white", lw=0.5)
            for i, v in enumerate(vals):
                if v >= 4:
                    ax2.text(left[i] + v / 2, i, f"{v:.0f}%", ha="center", va="center", fontsize=9,
                             color="white" if k == 0 else "black")
            left = [l + v for l, v in zip(left, vals)]
        ax2.set_yticks(range(n_all))
        ax2.set_yticklabels([f"{s.replace('_', '. ', 1)}\n({total[s]:,} genes)" for s in sp_order], style="italic")
        ax2.invert_yaxis(); ax2.set_xlim(0, 100); ax2.set_xlabel("% of protein-coding genes")
        ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False, fontsize=8.5)
        ax2.spines[["top", "right"]].set_visible(False)
        ax2.set_title("Genes by sharing category", fontsize=12)
        plt.tight_layout()
        for ext in ("png", "pdf"):
            plt.savefig(out / f"sharing_venn.{ext}", dpi=300, bbox_inches="tight")
    except ImportError:
        print("matplotlib not available: figure skipped")

    for s in sp_order:
        if not total[s]:
            print(f"{s:18s} no genes found in OrthoFinder results (check --groups)"); continue
        core = cat_genes[s]["core (all species)"]
        print(f"{s:18s} genes={total[s]:>6}  core={100 * core / total[s]:5.1f}%  "
              f"specific={100 * sum(v for k, v in cat_genes[s].items() if k.startswith('species')) / total[s]:5.1f}%")


if __name__ == "__main__":
    main()
