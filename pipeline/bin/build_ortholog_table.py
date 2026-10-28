#!/usr/bin/env python3
"""
build_ortholog_table.py - gene correspondence tables from OrthoFinder.

Answers "gene X in A. latus = which gene in A. spinulosporus / A. fumigatus?"
using OrthoFinder's gene-tree orthologues (more precise than orthogroups).
All IDs are locus_tags (transcript suffix -T1 removed).

Outputs (in --outdir):
  correspondence_<species>.tsv     one row per gene of that species, with its
                                   orthologues in every other species and the
                                   relation (1:1, 1:many, many:1, many:many)
  correspondence_A_latus.tsv       S and Q subgenome genes together, with the
                                   homeolog in the other subgenome
  correspondence_by_orthogroup.tsv one row per orthogroup, genes per species
  orthology_summary.tsv            counts for the paper (pairs and categories)

Usage:
  build_ortholog_table.py --of Results_XXX --outdir DIR \
      --products A_fumigatus=path.tsv A_latus=path.tsv ...
"""
import argparse
import csv
import re
from collections import defaultdict, Counter
from pathlib import Path

csv.field_size_limit(10**9)


def locus(x):
    return re.sub(r"-T\d+$", "", x.strip())


def split(cell):
    return [locus(g) for g in cell.split(",") if g.strip()]


def read_products(specs):
    prod = {}
    for s in specs or []:
        sp, path = s.split("=", 1)
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                t = r["transcript_id"]
                if t.endswith("-T1") or locus(t) not in prod:
                    prod[locus(t)] = r["product"]
    return prod


def relation(n_a, n_b):
    return f"{'1' if n_a == 1 else 'many'}:{'1' if n_b == 1 else 'many'}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--of", required=True, help="OrthoFinder Results_* directory")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--products", nargs="*")
    ap.add_argument("--subgenome-species", default="A_latus_S,A_latus_Q")
    ap.add_argument("--merged-name", default="A_latus")
    ap.add_argument("--id-map", help="A_latus.curated.id_map.tsv: true subgenome label per locus")
    a = ap.parse_args()
    of, out = Path(a.of), Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    prod = read_products(a.products)
    subS, subQ = a.subgenome_species.split(",")
    true_sub = {}
    if a.id_map:
        with open(a.id_map) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                true_sub[r["locus_tag"]] = r.get("subgenome", "")

    # ---- orthogroups (all genes, including unassigned) ----
    og_of = {}
    og_rows = []
    species = None
    for fn in ("Orthogroups/Orthogroups.tsv", "Orthogroups/Orthogroups_UnassignedGenes.tsv"):
        p = of / fn
        if not p.exists():
            continue
        with open(p) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            sp_cols = [c for c in rd.fieldnames if c != "Orthogroup"]
            species = species or sp_cols
            for r in rd:
                genes = {sp: split(r.get(sp, "") or "") for sp in species}
                og_rows.append((r["Orthogroup"], genes))
                for sp, gs in genes.items():
                    for g in gs:
                        og_of[g] = r["Orthogroup"]
    if not species:
        raise SystemExit("Orthogroups.tsv not found")
    gene_species = {g: sp for _, genes in og_rows for sp, gs in genes.items() for g in gs}

    # ---- pairwise orthologues ----
    pair = defaultdict(lambda: defaultdict(list))   # (A,B) -> geneA -> [geneB]
    rel_counts = Counter()
    for A in species:
        for B in species:
            if A == B:
                continue
            p = of / "Orthologues" / f"Orthologues_{A}" / f"{A}__v__{B}.tsv"
            if not p.exists():
                continue
            with open(p) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    ga, gb = split(r[A]), split(r[B])
                    if A < B:
                        rel_counts[(A, B, relation(len(ga), len(gb)))] += 1
                    for g in ga:
                        pair[(A, B)][g].extend(gb)

    # ---- one table per anchor species ----
    def anchor_rows(A):
        rows = []
        genes = sorted(g for g, s in gene_species.items() if s == A)
        for g in genes:
            r = {"locus_tag": g, "orthogroup": og_of.get(g, ""), "product": prod.get(g, "")}
            for B in species:
                if B == A:
                    continue
                hits = pair[(A, B)].get(g, [])
                if hits:
                    # how many A genes share these B genes?
                    back = {x for h in hits for x in pair[(B, A)].get(h, [])} or {g}
                    r[f"{B}_orthologs"] = ",".join(sorted(set(hits)))
                    r[f"{B}_relation"] = relation(len(back), len(set(hits)))
                    r[f"{B}_product"] = prod.get(sorted(set(hits))[0], "") if len(set(hits)) == 1 else ""
                else:
                    r[f"{B}_orthologs"] = ""
                    r[f"{B}_relation"] = "none"
                    r[f"{B}_product"] = ""
            rows.append(r)
        return rows

    def write(path, rows, cols=None):
        if not rows:
            return
        cols = cols or list(rows[0].keys())
        with open(path, "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    tables = {}
    for A in species:
        tables[A] = anchor_rows(A)
        write(out / f"correspondence_{A}.tsv", tables[A])

    # ---- merged hybrid table (S + Q) ----
    if subS in species and subQ in species:
        merged = []
        others = [s for s in species if s not in (subS, subQ)]
        cols = ["locus_tag", "subgenome", "orthogroup", "product",
                "homeolog", "homeolog_relation", "homeolog_product"]
        for B in others:
            cols += [f"{B}_orthologs", f"{B}_relation", f"{B}_product"]
        for X, Y, lab in ((subS, subQ, "S"), (subQ, subS, "Q")):
            for r in tables[X]:
                m = {k: v for k, v in r.items()}
                m["subgenome"] = true_sub.get(r["locus_tag"], lab)
                m["homeolog"] = r.get(f"{Y}_orthologs", "")
                m["homeolog_relation"] = r.get(f"{Y}_relation", "none")
                m["homeolog_product"] = r.get(f"{Y}_product", "")
                merged.append(m)
        merged.sort(key=lambda r: r["locus_tag"])
        write(out / f"correspondence_{a.merged_name}.tsv", merged, cols)

    # ---- by orthogroup ----
    with open(out / "correspondence_by_orthogroup.tsv", "w") as fo:
        fo.write("orthogroup\tcategory\t" + "\t".join(f"n_{s}" for s in species) + "\t"
                 + "\t".join(species) + "\n")
        cat_counts = Counter()
        for og, genes in og_rows:
            n = [len(genes[s]) for s in species]
            present = sum(1 for x in n if x)
            if present == len(species) and all(x == 1 for x in n):
                cat = "core_single_copy"
            elif present == len(species):
                cat = "core_multi_copy"
            elif present == 1:
                cat = "species_specific" if sum(n) > 1 else "singleton"
            else:
                cat = "shared_subset"
            cat_counts[cat] += 1
            fo.write(f"{og}\t{cat}\t" + "\t".join(map(str, n)) + "\t"
                     + "\t".join(",".join(genes[s]) for s in species) + "\n")

    # ---- summary ----
    with open(out / "orthology_summary.tsv", "w") as fo:
        fo.write("section\titem\tvalue\n")
        for s in species:
            gs = [g for g, sp in gene_species.items() if sp == s]
            with_orth = sum(1 for g in gs if any(pair[(s, B)].get(g) for B in species if B != s))
            fo.write(f"genes\t{s}\t{len(gs)}\n")
            fo.write(f"genes_with_orthologue_in_another_species\t{s}\t{with_orth}\n")
            for B in species:
                if B != s:
                    n1 = sum(1 for r in tables[s] if r[f"{B}_relation"] == "1:1")
                    fo.write(f"genes_1to1\t{s} -> {B}\t{n1}\n")
        for k, v in sorted(cat_counts.items()):
            fo.write(f"orthogroups\t{k}\t{v}\n")
        for (A, B, rel), v in sorted(rel_counts.items()):
            fo.write(f"orthologue_pairs\t{A} vs {B} {rel}\t{v}\n")
    print(f"Species: {species}")
    print("Orthogroup categories:", dict(cat_counts))


if __name__ == "__main__":
    main()
