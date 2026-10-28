#!/usr/bin/env python3
"""
functional_db.py - one functional annotation row per gene + enrichment gene sets.

Inputs: proteins.faa of the species (product= in the headers; official RefSeq
products for A. fumigatus), eggNOG-mapper .emapper.annotations (one protein per
gene, -T1), go-basic.obo, KEGG pathway names, id_map (subgenome), and the
orthology table of step 07 (optional).

GO: eggNOG GO terms are mapped to primary IDs (alt_id), obsolete terms dropped,
and propagated to all ancestors (is_a + part_of) so that over-representation
tests are valid for parent terms. Root terms are not used as gene sets.

Outputs:
  PREFIX.gene_annotation.tsv       one row per gene (all genes, also without hit)
  PREFIX.annotation_coverage.tsv   number / % of genes with each annotation type
  GENESET_DIR/<sp>.<coll>.TERM2GENE.tsv (term, gene) and .TERM2NAME.tsv (term, name)
     coll = GO_BP GO_MF GO_CC KEGG_pathway COG KO Pfam CAZy
"""
import argparse, csv, os, re
from collections import defaultdict

csv.field_size_limit(10**9)
COG = {
    "A": "RNA processing and modification", "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion", "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism", "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism", "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism", "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription", "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis", "N": "Cell motility",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "S": "Function unknown", "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms", "W": "Extracellular structures", "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}
NS = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}
ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}
locus = lambda t: re.sub(r"-T\d+$", "", t)
ok = lambda v: v not in (None, "", "-", "NA")


def parse_obo(path):
    terms, cur, alt = {}, None, {}
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("["):
            if cur and "id" in cur and "ns" in cur: terms[cur["id"]] = cur
            cur = {"is_a": [], "alt": [], "obsolete": False} if line == "[Term]" else None
            continue
        if cur is None or ": " not in line: continue
        k, v = line.split(": ", 1)
        if k == "id": cur["id"] = v
        elif k == "name": cur["name"] = v
        elif k == "namespace": cur["ns"] = NS.get(v, v)
        elif k == "is_a": cur["is_a"].append(v.split()[0])
        elif k == "relationship" and v.startswith("part_of "): cur["is_a"].append(v.split()[1])
        elif k == "alt_id": cur["alt"].append(v)
        elif k == "is_obsolete" and v == "true": cur["obsolete"] = True
    if cur and "id" in cur and "ns" in cur: terms[cur["id"]] = cur
    for t, d in terms.items():
        for x in d["alt"]: alt[x] = t
    return terms, alt


ap = argparse.ArgumentParser()
ap.add_argument("--species", required=True); ap.add_argument("--proteins", required=True)
ap.add_argument("--eggnog", required=True); ap.add_argument("--obo", required=True)
ap.add_argument("--kegg-names"); ap.add_argument("--id-map"); ap.add_argument("--correspondence")
ap.add_argument("--out", required=True); ap.add_argument("--geneset-dir", required=True)
a = ap.parse_args()
os.makedirs(a.geneset_dir, exist_ok=True)

# genes and products (T1 header)
genes, product = [], {}
for l in open(a.proteins):
    if l.startswith(">"):
        tid = l[1:].split()[0]; g = locus(tid)
        if g not in product: genes.append(g)
        m = re.search(r"product=(.*?)(?: \w+=|$)", l.rstrip("\n"))
        if tid.endswith("-T1") or g not in product:
            product[g] = m.group(1).strip() if m else ""
sub = {}
if a.id_map and os.path.exists(a.id_map):
    for r in csv.DictReader(open(a.id_map), delimiter="\t"):
        sub[r["locus_tag"]] = r.get("subgenome", "NA") or "NA"

# eggNOG
egg, hdr = {}, None
for line in open(a.eggnog):
    if line.startswith("#query"):
        hdr = line[1:].rstrip("\n").split("\t"); continue
    if line.startswith("#") or hdr is None: continue
    r = dict(zip(hdr, line.rstrip("\n").split("\t")))
    g = locus(r["query"])
    if g not in egg or r["query"].endswith("-T1"): egg[g] = r

terms, alt = parse_obo(a.obo)
memo = {}
def anc(t):
    if t not in memo:
        out = {t}
        for p in terms.get(t, {}).get("is_a", []):
            if p in terms: out |= anc(p)
        memo[t] = out
    return memo[t]

kegg_names = {}
if a.kegg_names and os.path.exists(a.kegg_names):
    for l in open(a.kegg_names):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 2:
            k = re.sub(r"^path:", "", f[0]); k = re.sub(r"^(ko|map)", "map", k)
            kegg_names[k] = re.sub(r" - .*$", "", f[1])

sets = defaultdict(lambda: defaultdict(set))
go_direct = defaultdict(set)
for g, r in egg.items():
    for go in r.get("GOs", "-").split(","):
        go = alt.get(go, go)
        if go in terms and not terms[go]["obsolete"]:
            go_direct[g].add(go)
            for x in anc(go):
                if x not in ROOTS: sets["GO_" + terms[x]["ns"]][g].add(x)
    for pw in r.get("KEGG_Pathway", "-").split(","):
        if pw.startswith("map"): sets["KEGG_pathway"][g].add(pw)
    for ko in r.get("KEGG_ko", "-").split(","):
        if ko.startswith("ko:"): sets["KO"][g].add(ko[3:])
    for c in r.get("COG_category", "-"):
        if c in COG: sets["COG"][g].add(c)
    for p in r.get("PFAMs", "-").split(","):
        if ok(p): sets["Pfam"][g].add(p)
    for c in r.get("CAZy", "-").split(","):
        if ok(c): sets["CAZy"][g].add(c)

names = {"COG": COG, "KEGG_pathway": kegg_names}
for ns in ("BP", "MF", "CC"):
    names[f"GO_{ns}"] = {t: d["name"] for t, d in terms.items()}
for coll in ("GO_BP", "GO_MF", "GO_CC", "KEGG_pathway", "COG", "KO", "Pfam", "CAZy"):
    t2g = defaultdict(set)
    for g, ts in sets[coll].items():
        for t in ts: t2g[t].add(g)
    with open(os.path.join(a.geneset_dir, f"{a.species}.{coll}.TERM2GENE.tsv"), "w") as fo:
        fo.write("term\tgene\n")
        for t in sorted(t2g):
            for g in sorted(t2g[t]): fo.write(f"{t}\t{g}\n")
    with open(os.path.join(a.geneset_dir, f"{a.species}.{coll}.TERM2NAME.tsv"), "w") as fo:
        fo.write("term\tname\n")
        for t in sorted(t2g): fo.write(f"{t}\t{names.get(coll, {}).get(t, t)}\n")

# orthologues from step 07
orth_cols, orth = [], {}
if a.correspondence and os.path.exists(a.correspondence):
    rd = csv.DictReader(open(a.correspondence), delimiter="\t")
    orth_cols = [c for c in rd.fieldnames if c.endswith("_orthologs") or c.endswith("_relation")]
    for r in rd: orth[r["locus_tag"]] = r

# product: official RefSeq product (A. fumigatus); otherwise the eggNOG description
src = {}
for g in genes:
    if product.get(g):
        src[g] = "RefSeq"
    elif ok(egg.get(g, {}).get("Description")):
        product[g] = egg[g]["Description"]; src[g] = "eggNOG"
    else:
        product[g] = "hypothetical protein"; src[g] = "none"

cols = ["locus_tag", "subgenome", "product", "product_source", "eggNOG_preferred_name", "eggNOG_description", "eggNOG_OG",
        "max_annot_lvl", "COG_category", "COG_description", "EC", "KEGG_ko", "KEGG_pathway", "KEGG_pathway_names",
        "GO_terms", "Pfam", "CAZy", "eggNOG_seed_ortholog", "eggNOG_evalue"] + orth_cols
with open(f"{a.out}.gene_annotation.tsv", "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t"); w.writerow(cols)
    for g in genes:
        r = egg.get(g, {})
        v = lambda k: r.get(k, "") if ok(r.get(k)) else ""
        og = v("eggNOG_OGs").split(",")[-1].split("|")[0] if v("eggNOG_OGs") else ""
        cogs = "".join(sorted(sets["COG"].get(g, ())))
        pws = sorted(sets["KEGG_pathway"].get(g, ()))
        row = [g, sub.get(g, "NA"), product.get(g, ""), src[g], v("Preferred_name"), v("Description"), og,
               v("max_annot_lvl"), cogs, "; ".join(COG[c] for c in cogs), v("EC"),
               ",".join(sorted(sets["KO"].get(g, ()))), ",".join(pws),
               "; ".join(kegg_names.get(p, p) for p in pws), ",".join(sorted(go_direct.get(g, ()))),
               ",".join(sorted(sets["Pfam"].get(g, ()))), ",".join(sorted(sets["CAZy"].get(g, ()))),
               v("seed_ortholog"), v("evalue")] + [orth.get(g, {}).get(c, "") for c in orth_cols]
        w.writerow(row)

n = len(genes)
hyp = lambda p: (not p) or re.search(r"hypothetical|uncharacteri[sz]ed", p, re.I)
cov = [("genes", n), ("eggNOG hit", sum(1 for g in genes if g in egg)),
       ("eggNOG preferred name", sum(1 for g in genes if ok(egg.get(g, {}).get("Preferred_name")))),
       ("eggNOG description", sum(1 for g in genes if ok(egg.get(g, {}).get("Description")))),
       ("named product (not hypothetical)", sum(1 for g in genes if not hyp(product.get(g, "")))),
       ("COG category", sum(1 for g in genes if sets["COG"].get(g))),
       ("COG S (function unknown) only", sum(1 for g in genes if sets["COG"].get(g) == {"S"})),
       ("GO (any)", sum(1 for g in genes if go_direct.get(g))),
       ("GO biological process", sum(1 for g in genes if sets["GO_BP"].get(g))),
       ("GO molecular function", sum(1 for g in genes if sets["GO_MF"].get(g))),
       ("GO cellular component", sum(1 for g in genes if sets["GO_CC"].get(g))),
       ("EC number", sum(1 for g in genes if ok(egg.get(g, {}).get("EC")))),
       ("KEGG KO", sum(1 for g in genes if sets["KO"].get(g))),
       ("KEGG pathway", sum(1 for g in genes if sets["KEGG_pathway"].get(g))),
       ("Pfam domain", sum(1 for g in genes if sets["Pfam"].get(g))),
       ("CAZy", sum(1 for g in genes if sets["CAZy"].get(g)))]
with open(f"{a.out}.annotation_coverage.tsv", "w") as fo:
    fo.write("species\tannotation\tgenes\tpercent_of_genes\n")
    for k, v in cov:
        fo.write(f"{a.species}\t{k}\t{v}\t{100 * v / n:.2f}\n")
with open(f"{a.out}.COG_counts.tsv", "w") as fo:
    fo.write("species\tCOG\tdescription\tgenes\tpercent_of_genes\n")
    for c in COG:
        k = sum(1 for g in genes if c in sets["COG"].get(g, ()))
        fo.write(f"{a.species}\t{c}\t{COG[c]}\t{k}\t{100 * k / n:.2f}\n")
for k, v in cov:
    print(f"  {k:36s} {v:>7}  {100 * v / n:5.1f}%")
