#!/usr/bin/env python3
"""
map_afum_afua.py - correspondence between the de novo BRAKER3 A. fumigatus genes
(AFUM_...) and the official Af293 genes (AFUA_...), by CDS overlap on the same
strand of the same genome. Used (1) to translate any result made with AFUM IDs
and (2) as a benchmark of the annotation pipeline against the curated reference.

class per AFUM gene:
  same_gene     >= 50 % of the CDS of both genes overlap (reciprocal)
  partial       overlaps an AFUA gene, but < 50 % reciprocal
  no_AFUA       no CDS overlap with any official gene (new / spurious model)
AFUA genes without any overlapping AFUM gene are listed as 'not_predicted'.

usage: map_afum_afua.py --afum A_fumigatus_BRAKER.gff3 --afua A_fumigatus.gff3 --out DIR
       [--translate table.tsv ...]   adds an AFUA column next to every AFUM ID column
outputs: AFUM_to_AFUA.tsv, AFUA_not_predicted.tsv, AFUM_vs_AFUA_summary.tsv
"""
import argparse, csv, os, re
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--afum", required=True); ap.add_argument("--afua", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--translate", nargs="*", default=[])
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
csv.field_size_limit(10**9)


def cds_by_gene(p):
    tx2g, cds = {}, defaultdict(list)
    rows = []
    for l in open(p):
        if l.startswith("##FASTA"): break
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) != 9: continue
        at = dict(x.split("=", 1) for x in f[8].strip(";").split(";") if "=" in x)
        if f[2] == "mRNA": tx2g[at["ID"]] = at["Parent"]
        elif f[2] == "CDS": rows.append((f, at["Parent"].split(",")[0]))
    for f, t in rows:
        cds[tx2g.get(t, re.sub(r"-T\d+$", "", t))].append((f[0], f[6], int(f[3]), int(f[4])))
    out = {}
    for g, iv in cds.items():
        iv = sorted(set(iv), key=lambda x: x[2]); merged = []
        for s, st, b, e in iv:
            if merged and b <= merged[-1][3] + 1: merged[-1][3] = max(merged[-1][3], e)
            else: merged.append([s, st, b, e])
        out[g] = merged
    return out


A, B = cds_by_gene(a.afum), cds_by_gene(a.afua)
blen = lambda iv: sum(e - b + 1 for _, _, b, e in iv)
bins = defaultdict(set)
for g, iv in B.items():
    for s, st, b, e in iv:
        for k in range(b // 10000, e // 10000 + 1): bins[(s, st, k)].add(g)


def overlap(x, y):
    return sum(max(0, min(e1, e2) - max(b1, b2) + 1) for s1, t1, b1, e1 in x for s2, t2, b2, e2 in y
               if s1 == s2 and t1 == t2)


rows, hit_afua = [], set()
for g, iv in A.items():
    cand = {h for s, st, b, e in iv for k in range(b // 10000, e // 10000 + 1) for h in bins[(s, st, k)]}
    hits = sorted(((overlap(iv, B[h]), h) for h in cand), reverse=True)
    hits = [(o, h) for o, h in hits if o > 0]
    hit_afua |= {h for _, h in hits}
    if hits:
        o, h = hits[0]
        fa, fb = o / blen(iv), o / blen(B[h])
        cls = "same_gene" if fa >= 0.5 and fb >= 0.5 else "partial"
        rows.append([g, h, cls, o, f"{fa:.3f}", f"{fb:.3f}", ",".join(x for _, x in hits)])
    else:
        rows.append([g, "", "no_AFUA", 0, "0", "0", ""])
rows.sort()
with open(os.path.join(a.out, "AFUM_to_AFUA.tsv"), "w") as fo:
    fo.write("AFUM\tAFUA\tclass\tCDS_overlap_bp\tfraction_of_AFUM_CDS\tfraction_of_AFUA_CDS\tall_overlapping_AFUA\n")
    for r in rows: fo.write("\t".join(map(str, r)) + "\n")
notpred = sorted(set(B) - hit_afua)
with open(os.path.join(a.out, "AFUA_not_predicted.tsv"), "w") as fo:
    fo.write("AFUA\n" + "".join(f"{g}\n" for g in notpred))
same = {r[1] for r in rows if r[2] == "same_gene"}
cnt = lambda c: sum(1 for r in rows if r[2] == c)
summ = [("official_AFUA_genes", len(B)), ("BRAKER_AFUM_genes", len(A)),
        ("AFUM_same_gene", cnt("same_gene")), ("AFUM_partial", cnt("partial")), ("AFUM_no_AFUA", cnt("no_AFUA")),
        ("AFUA_recovered_same_gene", len(same)), ("AFUA_recovered_percent", round(100 * len(same) / len(B), 2)),
        ("AFUA_not_predicted", len(notpred))]
with open(os.path.join(a.out, "AFUM_vs_AFUA_summary.tsv"), "w") as fo:
    fo.write("metric\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in summ))
for k, v in summ: print(f"  {k:28s} {v}")

# translate tables that use AFUM IDs
m = {r[0]: r[1] for r in rows}
for p in a.translate:
    if not os.path.exists(p): continue
    rd = list(csv.reader(open(p), delimiter="\t"))
    if not rd: continue
    cols = [i for i, c in enumerate(rd[0]) if any(re.match(r"AFUM_\d+", x or "") for x in (r[i] for r in rd[1:50] if i < len(r)))]
    outp = re.sub(r"\.tsv$", "", p) + ".with_AFUA.tsv"
    with open(outp, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        hdr = []
        for i, c in enumerate(rd[0]): hdr += [c] + ([f"{c}_AFUA"] if i in cols else [])
        w.writerow(hdr)
        for r in rd[1:]:
            new = []
            for i, x in enumerate(r):
                new.append(x)
                if i in cols:
                    new.append(",".join(m.get(re.sub(r"-T\d+$", "", y), "") for y in x.split(",") if y))
            w.writerow(new)
    print(f"  translated: {outp}")
