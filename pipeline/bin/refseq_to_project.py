#!/usr/bin/env python3
"""
refseq_to_project.py - official A. fumigatus Af293 annotation (NCBI RefSeq,
locus_tags AFUA_1G00100 ...) converted to the project reference format, on
the project's own genome copy (step 01), so the RNA-seq BAMs stay valid.

 1. Genome check: every RefSeq sequence that carries genes must be IDENTICAL
    (same sequence, case-insensitive) to one sequence of the project genome.
    Sequences are matched by MD5, so different FASTA names are fine.
    If not identical the script stops: the official coordinates would not fit.
 2. Protein-coding genes only (gene_biotype=protein_coding). Gene ID = official
    locus_tag (AFUA_...). Isoforms renamed <locus_tag>-T1, -T2 ... (T1 = longest
    protein); RefSeq IDs kept as attributes (transcript XM_, protein XP_).
 3. Proteins/CDS are translated from the genome; models with an internal stop
    are removed and listed (normally none).

usage: refseq_to_project.py --gff GCF_..._genomic.gff[.gz] --refseq-fasta GCF_..._genomic.fna[.gz]
                            --genome project.genome.fasta --out PREFIX
outputs (same names as curate_annotation.py, so step 06 continues unchanged):
  PREFIX.gff3 .proteins.faa .proteins.longest.faa .cds.fna .id_map.tsv .removed.tsv .stats.tsv
  PREFIX.seqid_map.tsv   RefSeq sequence -> project sequence
"""
import argparse, gzip, hashlib, re, sys
from collections import defaultdict, OrderedDict
from urllib.parse import unquote, quote

ap = argparse.ArgumentParser()
ap.add_argument("--gff", required=True); ap.add_argument("--refseq-fasta", required=True)
ap.add_argument("--genome", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
op = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_fasta(p):
    seqs, name, buf = OrderedDict(), None, []
    for l in op(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if name: seqs[name] = "".join(buf).upper()
            name, buf = l[1:].split()[0], []
        else:
            buf.append(l.strip())
    if name: seqs[name] = "".join(buf).upper()
    return seqs


md5 = lambda s: hashlib.md5(s.encode()).hexdigest()
genome = read_fasta(a.genome)
by_md5 = {md5(s): n for n, s in genome.items()}
ref = read_fasta(a.refseq_fasta)
seqmap = {n: by_md5.get(md5(s)) for n, s in ref.items()}

# ---------- parse GFF ----------
def attrs(s):
    d = {}
    for it in s.strip().strip(";").split(";"):
        if "=" in it:
            k, v = it.split("=", 1); d[k.strip()] = unquote(v.strip())
    return d


genes, tx_of_gene, feats_of_tx, tx_attr = OrderedDict(), defaultdict(list), defaultdict(list), {}
parent = {}
for l in op(a.gff):
    if l.startswith("##FASTA"): break
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if len(f) != 9: continue
    at = attrs(f[8]); fid = at.get("ID"); par = at.get("Parent", "").split(",")[0]
    if f[2] == "gene" and at.get("gene_biotype") == "protein_coding":
        genes[fid] = (f, at)
    elif f[2] in ("mRNA", "transcript") and par:
        parent[fid] = par; tx_of_gene[par].append(fid); tx_attr[fid] = (f, at)
    elif f[2] in ("exon", "CDS") and par:
        feats_of_tx[par].append((f, at))

used = {f[0] for f, _ in genes.values()}
missing = [s for s in used if not seqmap.get(s)]
with open(f"{a.out}.seqid_map.tsv", "w") as fo:
    fo.write("refseq_seqid\tproject_seqid\tlength\tidentical\n")
    for n, s in ref.items():
        fo.write(f"{n}\t{seqmap.get(n) or 'NA'}\t{len(s)}\t{'yes' if seqmap.get(n) else 'no'}\n")
if missing:
    sys.exit(f"ERROR: {len(missing)} RefSeq sequences with genes are not identical to any sequence of "
             f"{a.genome}: {', '.join(sorted(missing)[:5])} ...\n"
             f"The project genome is not the Af293 RefSeq assembly; official coordinates cannot be used. "
             f"See {a.out}.seqid_map.tsv")

CODE = {}
bases = "TCAG"; aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, b1 in enumerate(bases):
    for j, b2 in enumerate(bases):
        for k, b3 in enumerate(bases):
            CODE[b1 + b2 + b3] = aas[16 * i + 4 * j + k]
rc = lambda s: s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
translate = lambda s: "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
STOPS = {"TAA", "TAG", "TGA"}

models, removed, stats = defaultdict(list), [], defaultdict(int)
for gid, (gf, ga) in genes.items():
    lt = ga.get("locus_tag") or gid.replace("gene-", "")
    stats["input_genes"] += 1
    txs = tx_of_gene.get(gid) or [gid]          # CDS directly under the gene
    if gid not in tx_of_gene:
        tx_attr[gid] = (gf, ga)
        feats_of_tx[gid] = feats_of_tx.get(gid, [])
    for t in txs:
        fs = feats_of_tx.get(t, [])
        cds = [x for x in fs if x[0][2] == "CDS"]
        if not cds:
            continue
        strand = gf[6]; seq = genome[seqmap[gf[0]]]
        cds.sort(key=lambda x: int(x[0][3]), reverse=(strand == "-"))
        nt = "".join(seq[int(x[0][3]) - 1:int(x[0][4])] if strand == "+" else rc(seq[int(x[0][3]) - 1:int(x[0][4])])
                     for x in cds)
        ph = int(cds[0][0][7]) if cds[0][0][7] in "012" else 0
        nt = nt[ph:]
        p5 = ph > 0 or not nt.startswith("ATG") or "start_range" in tx_attr[t][1] or "start_range" in ga
        p3 = nt[-3:] not in STOPS
        aa = translate(nt[:len(nt) - len(nt) % 3])
        if aa.endswith("*"): aa = aa[:-1]
        if "*" in aa:
            removed.append((gid, t, gf[0], gf[3], gf[4], strand, "internal_stop")); stats["removed_internal_stop"] += 1
            continue
        cat = cds[0][1]
        models[lt].append(dict(t=t, aa=aa, nt=nt, p5=p5, p3=p3, cds=cds,
                               exons=[x for x in fs if x[0][2] == "exon"], gf=gf, ga=ga, gid=gid,
                               product=cat.get("product") or tx_attr[t][1].get("product") or "hypothetical protein",
                               protein_id=cat.get("protein_id", ""), tx_rs=t.replace("rna-", "")))
    if lt not in models:
        stats["removed_genes_no_valid_model"] += 1

enc = lambda v: quote(str(v), safe=" :/|.,()[]'+-_")
render = lambda d: ";".join(f"{k}={enc(v)}" for k, v in d.items())
wrap = lambda s: "\n".join(s[i:i + 60] for i in range(0, len(s), 60))
gff = open(f"{a.out}.gff3", "w"); faa = open(f"{a.out}.proteins.faa", "w")
lng = open(f"{a.out}.proteins.longest.faa", "w"); cdsf = open(f"{a.out}.cds.fna", "w")
mp = open(f"{a.out}.id_map.tsv", "w")
gff.write("##gff-version 3\n")
for n in genome:
    gff.write(f"##sequence-region {n} 1 {len(genome[n])}\n")
mp.write("old_gene_id\told_transcript_id\tlocus_tag\ttranscript_id\tseqid\tstart\tend\tstrand\tprotein_len\t"
         "partial5\tpartial3\tsubgenome\trefseq_protein_id\tproduct\n")
order = sorted(models, key=lambda lt: (list(genome).index(seqmap[models[lt][0]["gf"][0]]), int(models[lt][0]["gf"][3])))
for lt in order:
    ms = sorted(models[lt], key=lambda m: -len(m["aa"]))
    gf, ga = ms[0]["gf"], ms[0]["ga"]; sid = seqmap[gf[0]]; strand = gf[6]
    gat = OrderedDict(ID=lt, locus_tag=lt, gene_biotype="protein_coding")
    if ga.get("Name") and ga["Name"] != lt: gat["Name"] = ga["Name"]
    if ga.get("Dbxref"): gat["Dbxref"] = ga["Dbxref"]
    if any(m["p5"] or m["p3"] for m in ms): gat["partial"] = "true"
    gff.write("\t".join([sid, "RefSeq", "gene", gf[3], gf[4], ".", strand, ".", render(gat)]) + "\n")
    stats["genes"] += 1; stats["transcripts"] += len(ms)
    stats["partial_5prime"] += any(m["p5"] for m in ms); stats["partial_3prime"] += any(m["p3"] for m in ms)
    for i, m in enumerate(ms, 1):
        tid = f"{lt}-T{i}"
        s = min(int(x[0][3]) for x in m["cds"] + m["exons"]); e = max(int(x[0][4]) for x in m["cds"] + m["exons"])
        ta = OrderedDict(ID=tid, Parent=lt, locus_tag=lt, product=m["product"])
        if m["protein_id"]: ta["refseq_protein_id"] = m["protein_id"]
        ta["refseq_transcript_id"] = m["tx_rs"]
        if m["p5"] or m["p3"]: ta["partial"] = "true"
        gff.write("\t".join([sid, "RefSeq", "mRNA", str(s), str(e), ".", strand, ".", render(ta)]) + "\n")
        ex = m["exons"] or m["cds"]
        for k, x in enumerate(sorted(ex, key=lambda x: int(x[0][3])), 1):
            gff.write("\t".join([sid, "RefSeq", "exon", x[0][3], x[0][4], ".", strand, ".",
                                 render(OrderedDict(ID=f"{tid}.exon{k}", Parent=tid))]) + "\n")
        for x in sorted(m["cds"], key=lambda x: int(x[0][3])):
            gff.write("\t".join([sid, "RefSeq", "CDS", x[0][3], x[0][4], ".", strand, x[0][7],
                                 render(OrderedDict(ID=f"cds-{tid}", Parent=tid, protein_id=m["protein_id"] or tid,
                                                    product=m["product"]))]) + "\n")
        hdr = f"{tid} locus_tag={lt} product={m['product']}"
        if m["protein_id"]: hdr += f" refseq={m['protein_id']}"
        if m["p5"] or m["p3"]:
            hdr += f" partial={'5' if m['p5'] else ''}{'3' if m['p3'] else ''}"
        faa.write(f">{hdr}\n{wrap(m['aa'])}\n"); cdsf.write(f">{hdr}\n{wrap(m['nt'])}\n")
        if i == 1: lng.write(f">{tid}\n{wrap(m['aa'])}\n")
        mp.write("\t".join(map(str, [m["gid"], m["tx_rs"], lt, tid, sid, s, e, strand, len(m["aa"]),
                                     int(m["p5"]), int(m["p3"]), "NA", m["protein_id"], m["product"]])) + "\n")
for fh in (gff, faa, lng, cdsf, mp): fh.close()
with open(f"{a.out}.removed.tsv", "w") as fo:
    fo.write("old_gene_id\told_transcript_id\tseqid\tstart\tend\tstrand\treason\n")
    for r in removed: fo.write("\t".join(map(str, r)) + "\n")
with open(f"{a.out}.stats.tsv", "w") as fo:
    for k in ("input_genes", "genes", "transcripts", "partial_5prime", "partial_3prime",
              "removed_internal_stop", "removed_genes_no_valid_model"):
        fo.write(f"{k}\t{stats[k]}\n")
print(f"Genome check: {sum(1 for s in used if seqmap.get(s))}/{len(used)} RefSeq sequences with genes identical to the project genome")
print("Official annotation: " + ", ".join(f"{k}={stats[k]}" for k in ("input_genes", "genes", "transcripts",
      "partial_5prime", "partial_3prime", "removed_internal_stop")))
