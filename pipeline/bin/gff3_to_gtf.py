#!/usr/bin/env python3
"""
gff3_to_gtf.py - GFF3 -> GTF where gene_id IS the locus_tag.

Works for the curated project GFF3 (ALAT_000010 ...) and for RefSeq GFF3
(A. fumigatus Af293: gene_id becomes AFUA_1G00100, not "gene-AFUA_...").
This guarantees that featureCounts / Salmon count matrices use exactly the
same identifiers as the .gff3/.faa/.cds files and the annotation tables.

GTF attributes written: gene_id, transcript_id, gene_name (= locus_tag),
gene_biotype, product (when present), subgenome (A. latus).
Only transcripts with exons are written (mRNA, ncRNA, tRNA, rRNA ...).

Usage: gff3_to_gtf.py in.gff3 out.gtf [--protein-coding-only]
"""
import sys
from collections import OrderedDict, defaultdict
from urllib.parse import unquote


def attrs(s):
    d = {}
    for it in s.strip().strip(";").split(";"):
        if "=" in it:
            k, v = it.split("=", 1)
            d[k.strip()] = unquote(v.strip())
    return d


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    gin, gout = sys.argv[1], sys.argv[2]
    coding_only = "--protein-coding-only" in sys.argv

    feats = {}            # ID -> (type, attrs)
    parent = {}
    rows = []             # (seqid, src, type, s, e, strand, phase, parent_tx)
    order = OrderedDict()
    with open(gin) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            a = attrs(f[8])
            fid = a.get("ID")
            if fid:
                feats[fid] = (f[2], a, f)
                if "Parent" in a:
                    parent[fid] = a["Parent"].split(",")[0]
            if f[2] in ("exon", "CDS", "start_codon", "stop_codon"):
                for p in a.get("Parent", "").split(","):
                    if p:
                        rows.append((f, p))

    def gene_of(tx):
        g = parent.get(tx, tx)
        while g in parent and feats.get(g, ("",))[0] not in ("gene", "pseudogene"):
            g = parent[g]
        return g

    def locus(gid):
        a = feats.get(gid, (None, {}))[1]
        lt = a.get("locus_tag") or a.get("gene_id") or gid
        return lt[5:] if lt.startswith("gene-") else lt

    by_tx = defaultdict(list)
    for f, p in rows:
        by_tx[p].append(f)

    n_tx = n_gene = 0
    seen_genes = set()
    with open(gout, "w") as fo:
        for tx, lst in by_tx.items():
            ttype, ta, tf = feats.get(tx, ("mRNA", {}, None))
            gid = gene_of(tx)
            ga = feats.get(gid, (None, {}))[1]
            biotype = ga.get("gene_biotype") or ("protein_coding" if any(x[2] == "CDS" for x in lst) else ttype)
            if coding_only and biotype != "protein_coding":
                continue
            lt = locus(gid)
            tid = tx[4:] if tx.startswith("rna-") else tx
            extra = f'gene_name "{lt}"; gene_biotype "{biotype}";'
            prod = ta.get("product")
            if prod:
                extra += ' product "{}";'.format(prod.replace('"', "'"))
            sub = ga.get("subgenome")
            if sub:
                extra += f' subgenome "{sub}";'
            exons = [x for x in lst if x[2] == "exon"]
            cds = [x for x in lst if x[2] == "CDS"]
            if not exons:      # derive exons from CDS if absent
                exons = [x[:2] + ["exon"] + x[3:7] + ["."] + x[8:] for x in cds]
            if not exons:
                continue
            n_tx += 1
            if lt not in seen_genes:
                seen_genes.add(lt)
                n_gene += 1
            for x in sorted(exons, key=lambda r: int(r[3])):
                fo.write("\t".join([x[0], x[1], "exon", x[3], x[4], ".", x[6], ".",
                                    f'gene_id "{lt}"; transcript_id "{tid}"; {extra}']) + "\n")
            for x in sorted(cds, key=lambda r: int(r[3])):
                fo.write("\t".join([x[0], x[1], "CDS", x[3], x[4], ".", x[6], x[7],
                                    f'gene_id "{lt}"; transcript_id "{tid}"; {extra}']) + "\n")
    print(f"{gout}: {n_gene} genes, {n_tx} transcripts", file=sys.stderr)


if __name__ == "__main__":
    main()
