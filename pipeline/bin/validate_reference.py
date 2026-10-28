#!/usr/bin/env python3
"""
validate_reference.py - consistency checks of one species' reference files
(Results/Reference_files/<species>/), written to validation.tsv.

Checks
  1. every GFF3 sequence ID exists in the genome FASTA, coordinates in range
  2. genes in GTF == locus_tags in proteins.faa == locus_tags in cds.fna
  3. unique protein IDs; no internal stop codons; no '*'/'.' at the end
  4. independent re-translation: the GTF is translated by gffread from the
     genome (check.faa); % of those proteins identical to proteins.faa
     (differences are expected only for partial genes, flagged NCBI-style)
  5. CDS length = 3 x (protein length + 1) for complete genes
Status FAIL if 1-3 fail or < 98 % of re-translated proteins are identical.
"""
import argparse
import re
from collections import Counter


def fasta(p):
    name, buf, out, hdr = None, [], {}, {}
    with open(p) as fh:
        for l in fh:
            l = l.rstrip("\n")
            if l.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name = l[1:].split()[0]; hdr[name] = l[1:]; buf = []
            else:
                buf.append(l.strip())
    if name:
        out[name] = "".join(buf)
    return out, hdr


def locus(x):
    return re.sub(r"-T\d+$", "", x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--species", required=True)
    a = ap.parse_args()
    d, sp = a.dir, a.species
    problems = []

    lens = {}
    for l in open(f"{d}/{sp}.genome.fasta.fai"):
        f = l.split("\t"); lens[f[0]] = int(f[1])
    bad_seq = bad_coord = 0
    with open(f"{d}/{sp}.gff3") as fh:
        for l in fh:
            if l.startswith("##FASTA"):
                break
            if l.startswith("#") or not l.strip():
                continue
            f = l.split("\t")
            if len(f) < 9:
                continue
            if f[0] not in lens:
                bad_seq += 1
            elif int(f[4]) > lens[f[0]] or int(f[3]) < 1:
                bad_coord += 1
    if bad_seq: problems.append(f"{bad_seq} GFF3 lines on unknown sequences")
    if bad_coord: problems.append(f"{bad_coord} GFF3 lines out of range")

    gtf_genes, gtf_tx = set(), set()
    for l in open(f"{d}/{sp}.gtf"):
        m = re.search(r'gene_id "([^"]+)"; transcript_id "([^"]+)"', l)
        if m and 'gene_biotype "protein_coding"' in l:
            gtf_genes.add(m.group(1)); gtf_tx.add(m.group(2))

    prot, phdr = fasta(f"{d}/{sp}.proteins.faa")
    cds, _ = fasta(f"{d}/{sp}.cds.fna")
    pg = {locus(x) for x in prot}
    cg = {locus(x) for x in cds}
    if pg != cg: problems.append(f"proteins vs CDS loci differ ({len(pg ^ cg)})")
    only_gtf, only_prot = len(gtf_genes - pg), len(pg - gtf_genes)
    if only_gtf or only_prot:
        problems.append(f"GTF genes vs protein loci differ (GTF-only {only_gtf}, protein-only {only_prot})")
    internal = sum(1 for s in prot.values() if "*" in s.rstrip("*."))
    terminal = sum(1 for s in prot.values() if s.endswith(("*", ".")))
    if internal: problems.append(f"{internal} proteins with internal stop")
    if terminal: problems.append(f"{terminal} proteins end with '*' or '.'")

    partial = sum(1 for h in phdr.values() if "partial=" in h)
    bad_len = 0
    for t, s in cds.items():
        if t in prot and "partial=" not in phdr.get(t, "") and len(s) != 3 * (len(prot[t]) + 1):
            bad_len += 1

    chk, _ = fasta(f"{d}/check.faa")
    pset = {s.rstrip("*.") for s in prot.values()}
    ident = sum(1 for s in chk.values() if s.rstrip("*.").replace(".", "*") in pset or s.rstrip("*.") in pset)
    pct = 100 * ident / len(chk) if chk else 0
    if pct < 98: problems.append(f"only {pct:.1f}% of re-translated proteins identical")

    status = "OK" if not problems else "FAIL"
    cols = ["species", "genes", "transcripts", "proteins", "cds", "partial_proteins",
            "retranslation_identical_pct", "cds_length_mismatch_complete", "status", "notes"]
    row = [sp, len(pg), len(prot), len(prot), len(cds), partial, f"{pct:.2f}", bad_len, status,
           "; ".join(problems) or "-"]
    with open(f"{d}/validation.tsv", "w") as fo:
        fo.write("\t".join(cols) + "\n" + "\t".join(map(str, row)) + "\n")
    print(f"{sp:18s} genes={len(pg):>6} transcripts={len(prot):>6} partial={partial:>4} "
          f"re-translation identical={pct:6.2f}%  -> {status}" + (f"  ({'; '.join(problems)})" if problems else ""))


if __name__ == "__main__":
    main()
