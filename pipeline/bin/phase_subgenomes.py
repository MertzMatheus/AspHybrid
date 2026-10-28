#!/usr/bin/env python3
"""
phase_subgenomes.py  -  assigns every region of the hybrid genome (A. latus)
to one of the two parental subgenomes.

Logic (Steenwyk et al. 2020: homeologs differ by ~7 %):
  * the hybrid genome is cut into windows (default 5 kb);
  * each window is aligned (minimap2 -x asm20) separately to the genome of
    parent S (A. spinulosporus) and of parent Q (A. quadrilineatus, a close
    relative of the true parent);
  * identity = matches / aligned block length, of the best alignment;
  * window = S if idS - idQ >= delta; Q if idQ - idS >= delta;
    "ambiguous" if both are similar; "unassigned" if no good alignment;
  * consecutive windows with the same label are merged into blocks (BED);
    isolated ambiguous windows between two blocks of the same label are absorbed.

Subcommands:
  split  --fasta hybrid.fa --window 5000 --out windows.fa
  assign --windows-fai windows.fa.fai --paf-s S.paf --paf-q Q.paf --out prefix

Outputs of assign:
  prefix.windows.tsv   identity per window (for the histogram)
  prefix.blocks.bed    S / Q blocks (used by curate_annotation.py)
  prefix.contigs.tsv   summary per contig: S/Q fraction, subgenome switches
"""
import argparse
from collections import defaultdict


def cmd_split(a):
    def emit(name, seq, fo):
        L = len(seq)
        for s in range(0, L, a.window):
            e = min(s + a.window, L)
            if e - s < a.min_len:
                continue
            fo.write(f">{name}::{s}::{e}::{L}\n{seq[s:e]}\n")
    with open(a.fasta) as fi, open(a.out, "w") as fo:
        name, buf = None, []
        for line in fi:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    emit(name, "".join(buf), fo)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
        if name:
            emit(name, "".join(buf), fo)


def best_hits(paf, min_cov):
    best = {}
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            q, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            matches, block = int(f[9]), int(f[10])
            if block == 0:
                continue
            cov = (qe - qs) / qlen
            if cov < min_cov:
                continue
            ident = matches / block
            score = matches
            if q not in best or score > best[q][2]:
                best[q] = (ident, cov, score)
    return best


def cmd_assign(a):
    S = best_hits(a.paf_s, a.min_cov)
    Q = best_hits(a.paf_q, a.min_cov)
    wins = []
    with open(a.windows_fai) as fh:
        for line in fh:
            w = line.split("\t")[0]
            ctg, s, e, L = w.rsplit("::", 3)
            wins.append((ctg, int(s), int(e), int(L), w))

    rows = []
    for ctg, s, e, L, w in wins:
        iS = S.get(w, (0, 0, 0))[0]
        iQ = Q.get(w, (0, 0, 0))[0]
        if iS == 0 and iQ == 0:
            lab = "unassigned"
        elif iS - iQ >= a.delta:
            lab = "S"
        elif iQ - iS >= a.delta:
            lab = "Q"
        else:
            lab = "ambiguous"
        rows.append([ctg, s, e, L, iS, iQ, lab])

    with open(f"{a.out}.windows.tsv", "w") as fo:
        fo.write("contig\tstart\tend\tcontig_len\tident_S\tident_Q\tlabel\n")
        for r in rows:
            fo.write("\t".join([r[0], str(r[1]), str(r[2]), str(r[3]),
                                f"{r[4]:.4f}", f"{r[5]:.4f}", r[6]]) + "\n")

    by_ctg = defaultdict(list)
    for r in rows:
        by_ctg[r[0]].append(r)

    bed = open(f"{a.out}.blocks.bed", "w")
    summ = open(f"{a.out}.contigs.tsv", "w")
    summ.write("contig\tlength\tfrac_S\tfrac_Q\tfrac_ambiguous\tfrac_unassigned\tmajority\tn_switches\tflag\n")
    for ctg in sorted(by_ctg, key=lambda c: -by_ctg[c][0][3]):
        rs = sorted(by_ctg[ctg], key=lambda r: r[1])
        labs = [r[6] for r in rs]
        # absorb isolated ambiguous/unassigned windows between identical neighbours
        for i in range(1, len(labs) - 1):
            if labs[i] in ("ambiguous", "unassigned") and labs[i - 1] == labs[i + 1] and labs[i - 1] in "SQ":
                labs[i] = labs[i - 1]
        # short unlabelled contigs: take the majority label if there is one
        tot = defaultdict(int)
        for r, l in zip(rs, labs):
            tot[l] += r[2] - r[1]
        L = rs[0][3]
        known = {k: v for k, v in tot.items() if k in ("S", "Q")}
        maj = max(known, key=known.get) if known else "unassigned"
        # blocos
        blocks = []
        for r, l in zip(rs, labs):
            if blocks and blocks[-1][2] == l:
                blocks[-1][1] = r[2]
            else:
                blocks.append([r[1], r[2], l])
        if blocks:
            blocks[0][0] = 0
            blocks[-1][1] = L
        sq = [b[2] for b in blocks if b[2] in ("S", "Q")]
        switches = sum(1 for x, y in zip(sq, sq[1:]) if x != y)
        for s, e, l in blocks:
            bed.write(f"{ctg}\t{s}\t{e}\t{l}\n")
        flag = "OK"
        if switches:
            flag = "CHECK_switch(quimera_ou_recombinacao)"
        elif maj == "unassigned":
            flag = "CHECK_no_parent_hit"
        summ.write(f"{ctg}\t{L}\t{tot['S']/L:.3f}\t{tot['Q']/L:.3f}\t"
                   f"{tot['ambiguous']/L:.3f}\t{tot['unassigned']/L:.3f}\t{maj}\t{switches}\t{flag}\n")
    bed.close()
    summ.close()

    n = defaultdict(int)
    for r in rows:
        n[r[6]] += r[2] - r[1]
    total = sum(n.values())
    print("Genome fraction per label:")
    for k in ("S", "Q", "ambiguous", "unassigned"):
        print(f"  {k:12s} {n[k]:>12,d} bp  {100*n[k]/total:5.1f}%")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("split")
    p.add_argument("--fasta", required=True)
    p.add_argument("--window", type=int, default=5000)
    p.add_argument("--min-len", type=int, default=1000)
    p.add_argument("--out", required=True)
    p = sp.add_parser("assign")
    p.add_argument("--windows-fai", required=True)
    p.add_argument("--paf-s", required=True)
    p.add_argument("--paf-q", required=True)
    p.add_argument("--delta", type=float, default=0.02)
    p.add_argument("--min-cov", type=float, default=0.5)
    p.add_argument("--out", required=True)
    a = ap.parse_args()
    {"split": cmd_split, "assign": cmd_assign}[a.cmd](a)


if __name__ == "__main__":
    main()
