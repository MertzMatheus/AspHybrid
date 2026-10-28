#!/usr/bin/env python3
"""
clean_genome.py - clean copy of a genome FASTA + assembly statistics.

  * ID = first word of the header, characters outside [A-Za-z0-9._-] -> "_",
    max 50 characters, made unique
  * sequence in upper case, 80 nt per line (the original is not modified)
  * header map: new_id, length, original header

usage: clean_genome.py in.fasta[.gz] out.fasta header_map.tsv stats.tsv species strain
stats.tsv: one row (with header) - contigs, size, N50, L50, N90, largest, GC,
N, soft-masked fraction of the original, contigs < 200 bp, non-IUPAC
"""
import gzip
import os
import re
import sys

src, dst, mp, stats, sp, strain = sys.argv[1:7]
op = gzip.open if src.endswith(".gz") else open
iupac = set("ACGTNRYKMSWBDHV")
seen, lens = set(), []
gc = n = low = bad = 0


def flush(name, hdr, seq, fo, fm):
    global gc, n, low, bad
    s = "".join(seq)
    if not s:
        print(f"WARNING: empty sequence skipped: {hdr}", file=sys.stderr)
        return
    low += sum(c.islower() for c in s)
    u = s.upper()
    gc += u.count("G") + u.count("C")
    n += u.count("N")
    bad += sum(1 for c in u if c not in iupac)
    lens.append(len(u))
    fo.write(f">{name}\n")
    for i in range(0, len(u), 80):
        fo.write(u[i:i + 80] + "\n")
    fm.write(f"{name}\t{len(u)}\t{hdr}\n")


with op(src, "rt") as fi, open(dst, "w") as fo, open(mp, "w") as fm:
    fm.write("new_id\tlength\toriginal_header\n")
    name = hdr = None
    seq = []
    for line in fi:
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if name:
                flush(name, hdr, seq, fo, fm)
            hdr = line[1:].strip()
            base = re.sub(r"[^A-Za-z0-9._-]", "_", hdr.split()[0] if hdr else "seq")[:50]
            name, i = base, 2
            while name in seen:
                name = f"{base}_{i}"
                i += 1
            seen.add(name)
            seq = []
        elif line.strip():
            seq.append(line.strip())
    if name:
        flush(name, hdr, seq, fo, fm)

if not lens:
    sys.exit(f"ERROR: no sequences in {src}")
tot = sum(lens)
srt = sorted(lens, reverse=True)


def nx(frac):
    run = 0
    for k, L in enumerate(srt, 1):
        run += L
        if run >= tot * frac:
            return L, k
    return srt[-1], len(srt)


n50, l50 = nx(0.5)
n90, _ = nx(0.9)
acgt = tot - n
cols = ["species", "strain", "fasta", "contigs", "total_bp", "N50", "L50", "N90", "largest",
        "GC_percent", "N_bases", "N_percent", "softmasked_frac_original", "contigs_lt200bp", "non_IUPAC"]
vals = [sp, strain, os.path.basename(src), len(lens), tot, n50, l50, n90, srt[0],
        f"{100 * gc / acgt:.2f}" if acgt else "NA", n, f"{100 * n / tot:.3f}", f"{low / tot:.3f}",
        sum(L < 200 for L in lens), bad]
with open(stats, "w") as fo:
    fo.write("\t".join(cols) + "\n" + "\t".join(map(str, vals)) + "\n")
print(f"{sp}: {len(lens)} contigs, {tot / 1e6:.2f} Mb, N50 {n50 / 1e3:.1f} kb (L50 {l50}), "
      f"GC {vals[9]}%, N {vals[11]}%")
if bad:
    print(f"WARNING: {bad} non-IUPAC characters", file=sys.stderr)
