#!/usr/bin/env python3
"""
curate_annotation.py

Converts raw BRAKER3 / TSEBRA gene models (GFF3) into a curated GFF3 with
final locus_tags, ready for NCBI table2asn.

Steps:
  1. Reads the soft-masked genome and the raw GFF3.
  2. Rebuilds every transcript (exons, CDS, UTRs), translates the CDS and classifies it:
       - internal stop codon                      -> removed
       - protein < --min-aa                       -> removed
       - CDS masked fraction > --max-masked (likely TE) -> removed
       - complete CDS with length not a multiple of 3   -> removed
       - duplicated transcript (same CDS structure)     -> removed
       - no ATG / no stop codon                   -> kept, flagged as partial (NCBI)
  3. Sorts genes by contig order in the FASTA (not alphabetical), then by
     start, end and strand.
  4. Assigns locus_tags PREFIX_000010, PREFIX_000020 ... (steps of 10, as NCBI
     recommends, so genes can be inserted later without renumbering).
     Transcripts: PREFIX_000010-T1 (T1 = longest protein), -T2 ...
  5. Renames ALL features (exon, CDS, UTR); no original ID survives.
  6. Adds protein_id / transcript_id in gnl|DB|... format, product
     ("hypothetical protein" until functional annotation) and partial-gene
     markers (partial, start_range, end_range).
  7. Optional: adds the subgenome (BED from phase_subgenomes.py) to every gene.

Outputs (with --out OUTPUT_PREFIX):
  .gff3  .proteins.faa  .proteins.longest.faa  .cds.fna
  .id_map.tsv  .removed.tsv  .stats.tsv

Uses only the Python standard library.
"""
import argparse
import bisect
import re
import sys
from collections import defaultdict, OrderedDict
from urllib.parse import unquote

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTacgtNnRYKMSWrykmswBDHVbdhv", "TGCAtgcaNnYRMKSWyrmkswVHDBvhdb")

TX_TYPES = {"mRNA", "transcript"}
UTR5 = {"five_prime_UTR", "5'-UTR", "five_prime_utr"}
UTR3 = {"three_prime_UTR", "3'-UTR", "three_prime_utr"}


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3].upper(), "X") for i in range(0, len(nt) - 2, 3))


def read_fasta(path):
    seqs = OrderedDict()
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def parse_attrs(s):
    d = OrderedDict()
    for item in s.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
        else:
            k, v = item, ""
        d[k.strip()] = unquote(v.strip())
    return d


def esc(v):
    # GFF3: escapar caracteres reservados nos valores
    return (str(v).replace("%", "%25").replace(";", "%3B").replace("=", "%3D")
            .replace("&", "%26").replace(",", "%2C").replace("\t", "%09"))


def render(attrs):
    return ";".join(f"{k}={esc(v)}" for k, v in attrs.items() if v is not None and v != "")


def merge_intervals(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def wrap(seq, n=60):
    return "\n".join(seq[i:i + n] for i in range(0, len(seq), n))


class Tx:
    __slots__ = ("old_id", "gene", "seqid", "strand", "start", "end", "exons",
                 "cds", "utr5", "utr3", "source", "nt", "aa", "partial5",
                 "partial3", "masked_frac", "first_phase", "new_id")

    def __init__(self, old_id, gene, seqid, strand, start, end, source):
        self.old_id, self.gene = old_id, gene
        self.seqid, self.strand = seqid, strand
        self.start, self.end, self.source = start, end, source
        self.exons, self.cds, self.utr5, self.utr3 = [], [], [], []


def load_gff(path):
    genes = OrderedDict()   # old_gene -> dict(seqid,strand,source,txs=[])
    txs = OrderedDict()
    children = defaultdict(list)  # parent tx -> [(type,s,e,phase)]
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) != 9:
                continue
            seqid, source, ftype, s, e, _, strand, phase, attr = f
            a = parse_attrs(attr)
            s, e = int(s), int(e)
            if ftype == "gene":
                gid = a.get("ID") or sys.exit(f"Gene sem ID (linha {ln})")
                genes.setdefault(gid, {"seqid": seqid, "strand": strand, "source": source})
            elif ftype in TX_TYPES:
                tid = a.get("ID") or sys.exit(f"Transcrito sem ID (linha {ln})")
                gid = a.get("Parent", "").split(",")[0] or tid.rsplit(".t", 1)[0]
                txs[tid] = Tx(tid, gid, seqid, strand, s, e, source)
                genes.setdefault(gid, {"seqid": seqid, "strand": strand, "source": source})
            else:
                for p in a.get("Parent", "").split(","):
                    if p:
                        children[p].append((ftype, s, e, phase))
    for tid, items in children.items():
        if tid not in txs:
            continue  # children of genes (e.g. orphan introns) are ignored
        t = txs[tid]
        for ftype, s, e, ph in items:
            if ftype == "exon":
                t.exons.append((s, e))
            elif ftype == "CDS":
                t.cds.append((s, e))
            elif ftype in UTR5:
                t.utr5.append((s, e))
            elif ftype in UTR3:
                t.utr3.append((s, e))
    return genes, txs


def build_tx(t, genome):
    seq = genome[t.seqid]
    t.cds = merge_intervals(t.cds)
    if not t.exons:
        t.exons = merge_intervals(t.cds + t.utr5 + t.utr3)
    else:
        t.exons = merge_intervals(t.exons)
    if not t.cds:
        return "no_CDS"
    t.start = min(s for s, _ in t.exons)
    t.end = max(e for _, e in t.exons)
    parts = [seq[s - 1:e] for s, e in t.cds]
    nt = "".join(parts)
    if t.strand == "-":
        nt = revcomp(nt)
    elif t.strand != "+":
        return "no_strand"
    masked = sum(1 for c in nt if c.islower())
    t.masked_frac = masked / len(nt) if nt else 0
    nt = nt.upper()
    t.partial5 = not nt.startswith("ATG")
    t.partial3 = nt[-3:] not in STOPS
    # for a 5'-partial model, the reading frame is the one giving the longest ORF without stop
    offset = 0
    if t.partial5:
        best = None
        for o in (0, 1, 2):
            p = translate(nt[o:])
            body = p[:-1] if p.endswith("*") else p
            if "*" not in body and (best is None or len(p) > best[1]):
                best = (o, len(p))
        offset = best[0] if best else 0
    t.first_phase = offset
    t.nt = nt
    prot = translate(nt[offset:])
    if not t.partial3 and prot.endswith("*"):
        prot = prot[:-1]
    if "*" in prot:
        return "internal_stop"
    if not t.partial5 and not t.partial3 and len(nt) % 3 != 0:
        return "frameshift_len_not_mod3"
    t.aa = prot
    return None


def cds_phases(t):
    """GFF3 phase of every CDS segment, computed in the 5'->3' direction."""
    segs = sorted(t.cds, reverse=(t.strand == "-"))
    phases = {}
    done = 0  # nt before this segment, after the initial offset
    for i, (s, e) in enumerate(segs):
        if i == 0:
            phases[(s, e)] = t.first_phase
            done = (e - s + 1) - t.first_phase
        else:
            phases[(s, e)] = (3 - done % 3) % 3
            done += e - s + 1
    return phases


def load_bed(path):
    d = defaultdict(list)
    if not path:
        return d
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.rstrip().split("\t")
            d[f[0]].append((int(f[1]), int(f[2]), f[3]))
    for k in d:
        d[k].sort()
    return d


def bed_label(bed, seqid, pos):
    iv = bed.get(seqid)
    if not iv:
        return "unassigned"
    i = bisect.bisect_right([x[0] for x in iv], pos - 1) - 1
    if i >= 0 and iv[i][0] <= pos - 1 < iv[i][1]:
        return iv[i][2]
    return "unassigned"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", required=True, help="soft-masked genome used for gene prediction")
    ap.add_argument("--prefix", required=True, help="locus_tag prefix registered at NCBI")
    ap.add_argument("--out", required=True)
    ap.add_argument("--dbname", default=None, help="tag of gnl|DB|... (default: prefix)")
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--width", type=int, default=6)
    ap.add_argument("--min-aa", type=int, default=50)
    ap.add_argument("--max-masked", type=float, default=0.50)
    ap.add_argument("--subgenome-bed", default=None)
    ap.add_argument("--source", default=None, help="column 2 of the final GFF3 (default: keep)")
    args = ap.parse_args()

    if "_" in args.prefix:
        sys.exit("ERROR: the NCBI locus_tag prefix cannot contain '_' (e.g. ALAT, not A_LAT).")
    dbname = args.dbname or args.prefix

    genome = read_fasta(args.fasta)
    order = {k: i for i, k in enumerate(genome)}
    genes, txs = load_gff(args.gff)
    bed = load_bed(args.subgenome_bed)

    removed = []
    keep = defaultdict(list)
    for tid, t in txs.items():
        if t.seqid not in genome:
            sys.exit(f"ERROR: contig {t.seqid} of the GFF is not in the FASTA")
        why = build_tx(t, genome)
        if why is None and len(t.aa) < args.min_aa:
            why = f"short_protein_{len(t.aa)}aa"
        if why is None and t.masked_frac > args.max_masked:
            why = f"repeat_masked_{t.masked_frac:.2f}"
        if why:
            removed.append((t.gene, tid, t.seqid, t.start, t.end, t.strand, why))
        else:
            keep[t.gene].append(t)

    # remover isoformas redundantes (mesma CDS)
    for g, lst in keep.items():
        seen = {}
        uniq = []
        for t in sorted(lst, key=lambda x: (-len(x.aa), x.old_id)):
            key = (tuple(t.cds), t.strand)
            if key in seen:
                removed.append((g, t.old_id, t.seqid, t.start, t.end, t.strand,
                                f"duplicate_of_{seen[key]}"))
                continue
            seen[key] = t.old_id
            uniq.append(t)
        keep[g] = uniq   # already sorted: T1 = longest protein

    gene_recs = []
    for g, lst in keep.items():
        if not lst:
            continue
        gene_recs.append((order[lst[0].seqid], min(t.start for t in lst),
                          max(t.end for t in lst), lst[0].strand, g))
    gene_recs.sort()

    out_gff = open(f"{args.out}.gff3", "w")
    out_faa = open(f"{args.out}.proteins.faa", "w")
    out_long = open(f"{args.out}.proteins.longest.faa", "w")
    out_cds = open(f"{args.out}.cds.fna", "w")
    out_map = open(f"{args.out}.id_map.tsv", "w")
    out_map.write("old_gene_id\told_transcript_id\tlocus_tag\ttranscript_id\tseqid\tstart\tend\tstrand\tprotein_len\tpartial5\tpartial3\tsubgenome\n")

    out_gff.write("##gff-version 3\n")
    used = {r[0] for r in gene_recs}
    inv = {v: k for k, v in order.items()}
    for i in sorted(used):
        out_gff.write(f"##sequence-region {inv[i]} 1 {len(genome[inv[i]])}\n")

    stats = defaultdict(int)
    for n, (_, gs, ge, strand, g) in enumerate(gene_recs, start=1):
        lt = f"{args.prefix}_{n * args.step:0{args.width}d}"
        lst = keep[g]
        seqid = lst[0].seqid
        src = args.source or genes.get(g, {}).get("source", "BRAKER3")
        sub = bed_label(bed, seqid, (gs + ge) // 2) if bed else None
        if sub:
            stats[f"genes_subgenome_{sub}"] += 1
        stats["genes"] += 1
        g_p5 = any(t.partial5 for t in lst)
        g_p3 = any(t.partial3 for t in lst)
        ga = OrderedDict(ID=lt, locus_tag=lt)
        if sub:
            ga["subgenome"] = sub
        add_partial(ga, g_p5, g_p3, strand, gs, ge)
        out_gff.write("\t".join([seqid, src, "gene", str(gs), str(ge), ".", strand, ".", render(ga)]) + "\n")

        for k, t in enumerate(lst, start=1):
            tid = f"{lt}-T{k}"
            t.new_id = tid
            stats["transcripts"] += 1
            stats["partial_5prime"] += t.partial5
            stats["partial_3prime"] += t.partial3
            ta = OrderedDict(ID=tid, Parent=lt, locus_tag=lt,
                             product="hypothetical protein",
                             transcript_id=f"gnl|{dbname}|mrna.{tid}",
                             protein_id=f"gnl|{dbname}|{tid}")
            add_partial(ta, t.partial5, t.partial3, strand, t.start, t.end)
            out_gff.write("\t".join([seqid, src, "mRNA", str(t.start), str(t.end), ".", strand, ".", render(ta)]) + "\n")

            ex = sorted(t.exons, reverse=(strand == "-"))
            for j, (s, e) in enumerate(ex, start=1):
                ea = OrderedDict(ID=f"{tid}.exon{j}", Parent=tid, locus_tag=lt)
                out_gff.write("\t".join([seqid, src, "exon", str(s), str(e), ".", strand, ".", render(ea)]) + "\n")
            ph = cds_phases(t)
            for s, e in sorted(t.cds, reverse=(strand == "-")):
                ca = OrderedDict(ID=f"{tid}.cds", Parent=tid, locus_tag=lt,
                                 product="hypothetical protein",
                                 protein_id=f"gnl|{dbname}|{tid}",
                                 transcript_id=f"gnl|{dbname}|mrna.{tid}")
                out_gff.write("\t".join([seqid, src, "CDS", str(s), str(e), ".", strand, str(ph[(s, e)]), render(ca)]) + "\n")
            for lab, iv in (("five_prime_UTR", t.utr5), ("three_prime_UTR", t.utr3)):
                for j, (s, e) in enumerate(sorted(merge_intervals(iv), reverse=(strand == "-")), start=1):
                    ua = OrderedDict(ID=f"{tid}.{'utr5p' if lab.startswith('five') else 'utr3p'}{j}", Parent=tid)
                    out_gff.write("\t".join([seqid, src, lab, str(s), str(e), ".", strand, ".", render(ua)]) + "\n")

            hdr = f"{tid} locus_tag={lt} gene={lt} len={len(t.aa)}"
            if t.partial5 or t.partial3:
                hdr += f" partial={'5' if t.partial5 else ''}{'3' if t.partial3 else ''}"
            if sub:
                hdr += f" subgenome={sub}"
            out_faa.write(f">{hdr}\n{wrap(t.aa)}\n")
            out_cds.write(f">{hdr}\n{wrap(t.nt)}\n")
            if k == 1:
                out_long.write(f">{tid}\n{wrap(t.aa)}\n")
            out_map.write("\t".join(map(str, [g, t.old_id, lt, tid, seqid, t.start, t.end, strand,
                                              len(t.aa), int(t.partial5), int(t.partial3), sub or "NA"])) + "\n")

    for fh in (out_gff, out_faa, out_long, out_cds, out_map):
        fh.close()

    with open(f"{args.out}.removed.tsv", "w") as fo:
        fo.write("old_gene_id\told_transcript_id\tseqid\tstart\tend\tstrand\treason\n")
        for r in removed:
            fo.write("\t".join(map(str, r)) + "\n")
            stats["removed_" + re.sub(r"_(of_.*|[\d.]+(aa)?)$", "", r[-1])] += 1

    genes_lost = len({g for g in genes}) - stats["genes"]
    stats["input_genes"] = len(genes)
    stats["input_transcripts"] = len(txs)
    stats["genes_removed_entirely"] = genes_lost
    with open(f"{args.out}.stats.tsv", "w") as fo:
        for k in sorted(stats):
            fo.write(f"{k}\t{stats[k]}\n")
    for k in sorted(stats):
        print(f"{k:32s} {stats[k]}")


def add_partial(a, p5, p3, strand, s, e):
    if not (p5 or p3):
        return
    a["partial"] = "true"
    # start_range/end_range in genome coordinates (NCBI GFF3)
    left = p5 if strand == "+" else p3
    right = p3 if strand == "+" else p5
    if left:
        a["start_range"] = f".,{s}"
    if right:
        a["end_range"] = f"{e},."


if __name__ == "__main__":
    main()
