#!/usr/bin/env python3
"""
make_count_matrix.py - count matrices and mapping statistics for one species.

Inputs
  --featurecounts  featureCounts output (gene_id = locus_tag; .summary next to it)
  --sheet          rnaseq_samples.tsv: species sample R1 R2 strain condition replicate
  --species        species name as in the sheet
  --qc-dir         directory with <sample>.fastp.json and <sample>.hisat2.summary
  --annotation     optional <sp>.gene_annotation.tsv (adds product/GO-free columns)
  --out            output prefix

Outputs
  PREFIX.counts.tsv            raw counts, genes x samples (input for DESeq2/edgeR)
  PREFIX.TPM.tsv / .CPM.tsv    normalised expression
  PREFIX.counts.annotated.tsv  counts + product, Af293 orthologue, subgenome
  PREFIX.sample_metadata.tsv   colData for DESeq2
  PREFIX.mapping_stats.tsv     per-sample read / alignment / assignment table
"""
import argparse
import csv
import json
import os
import re


def parse_hisat(path):
    d = {}
    if not os.path.exists(path):
        return d
    txt = open(path).read()
    m = re.search(r"^(\d+) reads; of these:", txt, re.M)
    d["hisat2_input"] = int(m.group(1)) if m else None
    m = re.search(r"(\d+) \([\d.]+%\) aligned concordantly exactly 1 time", txt) or \
        re.search(r"(\d+) \([\d.]+%\) aligned exactly 1 time", txt)
    d["unique"] = int(m.group(1)) if m else None
    m = re.search(r"(\d+) \([\d.]+%\) aligned concordantly >1 times", txt) or \
        re.search(r"(\d+) \([\d.]+%\) aligned >1 times", txt)
    d["multi"] = int(m.group(1)) if m else None
    m = re.search(r"([\d.]+)% overall alignment rate", txt)
    d["overall_rate"] = float(m.group(1)) if m else None
    return d


def parse_fastp(path):
    if not os.path.exists(path):
        return {}
    j = json.load(open(path))
    s = j["summary"]
    return {"raw_reads": s["before_filtering"]["total_reads"],
            "clean_reads": s["after_filtering"]["total_reads"],
            "q30_after": round(100 * s["after_filtering"]["q30_rate"], 2),
            "gc_after": round(100 * s["after_filtering"]["gc_content"], 2),
            "dup_rate": round(100 * j.get("duplication", {}).get("rate", 0), 2)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--featurecounts", required=True)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--qc-dir", required=True)
    ap.add_argument("--annotation")
    ap.add_argument("--min-count", type=int, default=10)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # sample sheet: columns found by the header line ("# species sample R1 R2 strain condition replicate")
    meta, order, cols = {}, [], None
    for line in open(a.sheet):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            if cols is None and "sample" in [x.strip("# ") for x in f]:
                cols = [x.strip("# ") for x in f]
            continue
        if not f or f[0] != a.species:
            continue
        r = dict(zip(cols or ["species", "sample", "R1", "R2", "condition", "replicate"], f))
        strain, cond = r.get("strain", "NA") or "NA", r.get("condition", "NA") or "NA"
        meta[r["sample"]] = {"strain": strain, "condition": cond,
                             "group": f"{strain}_{cond}" if strain != "NA" else cond,
                             "replicate": r.get("replicate", "NA"),
                             "layout": "PE" if r.get("R2") else "SE"}
        order.append(r["sample"])

    with open(a.featurecounts) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    hdr = lines[0].rstrip("\n").split("\t")
    bam_cols = hdr[6:]
    names = [re.sub(r"\.bam$", "", os.path.basename(b)) for b in bam_cols]
    idx = {n: i for i, n in enumerate(names)}
    samples = [s for s in order if s in idx] + [n for n in names if n not in order]
    genes, length, counts = [], [], []
    for l in lines[1:]:
        f = l.rstrip("\n").split("\t")
        genes.append(f[0])
        length.append(int(f[5]))
        row = [int(x) for x in f[6:]]
        counts.append([row[idx[s]] for s in samples])

    tot = [sum(c[j] for c in counts) for j in range(len(samples))]
    rpk = [[c[j] / (L / 1000) if L else 0 for j in range(len(samples))] for c, L in zip(counts, length)]
    rpk_tot = [sum(r[j] for r in rpk) for j in range(len(samples))]

    def write_matrix(path, fn, fmt):
        with open(path, "w") as fo:
            fo.write("locus_tag\t" + "\t".join(samples) + "\n")
            for i, g in enumerate(genes):
                fo.write(g + "\t" + "\t".join(fmt(fn(i, j)) for j in range(len(samples))) + "\n")

    write_matrix(f"{a.out}.counts.tsv", lambda i, j: counts[i][j], str)
    write_matrix(f"{a.out}.TPM.tsv", lambda i, j: rpk[i][j] / rpk_tot[j] * 1e6 if rpk_tot[j] else 0, lambda x: f"{x:.3f}")
    write_matrix(f"{a.out}.CPM.tsv", lambda i, j: counts[i][j] / tot[j] * 1e6 if tot[j] else 0, lambda x: f"{x:.3f}")

    if a.annotation and os.path.exists(a.annotation):
        ann = {}
        with open(a.annotation) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            keep = [c for c in rd.fieldnames if c in ("subgenome", "product", "eggNOG_preferred_name")
                    or c.endswith("_ortholog") or c == "homeolog"]
            for r in rd:
                ann[r["locus_tag"]] = [r[c] for c in keep]
        with open(f"{a.out}.counts.annotated.tsv", "w") as fo:
            fo.write("locus_tag\tlength\t" + "\t".join(keep) + "\t" + "\t".join(samples) + "\n")
            for i, g in enumerate(genes):
                fo.write(f"{g}\t{length[i]}\t" + "\t".join(ann.get(g, [""] * len(keep))) + "\t"
                         + "\t".join(map(str, counts[i])) + "\n")

    with open(f"{a.out}.sample_metadata.tsv", "w") as fo:
        fo.write("sample\tspecies\tstrain\tcondition\tgroup\treplicate\tlayout\n")
        for s in samples:
            m = meta.get(s, {"strain": "NA", "condition": "NA", "group": "NA", "replicate": "NA", "layout": "NA"})
            fo.write(f"{s}\t{a.species}\t{m['strain']}\t{m['condition']}\t{m['group']}\t{m['replicate']}\t{m['layout']}\n")

    summ = {}
    sp = a.featurecounts + ".summary"
    if os.path.exists(sp):
        with open(sp) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            sn = [re.sub(r"\.bam$", "", os.path.basename(x)) for x in h[1:]]
            for l in fh:
                f = l.rstrip("\n").split("\t")
                for k, v in zip(sn, f[1:]):
                    summ.setdefault(k, {})[f[0]] = int(v)

    cols = ["species", "sample", "condition", "layout", "raw_reads", "clean_reads", "pct_clean",
            "q30_after", "gc_after", "dup_rate", "aligned_unique", "aligned_multi",
            "overall_alignment_pct", "assigned_to_genes", "pct_assigned_of_clean",
            "unassigned_multimapping", "unassigned_no_feature", "unassigned_ambiguous",
            f"genes_detected_ge{a.min_count}", "pct_genes_detected"]
    with open(f"{a.out}.mapping_stats.tsv", "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for j, s in enumerate(samples):
            fp = parse_fastp(os.path.join(a.qc_dir, f"{s}.fastp.json"))
            hs = parse_hisat(os.path.join(a.qc_dir, f"{s}.hisat2.summary"))
            fc = summ.get(s, {})
            layout = meta.get(s, {}).get("layout", "NA")
            # fastp counts reads; HISAT2/featureCounts count pairs in PE mode
            clean_units = (fp.get("clean_reads", 0) // 2) if layout == "PE" else fp.get("clean_reads", 0)
            det = sum(1 for c in counts if c[j] >= a.min_count)

            def pct(x, y):
                return f"{100 * x / y:.2f}" if x is not None and y else "NA"
            row = [a.species, s, meta.get(s, {}).get("condition", "NA"), layout,
                   fp.get("raw_reads", "NA"), fp.get("clean_reads", "NA"),
                   pct(fp.get("clean_reads"), fp.get("raw_reads")),
                   fp.get("q30_after", "NA"), fp.get("gc_after", "NA"), fp.get("dup_rate", "NA"),
                   hs.get("unique", "NA"), hs.get("multi", "NA"), hs.get("overall_rate", "NA"),
                   fc.get("Assigned", tot[j]), pct(fc.get("Assigned", tot[j]), clean_units),
                   fc.get("Unassigned_MultiMapping", "NA"), fc.get("Unassigned_NoFeatures", "NA"),
                   fc.get("Unassigned_Ambiguity", "NA"), det, pct(det, len(genes))]
            fo.write("\t".join(map(str, row)) + "\n")
    print(f"{a.species}: {len(genes)} genes x {len(samples)} samples -> {a.out}.counts.tsv")


if __name__ == "__main__":
    main()
