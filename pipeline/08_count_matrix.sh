#!/usr/bin/env bash
# ============================================================
# 08_count_matrix.sh <species>              GOAL 5: count matrices
#
# Counts read pairs per gene (locus_tag) with featureCounts on the BAMs of
# step 03 and the GTF of step 06, one matrix per species.
#  - library strandedness detected automatically (unstranded / forward / reverse)
#  - paired-end: fragments counted once, both mates mapped (-B), same
#    chromosome (-C); only uniquely placed reads (MAPQ >= 10, primary).
#    For A. latus this keeps reads that cannot be assigned to one homeolog
#    out of the counts; they are reported as 'Unassigned_MultiMapping'.
#
# Needs: steps 03 and 06 for the species (step 07 optional: adds orthologue
#        columns to the annotated matrix).
# Output: final_results/08_count_matrices/<sp>/
#   <sp>.counts.tsv            raw counts, genes x samples (input for DESeq2/edgeR)
#   <sp>.TPM.tsv  <sp>.CPM.tsv normalised expression
#   <sp>.counts.annotated.tsv  counts + subgenome + orthologues
#   <sp>.sample_metadata.tsv   sample, group, replicate
#   <sp>.mapping_stats.tsv     reads, alignment, assignment per sample
#   <sp>.featureCounts.txt(.summary), strandedness.txt
#   figures/Fig_RNAseq_QC_<sp>.{png,pdf} (+ single panels in figures/RNAseq_<sp>/)
# Time: ~10-30 min.
# ============================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"; source "$(dirname "$0")/lib/common.sh"
species_arg "${RNASEQ_SPECIES[@]}" -- "$@"
activate "${ENV_TOOLS}"
banner "08 count matrix: ${SP}"
[[ -s "${SAMPLES}" ]] || die "sample sheet not found: ${SAMPLES}"

R="$(ref_dir "${SP}")"; GTF="${R}/${SP}.gtf"; need "${GTF}"
C="${D08}/${SP}"; mkdir -p "${C}"
BAMS=()
while IFS=$'\t' read -r s r1 r2 grp rep; do
    need "${D03}/${SP}/bam/${s}.bam"; BAMS+=("${D03}/${SP}/bam/${s}.bam")
done < <(sheet_samples "${SP}")
echo "${#BAMS[@]} samples"
FC=(-p --countReadPairs -B -C -Q 10 --primary -t exon -g gene_id)

# ---------- strandedness (first sample) ----------
if [[ ! -s "${C}/strandedness.txt" ]]; then
    for s in 0 1 2; do
        featureCounts -T 8 "${FC[@]}" -s "${s}" -a "${GTF}" -o "${C}/strand_s${s}.txt" "${BAMS[0]}" > /dev/null 2>&1
    done
    python3 - "${C}" > "${C}/strandedness.txt" <<'PY'
import sys
a = {s: int(next(l for l in open(f"{sys.argv[1]}/strand_s{s}.txt.summary") if l.startswith("Assigned")).split()[1]) for s in (0, 1, 2)}
print(1 if a[1] > 3 * a[2] else 2 if a[2] > 3 * a[1] else 0)
PY
    rm -f "${C}"/strand_s*.txt*
fi
STRAND="$(cat "${C}/strandedness.txt")"
echo "Library: $( case ${STRAND} in 0) echo unstranded;; 1) echo stranded, forward;; 2) echo stranded, reverse;; esac ) (-s ${STRAND})"

featureCounts -T "${THREADS}" "${FC[@]}" -s "${STRAND}" -a "${GTF}" -o "${C}/${SP}.featureCounts.txt" "${BAMS[@]}" \
    > "${C}/featureCounts.log" 2>&1 || { tail -20 "${C}/featureCounts.log"; die "featureCounts failed"; }

# ---------- gene annotation columns for the annotated matrix ----------
python3 - "${R}/${SP}.id_map.tsv" "${D07}/tables/correspondence_${SP}.tsv" "${C}/${SP}.gene_info.tsv" <<'PY'
import csv, os, sys
idmap, corr, out = sys.argv[1:]
genes = {}
for r in csv.DictReader(open(idmap), delimiter="\t"):
    genes.setdefault(r["locus_tag"], {"locus_tag": r["locus_tag"], "contig": r["seqid"], "start": r["start"],
                                      "end": r["end"], "strand": r["strand"], "subgenome": r.get("subgenome", "NA")})
cols = ["locus_tag", "contig", "start", "end", "strand", "subgenome"]
if os.path.exists(corr):
    rd = csv.DictReader(open(corr), delimiter="\t")
    oc = [c for c in rd.fieldnames if c.endswith("_orthologs")]
    for r in rd:
        if r["locus_tag"] in genes:
            for c in oc: genes[r["locus_tag"]][c[:-1]] = r[c]
            if "homeolog" in r: genes[r["locus_tag"]]["homeolog"] = r["homeolog"]
    cols += [c[:-1] for c in oc] + (["homeolog"] if "homeolog" in rd.fieldnames else [])
with open(out, "w", newline="") as fo:
    w = csv.DictWriter(fo, cols, delimiter="\t", extrasaction="ignore", restval=""); w.writeheader()
    w.writerows(genes.values())
PY

python3 "${BIN}/make_count_matrix.py" --featurecounts "${C}/${SP}.featureCounts.txt" \
    --sheet "${SAMPLES}" --species "${SP}" --qc-dir "${D03}/${SP}/qc" \
    --annotation "${C}/${SP}.gene_info.tsv" --out "${C}/${SP}"

python3 "${BIN}/rnaseq_figures.py" --prefix "${C}/${SP}" --species "${SP}" --outdir "${FIG}/RNAseq_${SP}"
for ext in png pdf; do cp -f "${FIG}/RNAseq_${SP}/QC_panel.${ext}" "${FIG}/Fig_RNAseq_QC_${SP}.${ext}"; done
mark_done "${C}"
echo
column -t -s $'\t' "${C}/${SP}.mapping_stats.tsv" | cut -c1-200
echo "Count matrix: ${C}/${SP}.counts.tsv     Figure: ${FIG}/Fig_RNAseq_QC_${SP}.png"
