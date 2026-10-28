#!/usr/bin/env bash
# ============================================================
# 11_differential_expression.sh <species>   DEGs + enrichment + figures + tables
#
# For every strain of the species: fungus + macrophage vs fungus alone in DMEN
# (control), from the design table METADATA (config/RNAseq_metadata.csv).
#  - DESeq2: one model per species, design ~ strain_treatment, Wald test,
#    apeglm-shrunken log2FC; DEG = padj < DE_PADJ and |log2FC| >= DE_LFC
#  - enrichment of up- and down-regulated DEGs: GO BP/MF/CC, KEGG pathways,
#    COG (eggNOG annotation of step 10), hypergeometric test + BH
#  - figures: PCA, DEG bar graph, volcano per strain, Venn of DEGs shared by
#    the strains (up and down), enrichment dot plot
#  - Supplementary workbook per species
#
# Needs: steps 08 (counts) and 10 (eggNOG) for the species.
# Output: final_results/11_differential_expression/<sp>/
#   deseq2/  DE_<strain>.tsv, vst.tsv, normalized_counts.tsv, size_factors.tsv
#   DE_<strain>.annotated.tsv   all tested genes + eggNOG annotation + orthologues
#   DEGs_<strain>.tsv           significant genes only
#   DE_summary.tsv, DEG_sharing.tsv, DEG_strain_membership.tsv
#   enrichment_significant.tsv, enrichment_all_tests.tsv, enrichment_summary.tsv
#   Supplementary_DE_<sp>.xlsx
# figures/Fig_DE_<sp>.{png,pdf}, Fig_enrichment_<sp>.{png,pdf}
# Time: ~5-15 min per species.
# ============================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"; source "$(dirname "$0")/lib/common.sh"
# (also defined in lib/common.sh; kept here so an older common.sh still works)
D10="${D10:-${OUT}/10_functional_annotation}"; D11="${D11:-${OUT}/11_differential_expression}"
species_arg "${RNASEQ_SPECIES[@]}" -- "$@"
banner "11 differential expression: ${SP} (each strain: Macrophage vs DMEN)"

CNT="${D08}/${SP}/${SP}.counts.tsv"; need "${CNT}" "${METADATA}"
E="${D10}/${SP}"
[[ -s "${E}/${SP}.gene_annotation.tsv" ]] || echo "WARNING: no eggNOG annotation yet (run 10_eggnog.sh ${SP}); DEGs without enrichment"
O="${D11}/${SP}"; mkdir -p "${O}/deseq2"

# metadata CSV -> check it covers every sample of the species
activate "${ENV_TOOLS}"
python3 - "${METADATA}" "${CNT}" "${SP}" <<'PY'
import csv, sys
meta, cnt, sp = sys.argv[1:]
rows = [r for r in csv.DictReader(open(meta, encoding="utf-8-sig")) if r["species_id"] == sp]
cols = open(cnt).readline().rstrip("\n").split("\t")[1:]
need = {"sample_old", "sample_id", "species_id", "strain", "treatment", "replicate"}
miss = need - set(rows[0].keys()) if rows else need
if miss: sys.exit(f"ERROR: metadata lacks columns {miss}")
names = {r["sample_old"] for r in rows}
if names - set(cols): sys.exit(f"ERROR: in metadata but not in counts: {sorted(names - set(cols))}")
if set(cols) - names: print(f"NOTE: count columns not in metadata (ignored): {sorted(set(cols) - names)}")
from collections import Counter
for (s, t), n in sorted(Counter((r["strain"], r["treatment"]) for r in rows).items()):
    print(f"  strain {s:>6}  {t:<10} {n} replicates")
PY

# ---------- DESeq2 (R, DEtools environment) ----------
activate "${ENV_DE}"
Rscript "${BIN}/deseq2_strains.R" "${CNT}" "${METADATA}" "${SP}" "${O}/deseq2" "${DE_MIN_COUNT}" \
    2> "${O}/deseq2/deseq2.log" || { tail -20 "${O}/deseq2/deseq2.log"; die "DESeq2 failed (${O}/deseq2/deseq2.log)"; }

# ---------- calls, enrichment, figures ----------
activate "${ENV_TOOLS}"
python3 "${BIN}/de_analysis.py" --species "${SP}" --de-dir "${O}/deseq2" \
    --annotation "${E}/${SP}.gene_annotation.tsv" --genesets "${E}/genesets" \
    --out "${O}" --fig-dir "${FIG}" --padj "${DE_PADJ}" --lfc "${DE_LFC}" \
    --enrich-padj "${ENRICH_PADJ}" --min-set "${ENRICH_MIN_SET}" --max-set "${ENRICH_MAX_SET}"

# ---------- supplementary workbook ----------
SHEETS=("README=${O}/README_tables.tsv" "Samples=${O}/deseq2/samples_used.tsv" "DE summary=${O}/DE_summary.tsv"
        "DEG sharing=${O}/DEG_sharing.tsv" "DEGs by strain=${O}/DEG_strain_membership.tsv")
for f in "${O}"/DEGs_*.tsv;  do st="$(basename "${f}" .tsv)"; SHEETS+=("${st}=${f}"); done
for f in "${O}"/DE_*.annotated.tsv; do st="$(basename "${f}" .annotated.tsv)"; SHEETS+=("${st} all genes=${f}"); done
SHEETS+=("Enrichment significant=${O}/enrichment_significant.tsv" "Enrichment summary=${O}/enrichment_summary.tsv"
         "Gene annotation=${E}/${SP}.gene_annotation.tsv")
python3 "${BIN}/tables_to_excel.py" "${O}/Supplementary_DE_${SP}.xlsx" "${SHEETS[@]}"
mark_done "${O}"
echo
cat "${O}/DE_summary.tsv" | cut -f2,6-9 | column -t
echo "Figures: ${FIG}/Fig_DE_${SP}.png  ${FIG}/Fig_enrichment_${SP}.png"
echo "Tables : ${O}/Supplementary_DE_${SP}.xlsx"
