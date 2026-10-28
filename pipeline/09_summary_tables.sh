#!/usr/bin/env bash
# ============================================================
# 09_summary_tables.sh                      all tables for the paper (no argument)
#
# Collects every statistics table into one Excel workbook (one sheet per
# table) and lists the figures. Run any time; missing tables are skipped.
# Output: final_results/09_tables/Supplementary_Tables.xlsx
#         final_results/09_tables/software_versions.tsv
# ============================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"; source "$(dirname "$0")/lib/common.sh"
# (also defined in lib/common.sh; kept here so an older common.sh still works)
D10="${D10:-${OUT}/10_functional_annotation}"; D11="${D11:-${OUT}/11_differential_expression}"
activate "${ENV_TOOLS}"
banner "09 summary tables"
mkdir -p "${D09}"

# software versions (Methods section)
{
  printf "tool\tversion\n"
  printf "BUSCO\t%s\n"        "$(busco --version 2>&1 | head -1)"
  printf "BUSCO lineage\t%s\n" "$(basename "${BUSCO_LINEAGE_DIR}")"
  printf "HISAT2\t%s\n"       "$(hisat2 --version | grep -m1 -o 'version [0-9.]*')"
  printf "samtools\t%s\n"     "$(samtools --version | head -1)"
  printf "fastp\t%s\n"        "$(fastp --version 2>&1)"
  printf "featureCounts\t%s\n" "$(featureCounts -v 2>&1 | grep -o 'v[0-9.]*' | head -1)"
  printf "minimap2\t%s\n"     "$(minimap2 --version)"
  printf "gffread\t%s\n"      "$(gffread --version 2>&1)"
  printf "seqkit\t%s\n"       "$(seqkit version)"
  printf "BRAKER3 image\t%s\n" "${IMG_BRAKER}"
  printf "TE Tools image\t%s\n" "${IMG_TETOOLS}"
  printf "OrthoFinder\t%s\n"  "$(micromamba list -n "${ENV_ORTHO}" 2>/dev/null | awk '$1=="orthofinder"{print $2}')"
  printf "OrthoDB\t%s\n"      "$(basename "${ORTHODB_FUNGI}")"
  printf "eggNOG-mapper\t%s\n" "$(cat "${D10}"/*/emapper_version.txt 2>/dev/null | head -1)"
  printf "GO (go-basic.obo)\t%s\n" "$(grep -m1 '^data-version' "${GO_OBO}" 2>/dev/null | cut -d' ' -f2)"
  printf "DESeq2\t%s\n"       "$(head -2 "${D11}"/*/deseq2/deseq2_info.txt 2>/dev/null | grep -m1 DESeq2 | cut -d' ' -f2)"
  printf "log2FC shrinkage\t%s\n" "$(grep -h -m1 shrinkage "${D11}"/*/deseq2/deseq2_info.txt 2>/dev/null | head -1 | cut -d' ' -f2-)"
} > "${D09}/software_versions.tsv"

python3 "${BIN}/tables_to_excel.py" "${D09}/Supplementary_Tables.xlsx" \
  "S1 Genome QC=${D01}/genome_qc_summary.tsv" \
  "S2 Repeats=${D02}/*/*.repeat_summary.tsv" \
  "S3 A latus subgenomes=${D05}/subgenome_summary.tsv" \
  "S4 A latus contigs=${D05}/A_latus.contigs.tsv" \
  "S5 Annotation summary=${D06}/annotation_summary_all.tsv" \
  "S6 Validation=${D06}/*/validation.tsv" \
  "S7 Orthology statistics=${D07}/ortholog_statistics.tsv" \
  "S8 Orthologs A fumigatus=${D07}/Table_orthologs_Afumigatus.tsv" \
  "S8b Pairs Afum-Alat=${D07}/Table_pairs_Afumigatus_Alatus.tsv" \
  "S9 Gene sharing=${D07}/sharing/sharing_by_species.tsv" \
  "S10 Sharing patterns=${D07}/sharing/sharing_patterns.tsv" \
  "S11 Read mapping=${D03}/*/*.mapping_rates.tsv" \
  "S12 Counting stats=${D08}/*/*.mapping_stats.tsv" \
  "S13 Functional annotation=${D10}/functional_annotation_summary.tsv" \
  "S14 COG categories=${D10}/COG_summary.tsv" \
  "S15 DE summary=${D11}/*/DE_summary.tsv" \
  "S16 DEG sharing strains=${D11}/*/DEG_sharing.tsv" \
  "S17 Enrichment summary=${D11}/*/enrichment_summary.tsv" \
  "S18 Enrichment significant=${D11}/*/enrichment_significant.tsv" \
  "S19 Software=${D09}/software_versions.tsv"
echo
echo "Figures:"; ls -1 "${FIG}"/*.png 2>/dev/null | sed 's|^|  |'
