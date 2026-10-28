#!/usr/bin/env bash
# ============================================================
# 10_eggnog.sh <species>                    functional annotation (eggNOG-mapper)
#
# eggNOG-mapper v2 (DIAMOND, local database in db/eggnog) on ONE protein per
# gene (<locus_tag>-T1, longest isoform) of A_fumigatus (official AFUA_ genes),
# A_spinulosporus and A_latus. All three species are annotated the same way,
# so functional categories and enrichment are comparable between species.
#
# Then one table per gene (locus_tag) with: product, eggNOG name/description,
# COG category, GO, EC, KEGG KO/pathway, CAZy, Pfam, orthologues (step 07),
# subgenome (A. latus), and the gene sets used for enrichment in step 11:
#   GO_BP / GO_MF / GO_CC  (eggNOG GO terms propagated to all ancestors,
#                           is_a + part_of, go-basic.obo)
#   KEGG_pathway, COG, KO, Pfam, CAZy
#
# Needs: step 06 (step 07 optional: orthologue columns)
# Output: final_results/10_functional_annotation/<sp>/
#   <sp>.emapper.annotations          raw eggNOG-mapper table
#   <sp>.gene_annotation.tsv          one row per gene  (Supplementary table)
#   <sp>.annotation_coverage.tsv      % genes with each annotation type
#   genesets/<sp>.<collection>.TERM2GENE.tsv / TERM2NAME.tsv
# + 10_functional_annotation/functional_annotation_summary.tsv and
#   figures/Fig_functional_annotation.{png,pdf} (all species done so far)
# Time: ~30-90 min per species with EGGNOG_DBMEM=1.
#   nohup bash 10_eggnog.sh A_latus > logs_10_A_latus.log 2>&1 &
# ============================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"; source "$(dirname "$0")/lib/common.sh"
# (also defined in lib/common.sh; kept here so an older common.sh still works)
D10="${D10:-${OUT}/10_functional_annotation}"; D11="${D11:-${OUT}/11_differential_expression}"
species_arg "${RNASEQ_SPECIES[@]}" -- "$@"
activate "${ENV_TOOLS}"
banner "10 eggNOG functional annotation: ${SP}"

R="$(ref_dir "${SP}")"; need "${R}/${SP}.proteins.longest.faa" "${R}/${SP}.proteins.faa"
[[ -s "${EGGNOG_DATA_DIR}/eggnog.db" && -s "${EGGNOG_DATA_DIR}/eggnog_proteins.dmnd" ]] \
    || die "eggNOG database incomplete in ${EGGNOG_DATA_DIR} (needs eggnog.db and eggnog_proteins.dmnd)"
E="${D10}/${SP}"; mkdir -p "${E}/genesets"

# ---------- reference vocabularies (downloaded once) ----------
if [[ ! -s "${GO_OBO}" ]]; then
    mkdir -p "$(dirname "${GO_OBO}")"
    echo "Downloading go-basic.obo (once) ..."
    wget -q -O "${GO_OBO}.tmp" https://purl.obolibrary.org/obo/go/go-basic.obo && mv "${GO_OBO}.tmp" "${GO_OBO}" \
        || die "could not download go-basic.obo; put it in ${GO_OBO}"
fi
if [[ ! -s "${KEGG_PATHWAY_NAMES}" ]]; then
    mkdir -p "$(dirname "${KEGG_PATHWAY_NAMES}")"
    echo "Downloading KEGG pathway names (once) ..."
    wget -q -O "${KEGG_PATHWAY_NAMES}.tmp" https://rest.kegg.jp/list/pathway && mv "${KEGG_PATHWAY_NAMES}.tmp" "${KEGG_PATHWAY_NAMES}" \
        || echo "WARNING: KEGG names not downloaded - pathways will be shown by ID only"
fi
echo "GO release: $(grep -m1 '^data-version' "${GO_OBO}" | cut -d' ' -f2)"

# ---------- eggNOG-mapper (re-run only if the proteins changed) ----------
STAMP="$(md5sum "${R}/${SP}.proteins.longest.faa" | cut -c1-32)"
if [[ "$(cat "${E}/.input_md5" 2>/dev/null)" != "${STAMP}" || ! -s "${E}/${SP}.emapper.annotations" ]]; then
    MEM=(); [[ "${EGGNOG_DBMEM}" == 1 ]] && MEM=(--dbmem)
    echo "eggNOG-mapper on $(grep -c '^>' "${R}/${SP}.proteins.longest.faa") proteins (one per gene) ..."
    mkdir -p "${TMPDIR}/emapper_${SP}"
    emapper.py -i "${R}/${SP}.proteins.longest.faa" --itype proteins -m diamond \
        --data_dir "${EGGNOG_DATA_DIR}" --tax_scope "${EGGNOG_TAX_SCOPE}" --go_evidence "${EGGNOG_GO_EVIDENCE}" \
        --cpu "${THREADS}" "${MEM[@]}" --temp_dir "${TMPDIR}/emapper_${SP}" \
        --output_dir "${E}" -o "${SP}" --override > "${E}/emapper.log" 2>&1 \
        || { tail -30 "${E}/emapper.log"; die "eggNOG-mapper failed (log: ${E}/emapper.log)"; }
    rm -rf "${TMPDIR}/emapper_${SP}"
    echo "${STAMP}" > "${E}/.input_md5"
fi
emapper.py --version 2>&1 | grep -m1 -o 'emapper-[0-9.]*' > "${E}/emapper_version.txt" || true

# ---------- gene table + gene sets ----------
CORR="${D07}/tables/correspondence_${SP}.tsv"
python3 "${BIN}/functional_db.py" --species "${SP}" --proteins "${R}/${SP}.proteins.faa" \
    --eggnog "${E}/${SP}.emapper.annotations" --obo "${GO_OBO}" --kegg-names "${KEGG_PATHWAY_NAMES}" \
    --id-map "${R}/${SP}.id_map.tsv" $( [[ -s "${CORR}" ]] && echo "--correspondence ${CORR}" ) \
    --out "${E}/${SP}" --geneset-dir "${E}/genesets"
mark_done "${E}"

python3 "${BIN}/functional_figure.py" --dir "${D10}" --fig-dir "${FIG}" --species "${RNASEQ_SPECIES[@]}"
echo "Gene table: ${E}/${SP}.gene_annotation.tsv      Figure: ${FIG}/Fig_functional_annotation.png"
