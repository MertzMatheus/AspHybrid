#!/usr/bin/env bash
# ============================================================
# 07_orthology.sh                           GOAL 3: orthologues (no argument)
#
# "Gene X in A. fumigatus is gene Y in A. spinulosporus and gene(s) Z in A. latus."
#  - OrthoFinder (DIAMOND, gene-tree orthologues) on one protein per gene
#    (longest isoform, -T1) of three species: A_fumigatus, A_spinulosporus, A_latus.
#    A. quadrilineatus is NOT included.
#  - A. latus is one proteome. Being a hybrid, most A. fumigatus genes have
#    TWO A. latus orthologues (the two homeologs); the table lists all of them
#    (relation 1:2), with the subgenome of each one if step 05/06 labelled it.
#  - shared / species-specific genes (% per species, Venn of orthogroups)
#
# Needs: step 06 for A_fumigatus, A_spinulosporus and A_latus.
# Output: final_results/07_orthology/
#   Table_orthologs_Afumigatus.tsv    one row per A. fumigatus gene (AFUA_) -> ASPN gene(s), ALAT gene(s)
#   Table_pairs_Afumigatus_Alatus.tsv one row per AFUA-ALAT pair (easy to join to count matrices)
#   Table_pairs_Afumigatus_Aspinulosporus.tsv  same for ASPN
#   ortholog_statistics.tsv           % of A. fumigatus genes per relation type
#   tables/correspondence_<sp>.tsv    same, anchored on each species
#   tables/correspondence_by_orthogroup.tsv, tables/orthology_summary.tsv
#   sharing/  sharing_by_species.tsv, sharing_patterns.tsv, unique_genes_*.tsv
#   figures/Fig_orthologs.*, figures/Fig_gene_sharing.*
# Time: ~30-90 min.
# ============================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"; source "$(dirname "$0")/lib/common.sh"
activate "${ENV_TOOLS}"
banner "07 orthology: A. fumigatus vs A. spinulosporus vs A. latus"

ORTHO_SPECIES=(A_fumigatus A_spinulosporus A_latus)
O="${D07}"; PD="${O}/proteomes"; mkdir -p "${PD}" "${O}/tables" "${O}/sharing"
rm -f "${PD}"/*.faa "${PD}"/*.fa "${PD}"/*.fasta          # OrthoFinder uses EVERY file in this folder
for sp in "${ORTHO_SPECIES[@]}"; do
    f="$(ref_dir "${sp}")/${sp}.proteins.longest.faa"; need "${f}"
    seqkit replace -p '\s.*' -r '' "${f}" > "${PD}/${sp}.faa"
    printf "  %-18s %6d proteins\n" "${sp}" "$(grep -c '^>' "${PD}/${sp}.faa")"
done

# ---------- OrthoFinder: re-used only if run on exactly these proteins ----------
# (protein IDs compared with the run's SequenceIDs.txt; a changed proteome, e.g.
#  A. fumigatus AFUM -> official AFUA, triggers a new run; the old one is archived)
OFDIR="${O}/OrthoFinder"
OF="$(ls -d "${OFDIR}"/Results_* 2>/dev/null | sort | tail -1 || true)"
ids_now="$(cat "${PD}"/*.faa | grep '^>' | sed 's/^>//; s/\s.*//' | sort | md5sum | cut -c1-32)"
ids_run=""
[[ -n "${OF}" && -s "${OF}/WorkingDirectory/SequenceIDs.txt" ]] && \
    ids_run="$(sed 's/^[0-9_]*: //; s/\s.*//' "${OF}/WorkingDirectory/SequenceIDs.txt" | sort | md5sum | cut -c1-32)"
if [[ "${ids_now}" != "${ids_run}" ]]; then
    if [[ -d "${OFDIR}" ]]; then
        ARCH="${O}/previous_$(date +%Y%m%d_%H%M)"
        echo "Proteomes changed since the last OrthoFinder run - old results kept in ${ARCH}"
        mkdir -p "${ARCH}"
        mv "${OFDIR}" "${ARCH}/"
        for f in "${O}"/Table_*.tsv "${O}"/ortholog_statistics.tsv "${O}/tables" "${O}/sharing"; do
            [[ -e "${f}" ]] && mv "${f}" "${ARCH}/"
        done
        mkdir -p "${O}/tables" "${O}/sharing"
    fi
    activate "${ENV_ORTHO}"
    orthofinder -f "${PD}" -o "${OFDIR}" -S diamond -t "${THREADS}" -a "$(( THREADS / 4 ))" > "${O}/orthofinder.log" 2>&1 \
        || { tail -20 "${O}/orthofinder.log"; die "OrthoFinder failed (${O}/orthofinder.log)"; }
    activate "${ENV_TOOLS}"
    OF="$(ls -d "${OFDIR}"/Results_* | sort | tail -1)"
fi
echo "OrthoFinder results: ${OF}"
grep -E "Percentage of genes in orthogroups|Number of orthogroups with all species present|single-copy" \
     "${OF}/Comparative_Genomics_Statistics/Statistics_Overall.tsv" 2>/dev/null | head -5 || true

python3 "${BIN}/build_ortholog_table.py" --of "${OF}" --outdir "${O}/tables" \
        --id-map "$(ref_dir A_latus)/A_latus.id_map.tsv" > "${O}/tables.log"
python3 "${BIN}/ortholog_summary.py" --tables "${O}/tables" --out "${O}" --fig-dir "${FIG}" \
        --id-map "$(ref_dir A_latus)/A_latus.id_map.tsv"
python3 "${BIN}/gene_sharing.py" --of "${OF}" \
        --groups A_fumigatus=A_fumigatus A_spinulosporus=A_spinulosporus A_latus=A_latus \
        --short A_fumigatus=Af A_spinulosporus=As A_latus=Al --outdir "${O}/sharing"
for ext in png pdf; do cp -f "${O}/sharing/sharing_venn.${ext}" "${FIG}/Fig_gene_sharing.${ext}"; done
mark_done "${O}"
echo
column -t -s $'\t' "${O}/ortholog_statistics.tsv"
echo "Main table: ${O}/Table_orthologs_Afumigatus.tsv"
echo "Pairs     : ${O}/Table_pairs_Afumigatus_Alatus.tsv"
