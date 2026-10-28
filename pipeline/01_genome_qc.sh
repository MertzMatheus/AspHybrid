#!/usr/bin/env bash
# ============================================================
# 01_genome_qc.sh <species>                 GOAL 1: genome quality
#
#  - clean copy of the genome (IDs = first word of the header, upper case);
#    the original FASTA in Genomes/ is never modified
#  - assembly statistics: contigs, size, N50, L50, N90, largest, GC, gaps (N)
#  - BUSCO (genome mode, eurotiales_odb12.2): completeness and duplication.
#    Ploidy check: A. latus (allodiploid) should show ~2x genome size and a
#    high "Duplicated" fraction; haploid genomes mostly single-copy.
#  - after each species, the cross-species table and figure are updated
#
# Output: final_results/01_genome_qc/<sp>/  <sp>.genome.fasta, <sp>.genome_stats.tsv,
#                                           busco/, <sp>.header_map.tsv
#         final_results/01_genome_qc/genome_qc_summary.tsv
#         final_results/figures/Fig_genome_qc.{png,pdf}
# Run for all four species (A. quadrilineatus is needed for step 05).
# Time: ~10-20 min per species.
# ============================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"; source "$(dirname "$0")/lib/common.sh"
species_arg "${ALL_SPECIES[@]}" -- "$@"
activate "${ENV_TOOLS}"
banner "01 genome QC: ${SP} (${STRAIN[$SP]})"

O="${D01}/${SP}"; mkdir -p "${O}"
mapfile -t FA < <(find "${GENOMES}/${SP}" -maxdepth 1 -type f \( -iname '*.fasta' -o -iname '*.fa' -o -iname '*.fna' -o -iname '*.fasta.gz' -o -iname '*.fa.gz' -o -iname '*.fna.gz' \) | sort)
[[ ${#FA[@]} -eq 1 ]] || die "${GENOMES}/${SP} must contain exactly one FASTA (found ${#FA[@]})"
echo "Input: ${FA[0]}"

python3 "${BIN}/clean_genome.py" "${FA[0]}" "$(genome_fa "${SP}")" "${O}/${SP}.header_map.tsv" \
        "${O}/${SP}.genome_stats.tsv" "${SP}" "${STRAIN[$SP]}"
samtools faidx "$(genome_fa "${SP}")"

# ---------- BUSCO genome mode ----------
busco_lineage
STAMP="$(md5sum "$(genome_fa "${SP}")" | cut -c1-32)"
if [[ "$(cat "${O}/busco/.input_md5" 2>/dev/null)" != "${STAMP}" ]]; then
    rm -rf "${O}/busco"
    busco -i "$(genome_fa "${SP}")" -m genome -l "${BUSCO_LINEAGE_DIR}" --offline \
          -c "${THREADS}" -o busco --out_path "${O}" --force > "${O}/busco.log" 2>&1 \
        || { tail -20 "${O}/busco.log"; die "BUSCO failed (${O}/busco.log)"; }
    echo "${STAMP}" > "${O}/busco/.input_md5"
fi
echo "BUSCO genome: $(grep -h -o 'C:.*n:[0-9]*' "${O}"/busco/short_summary*.txt | head -1)"
mark_done "${O}"

# ---------- cross-species table + figure (species done so far) ----------
python3 "${BIN}/qc_figures.py" genome --qc-dir "${D01}" --fig-dir "${FIG}" --species "${ALL_SPECIES[@]}"
echo
column -t -s $'\t' "${D01}/genome_qc_summary.tsv"
echo "Figure: ${FIG}/Fig_genome_qc.png"
