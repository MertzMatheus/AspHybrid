#!/usr/bin/env bash
# ============================================================
# 05_latus_subgenomes.sh                    (A. latus only, no argument)
#
# Assigns every region of the allodiploid A. latus genome to a parental
# subgenome: S = A. spinulosporus-derived, Q = A. quadrilineatus-like.
# Homeologs differ by ~7 % (Steenwyk et al. 2020, Mol Biol Evol), so 5-kb
# windows of A. latus aligned with minimap2 (asm20) to both parental genomes
# have a clearly higher identity to their parent of origin:
#   S if identity_S - identity_Q >= 2 %; Q if the reverse; else ambiguous.
# Every A. latus gene gets the label of its window in step 06 (subgenome=S|Q).
#
# Needs: step 01 for A_latus, A_spinulosporus and A_quadrilineatus.
# Output: final_results/05_latus_subgenomes/
#   A_latus.blocks.bed, A_latus.windows.tsv, A_latus.contigs.tsv,
#   A_latus_S.fasta / A_latus_Q.fasta + BUSCO of each subgenome,
#   subgenome_summary.tsv, figures/Fig_latus_subgenomes.{png,pdf}
# Time: ~20-40 min.
# ============================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"; source "$(dirname "$0")/lib/common.sh"
activate "${ENV_TOOLS}"
banner "05 A. latus subgenomes (S = ${LATUS_PARENT_S}, Q ~ ${LATUS_PARENT_Q})"

LAT="$(genome_fa A_latus)"; RS="$(genome_fa "${LATUS_PARENT_S}")"; RQ="$(genome_fa "${LATUS_PARENT_Q}")"
need "${LAT}" "${RS}" "${RQ}"
O="${D05}"; mkdir -p "${O}"

python3 "${BIN}/phase_subgenomes.py" split --fasta "${LAT}" --window 5000 --min-len 1000 --out "${O}/latus_windows.fa"
samtools faidx "${O}/latus_windows.fa"
for P in S Q; do
    REF="${RS}"; [[ "${P}" == Q ]] && REF="${RQ}"
    [[ -s "${O}/vs_${P}.paf" ]] || minimap2 -x asm20 -c --secondary=no -t "${THREADS}" "${REF}" "${O}/latus_windows.fa" \
        > "${O}/vs_${P}.paf" 2> "${O}/minimap2_${P}.log"
done
python3 "${BIN}/phase_subgenomes.py" assign --windows-fai "${O}/latus_windows.fa.fai" \
    --paf-s "${O}/vs_S.paf" --paf-q "${O}/vs_Q.paf" --delta 0.02 --out "${O}/A_latus"

busco_lineage
for P in S Q; do
    awk -v p="${P}" -F'\t' 'NR>1 && $7==p {print $1}' "${O}/A_latus.contigs.tsv" > "${O}/contigs_${P}.txt"
    samtools faidx "${LAT}" -r "${O}/contigs_${P}.txt" > "${O}/A_latus_${P}.fasta"
    compgen -G "${O}/busco_${P}/short_summary*.txt" >/dev/null || \
    busco -i "${O}/A_latus_${P}.fasta" -m genome -l "${BUSCO_LINEAGE_DIR}" --offline -c "${THREADS}" \
          -o "busco_${P}" --out_path "${O}" --force > "${O}/busco_${P}.log" 2>&1
done

python3 "${BIN}/subgenome_figure.py" --dir "${O}" --fig-dir "${FIG}"
mark_done "${O}"
echo
column -t -s $'\t' "${O}/subgenome_summary.tsv"
echo
echo "Contigs to inspect (S/Q switch inside the contig, or no parental hit):"
awk -F'\t' 'NR==1 || $9 ~ /CHECK/' "${O}/A_latus.contigs.tsv" | column -t -s $'\t' | head -20
echo "Figure: ${FIG}/Fig_latus_subgenomes.png   (expect two clear peaks, ~50 % S / 50 % Q)"
