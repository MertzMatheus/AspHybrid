#!/usr/bin/env bash
# ============================================================
# 02_repeat_masking.sh <species>            GOAL 2 (part 1): repeats
#
# De novo repeat library (RepeatModeler2 with LTR module) + RepeatMasker
# soft-masking (repeats in lower case), so that BRAKER3 does not call
# transposable elements genes.
# If the same genome was already masked in the earlier run (Results/Masking),
# that result is reused after checking that sequence names and lengths match.
#
# Output: final_results/02_repeat_masking/<sp>/<sp>.softmasked.fasta
#         <sp>-families.fa (library), rm_out/<sp>.fa.tbl (repeat classes),
#         <sp>.repeat_summary.tsv
# Time: 10-24 h per genome (or seconds when reused).
#   nohup bash 02_repeat_masking.sh A_latus > logs/02_A_latus.log 2>&1 &
# ============================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"; source "$(dirname "$0")/lib/common.sh"
species_arg "${RNASEQ_SPECIES[@]}" -- "$@"
activate "${ENV_TOOLS}"
banner "02 repeat masking: ${SP}"

G="$(genome_fa "${SP}")"; need "${G}" "${G}.fai"
W="${D02}/${SP}"; mkdir -p "${W}"
FINAL="$(masked_fa "${SP}")"
same_genome() { [[ -s "$1.fai" ]] && diff -q <(cut -f1,2 "$1.fai") <(cut -f1,2 "${G}.fai") >/dev/null; }

if is_done "${W}" && same_genome "${FINAL}"; then
    echo "Already done: ${FINAL}"
else
    PREV="${PREVIOUS_RESULTS}/Masking/${SP}"
    if [[ -s "${PREV}/${SP}.softmasked.fasta" ]] && { [[ -s "${PREV}/${SP}.softmasked.fasta.fai" ]] || samtools faidx "${PREV}/${SP}.softmasked.fasta"; } \
       && same_genome "${PREV}/${SP}.softmasked.fasta"; then
        echo "Reusing the masking of the earlier run: ${PREV}"
        cp "${PREV}/${SP}.softmasked.fasta" "${FINAL}"
        cp "${PREV}/${SP}-families.fa" "${W}/" 2>/dev/null || true
        [[ -d "${PREV}/rm_out" ]] && cp -r "${PREV}/rm_out" "${W}/"
    else
        cp -f "${G}" "${W}/${SP}.fa"
        docker run --rm --user "$(id -u):$(id -g)" -e HOME=/work -v "${W}:/work" -w /work "${IMG_TETOOLS}" bash -c "
            set -e
            if [[ ! -s ${SP}-families.fa ]]; then
                BuildDatabase -name ${SP} ${SP}.fa
                RepeatModeler -database ${SP} -threads ${THREADS} -LTRStruct
            fi
            RepeatMasker -pa $(( THREADS / 4 )) -xsmall -gff -lib ${SP}-families.fa -dir rm_out ${SP}.fa
        "
        if [[ -s "${W}/rm_out/${SP}.fa.masked" ]]; then cp "${W}/rm_out/${SP}.fa.masked" "${FINAL}"
        else echo "WARNING: no repeats found; unmasked genome used."; cp "${W}/${SP}.fa" "${FINAL}"; fi
        rm -rf "${W}"/RM_* "${W}/${SP}.fa"
    fi
    samtools faidx "${FINAL}"
    same_genome "${FINAL}" || die "masked genome does not match ${G}"
    mark_done "${W}"
fi

python3 - "${FINAL}" "${W}/${SP}-families.fa" "${SP}" "${W}/${SP}.repeat_summary.tsv" <<'PY'
import os, sys
f, fam, sp, out = sys.argv[1:]
low = tot = 0
for l in open(f):
    if l[0] != ">":
        s = l.strip(); tot += len(s); low += sum(c.islower() for c in s)
nf = sum(1 for l in open(fam) if l.startswith(">")) if os.path.exists(fam) else "NA"
open(out, "w").write(f"species\tgenome_bp\tmasked_bp\tmasked_percent\trepeat_families\n{sp}\t{tot}\t{low}\t{100*low/tot:.2f}\t{nf}\n")
print(f"{sp}: {100*low/tot:.2f}% soft-masked ({nf} de novo repeat families). Typical for Aspergillus: 3-10 %.")
PY
