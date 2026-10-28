#!/usr/bin/env bash
# ============================================================
# 04_gene_prediction.sh <species>           GOAL 2 (part 2): gene models
#
# BRAKER3 in ETP mode (Docker image pinned in config.sh):
#   genome      : soft-masked genome (step 02)
#   RNA-seq     : BAMs from step 03; one replicate per group is merged
#                 (BRAKER_EVIDENCE_REPLICATE) - every condition is
#                 represented, the file stays small
#   proteins    : OrthoDB v12 Fungi partition (db/orthodb)
#   --fungus    : GeneMark fungal branch-point model
# The trained AUGUSTUS model is kept in 04_gene_prediction/augustus_config.
#
# Output: final_results/04_gene_prediction/<sp>/braker/braker.gff3 (raw models)
# Time: ~6-15 h per genome.
#   nohup bash 04_gene_prediction.sh A_latus > logs/04_A_latus.log 2>&1 &
# ============================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"; source "$(dirname "$0")/lib/common.sh"
species_arg "${RNASEQ_SPECIES[@]}" -- "$@"
activate "${ENV_TOOLS}"
banner "04 gene prediction (BRAKER3 ETP): ${SP}"
[[ -s "${SAMPLES}" ]] || die "sample sheet not found: ${SAMPLES}"

G="$(masked_fa "${SP}")"; need "${G}"
P="${D04}/${SP}"; WD="${P}/braker"; mkdir -p "${WD}"
if is_done "${P}" && [[ -s "${WD}/braker.gff3" ]]; then echo "Already done: ${WD}/braker.gff3"; exit 0; fi

# ---------- protein evidence (built once, shared by all species) ----------
need "${ORTHODB_FUNGI}"
PROT="${D04}/OrthoDB12_Fungi.proteins.faa"
if [[ ! -s "${PROT}" ]]; then
    echo "Preparing protein evidence (headers = ID only, no duplicates) ..."
    seqkit seq -i -w 0 "${ORTHODB_FUNGI}" | seqkit rmdup -s -j 8 > "${PROT}.tmp" && mv "${PROT}.tmp" "${PROT}"
fi
echo "Protein evidence: $(grep -c '^>' "${PROT}") sequences"

# ---------- RNA-seq evidence: one replicate per group ----------
BAMS=()
while IFS=$'\t' read -r s r1 r2 grp rep; do
    [[ "${rep}" == "${BRAKER_EVIDENCE_REPLICATE}" ]] || continue
    need "${D03}/${SP}/bam/${s}.bam"
    BAMS+=("${D03}/${SP}/bam/${s}.bam")
done < <(sheet_samples "${SP}")
(( ${#BAMS[@]} > 0 )) || die "no BAMs of replicate ${BRAKER_EVIDENCE_REPLICATE} - run 03_map_rnaseq.sh ${SP} first"
EV="${P}/${SP}.rnaseq_evidence.bam"
# Only mapped, primary alignments; each sample capped at EVIDENCE_MAX_PAIRS read pairs
# (random, pair-consistent subsampling, fixed seed). Gene-structure evidence (introns,
# transcript ends) saturates far below this depth, and a smaller BAM avoids the
# BRAKER-internal 'samtools sort' failure on very large inputs.
CAP="${EVIDENCE_MAX_PAIRS:-0}"
STAMP="${EV}.params"
if [[ ! -s "${EV}.bai" || "$(cat "${STAMP}" 2>/dev/null)" != "cap=${CAP}" ]]; then
    echo "Building RNA-seq evidence from ${#BAMS[@]} BAMs (one per group), cap ${CAP} pairs/sample (0 = no cap)"
    rm -f "${EV}" "${EV}.bai" "${STAMP}"
    T="${P}/evidence_tmp"; rm -rf "${T}"; mkdir -p "${T}"
    printf "sample\tmapped_pairs\tfraction_kept\n" > "${P}/evidence_subsampling.tsv"
    for b in "${BAMS[@]}"; do
        s="$(basename "${b}" .bam)"
        n="$(samtools view -c -@ 8 -f 0x40 -F 0x904 "${b}")"        # mapped primary first mates = pairs
        sub=()
        if (( CAP > 0 && n > CAP )); then
            frac="$(awk -v c="${CAP}" -v n="${n}" 'BEGIN{printf "%.6f", c/n}')"
            sub=(-s "$(awk -v f="${frac}" 'BEGIN{printf "42%s", substr(f, 2)}')")   # seed 42 + fraction, e.g. 42.531
        else
            frac=1
        fi
        printf "%s\t%s\t%s\n" "${s}" "${n}" "${frac}" >> "${P}/evidence_subsampling.tsv"
        echo "  ${s}: ${n} mapped pairs, keep ${frac}"
        samtools view -b -@ 8 -F 0x904 "${sub[@]}" -o "${T}/${s}.bam" "${b}"
    done
    samtools merge -f -@ 8 -o "${EV}" "${T}"/*.bam
    samtools index "${EV}"
    rm -rf "${T}"
    echo "cap=${CAP}" > "${STAMP}"
fi
echo "Evidence BAM: $(du -h "${EV}" | cut -f1)"
samtools quickcheck "${EV}" || die "evidence BAM is damaged: ${EV}"

# ---------- BRAKER3 ----------
cpath() { echo "/data${1#${BASE}}"; }
AUG="${D04}/augustus_config"
if [[ ! -d "${AUG}/species" ]]; then
    mkdir -p "${AUG}"
    docker run --rm --user "$(id -u):$(id -g)" -v "${BASE}:/data" "${IMG_BRAKER}" \
        bash -c "cp -r \"\${AUGUSTUS_CONFIG_PATH:-/opt/Augustus/config}/.\" \"$(cpath "${AUG}")/\""
fi
MODEL="${SP}_final"
rm -rf "${AUG}/species/${MODEL}" "${WD:?}"/*                 # clean restart
# ---------- spaln watchdog ----------
# ProtHint aligns each protein to its genome region with spaln (seconds per pair). A few
# pairs (very long proteins, e.g. PKS/NRPS) can make spaln run for days on one core and
# block BRAKER. Any single spaln process older than SPALN_MAX_HOURS is killed: only the
# hints of that one protein are lost; ProtHint and BRAKER continue normally.
CNAME="braker_${SP}"
docker rm -f "${CNAME}" > /dev/null 2>&1 || true
WDOG="${P}/spaln_watchdog.log"
(
    lim=$(( ${SPALN_MAX_HOURS:-2} * 3600 ))
    sleep 120
    while docker inspect "${CNAME}" > /dev/null 2>&1; do
        docker top "${CNAME}" -eo pid,etimes,comm 2>/dev/null | awk -v l="${lim}" '$3=="spaln" && $2>l {print $1, $2}' |
        while read -r pid el; do
            echo "$(date '+%F %T')  killed spaln PID ${pid} after $(( el / 3600 )) h (one protein, hints skipped)" >> "${WDOG}"
            kill "${pid}" 2>/dev/null || true
        done
        sleep 600
    done
) &
WATCHDOG_PID=$!
trap 'kill ${WATCHDOG_PID} 2>/dev/null || true' EXIT

echo "BRAKER3 running; log: ${WD}/braker.log   (container: ${CNAME})"
if ! docker run --rm --name "${CNAME}" --user "$(id -u):$(id -g)" -e HOME=/tmp -e TMPDIR="$(cpath "${WD}")" \
        -e AUGUSTUS_CONFIG_PATH="$(cpath "${AUG}")" -v "${BASE}:/data" "${IMG_BRAKER}" \
        braker.pl --genome="$(cpath "${G}")" --prot_seq="$(cpath "${PROT}")" --bam="$(cpath "${EV}")" \
                  --fungus --softmasking --gff3 --species="${MODEL}" \
                  --AUGUSTUS_CONFIG_PATH="$(cpath "${AUG}")" \
                  --threads="${THREADS}" --workingdir="$(cpath "${WD}")"; then
    echo; echo "---- BRAKER3 error files (non-empty) ----"
    for f in "${WD}"/errors/*; do [[ -s "${f}" ]] && { echo "== ${f}"; tail -15 "${f}"; }; done
    die "BRAKER3 failed - send the lines above"
fi
need "${WD}/braker.gff3"
[[ -s "${WDOG}" ]] && { echo "spaln watchdog:"; cat "${WDOG}"; }
mark_done "${P}"
echo "Raw models: $(awk -F'\t' '$3=="gene"' "${WD}/braker.gff3" | wc -l) genes, $(awk -F'\t' '$3=="mRNA"' "${WD}/braker.gff3" | wc -l) transcripts"
