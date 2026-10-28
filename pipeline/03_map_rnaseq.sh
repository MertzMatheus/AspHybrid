#!/usr/bin/env bash
# ============================================================
# 03_map_rnaseq.sh <species>                GOAL 4: RNA-seq mapping
#
# Maps every RNA-seq sample of the species to its OWN genome.
# Reads are already clean (adapters and rRNA removed), so they are NOT
# trimmed again; fastp runs in report-only mode for the QC table.
#   HISAT2 (spliced; --dta; max intron 5 kb, fungal introns are short)
#   -> coordinate-sorted, indexed BAM per sample + samtools flagstat
# A. latus reads are mapped to the full hybrid genome (S + Q subgenomes),
# never to one parent only, so each read can find its own homeolog.
# These BAMs are used twice: as RNA-seq evidence for BRAKER3 (step 04)
# and for counting reads per gene (step 08).
#
# Output: final_results/03_rnaseq_mapping/<sp>/bam/<sample>.bam (+.bai)
#         .../qc/<sample>.fastp.json|html, <sample>.hisat2.summary, <sample>.flagstat
#         .../<sp>.mapping_rates.tsv
# Time: ~10-20 min per sample.
#   nohup bash 03_map_rnaseq.sh A_latus > logs/03_A_latus.log 2>&1 &
# ============================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"; source "$(dirname "$0")/lib/common.sh"
species_arg "${RNASEQ_SPECIES[@]}" -- "$@"
activate "${ENV_TOOLS}"
banner "03 RNA-seq mapping: ${SP}"

G="$(genome_fa "${SP}")"; need "${G}"
M="${D03}/${SP}"; mkdir -p "${M}/bam" "${M}/qc" "${M}/index"
# only one run per species at a time (two runs would overwrite each other's files)
exec 9> "${M}/.lock"
flock -n 9 || die "03_map_rnaseq.sh ${SP} is already running (lock: ${M}/.lock). Check with: ps aux | grep 03_map_rnaseq"

# ---------- check the sample sheet ----------
check_sheet "${SP}"
mapfile -t ROWS < <(sheet_samples "${SP}")
echo "${#ROWS[@]} paired-end samples"

# ---------- index ----------
if [[ ! -s "${M}/index/genome.1.ht2" ]]; then
    echo "Building HISAT2 index ..."
    hisat2-build -p "${THREADS}" "${G}" "${M}/index/genome" > "${M}/index/build.log" 2>&1
fi

# ---------- per sample ----------
for row in "${ROWS[@]}"; do
    IFS=$'\t' read -r s r1 r2 grp rep <<< "${row}"
    BAM="${M}/bam/${s}.bam"
    if [[ -s "${BAM}.bai" && -s "${M}/qc/${s}.hisat2.summary" ]]; then echo "  ${s}: done"; continue; fi
    echo "  ${s} (${grp}, rep ${rep})"
    mkdir -p "${M}/bam" "${M}/qc"
    rm -f "${BAM}.tmp.bam" "${TMPDIR}/${SP}_${s}".*.bam        # leftovers of an interrupted run
    fastp -i "${r1}" -I "${r2}" -w 8 \
          --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering --disable_trim_poly_g \
          -j "${M}/qc/${s}.fastp.json" -h "${M}/qc/${s}.fastp.html" 2> "${M}/qc/${s}.fastp.log"
    hisat2 -p "${THREADS}" --dta --max-intronlen 5000 -x "${M}/index/genome" -1 "${r1}" -2 "${r2}" \
           --rg-id "${s}" --rg "SM:${s}" --summary-file "${M}/qc/${s}.hisat2.summary" 2> "${M}/qc/${s}.hisat2.log" \
      | samtools sort -@ 8 -m 1500M -T "${TMPDIR}/${SP}_${s}" -o "${BAM}.tmp.bam" - \
      || die "${s}: mapping/sorting failed. Does ${M}/bam still exist? (ls -ld ${M}/bam; df -h ${M})"
    mv "${BAM}.tmp.bam" "${BAM}"
    samtools index "${BAM}"
    samtools quickcheck "${BAM}" || die "${s}: BAM is damaged"
    samtools flagstat -@ 4 "${BAM}" > "${M}/qc/${s}.flagstat"
    echo "    $(grep overall "${M}/qc/${s}.hisat2.summary")"
done

# ---------- mapping-rate table ----------
python3 - "${M}" "${SP}" "${ROWS[@]}" <<'PY'
import json, re, sys
M, sp, rows = sys.argv[1], sys.argv[2], sys.argv[3:]
out = open(f"{M}/{sp}.mapping_rates.tsv", "w")
out.write("species\tsample\tgroup\treplicate\tread_pairs\tQ30_percent\tGC_percent\tdup_percent\t"
          "aligned_pairs_unique\taligned_pairs_multi\toverall_alignment_percent\n")
for row in rows:
    s, r1, r2, grp, rep = row.split("\t")
    j = json.load(open(f"{M}/qc/{s}.fastp.json"))["summary"]["before_filtering"]
    dup = json.load(open(f"{M}/qc/{s}.fastp.json")).get("duplication", {}).get("rate", 0)
    t = open(f"{M}/qc/{s}.hisat2.summary").read()
    g = lambda p: (re.search(p, t) or [None, "NA"])[1]
    uniq = g(r"(\d+) \S+ aligned concordantly exactly 1 time")
    multi = g(r"(\d+) \S+ aligned concordantly >1 times")
    rate = g(r"([\d.]+)% overall alignment rate")
    vals = [sp, s, grp, rep, j["total_reads"] // 2, f"{100 * j['q30_rate']:.2f}", f"{100 * j['gc_content']:.2f}",
            f"{100 * dup:.2f}", uniq, multi, rate]
    out.write("\t".join(map(str, vals)) + "\n")
out.close()
PY
mark_done "${M}"
echo
column -t -s $'\t' "${M}/${SP}.mapping_rates.tsv"
echo "Expected: overall alignment > 85-90 % (reads of the right species on its own genome)."
