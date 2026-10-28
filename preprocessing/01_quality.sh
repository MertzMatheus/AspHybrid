#!/usr/bin/env bash
#SBATCH --job-name=fastp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=180:00:00
#SBATCH --mem=100G
#SBATCH --output=fastp_%j.out
#SBATCH --error=fastp_%j.err
# =============================================================================
# 01_quality.sh - adapter and quality trimming of paired-end reads (fastp)
#
#   adapters detected automatically for paired-end data
#   base qualified at Phred >= 20 (reads with > 40 % unqualified bases removed)
#   reads shorter than 30 nt after trimming removed (with their mate)
#
# Input : ${BASE}/01.RawData/<sample>/<sample>_1.fq.gz and _2.fq.gz
# Output: ${BASE}/02_qc/<sample>_R1.trim.fastq.gz, <sample>_R2.trim.fastq.gz
#         ${BASE}/02_qc/reports/<sample>.html / .json
# Run   : BASE=/path/to/project bash 01_quality.sh     (or sbatch on a SLURM cluster)
# Samples already processed are skipped.
# =============================================================================
set -eo pipefail          # -u is set after the environment is activated

BASE="${BASE:-/path/to/project}"          # project folder (edit or export BASE)
ENV_NAME="${ENV_NAME:-rnaseq_sa_ca}"      # micromamba environment with fastp (envs/preprocessing.yml)
FASTP_THREADS="${FASTP_THREADS:-4}"

raw_data="${BASE}/01.RawData"
raw_trimmed="${BASE}/02_qc"
fastp_reports="${raw_trimmed}/reports"

eval "$(micromamba shell hook --shell bash)"
micromamba activate "${ENV_NAME}"
set -u
mkdir -p "${fastp_reports}" "${raw_trimmed}"

shopt -s nullglob
for R1_file in "${raw_data}"/*/*_1.fq.gz; do
    dir=$(dirname "${R1_file}")
    base=$(basename "${R1_file}" _1.fq.gz)
    R2_file="${dir}/${base}_2.fq.gz"
    [[ -f "${R2_file}" ]] || { echo "ERROR: missing ${R2_file}"; exit 1; }

    out_R1="${raw_trimmed}/${base}_R1.trim.fastq.gz"
    out_R2="${raw_trimmed}/${base}_R2.trim.fastq.gz"
    out_html="${fastp_reports}/${base}.html"
    out_json="${fastp_reports}/${base}.json"
    out_unp1="${raw_trimmed}/${base}_unpaired_1.fq.gz"
    out_unp2="${raw_trimmed}/${base}_unpaired_2.fq.gz"

    if [[ -f "${out_R1}" && -f "${out_R2}" && -f "${out_html}" ]]; then
        echo "SKIP ${base}: output already exists"
        continue
    fi

    echo "Processing ${base}"
    fastp \
        -i "${R1_file}" -I "${R2_file}" \
        -o "${out_R1}" -O "${out_R2}" \
        --unpaired1 "${out_unp1}" --unpaired2 "${out_unp2}" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 30 \
        --thread "${FASTP_THREADS}" \
        --html "${out_html}" --json "${out_json}" \
        --report_title "fastp report for ${base}"
    echo "Done ${base}"
done

echo "All samples processed."
