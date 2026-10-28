#!/usr/bin/env bash
#SBATCH --job-name=rmRNA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=180:00:00
#SBATCH --mem=100G
#SBATCH --output=rmRNA_%j.out
#SBATCH --error=rmRNA_%j.err
# =============================================================================
# 02_remove_rRNA.sh - removal of ribosomal RNA reads (SortMeRNA)
#
#   databases: SILVA-based rRNA sequences distributed with SortMeRNA
#              bacterial 16S (id90), 23S (id98); eukaryotic 18S (id95), 28S (id98)
#   --paired_out: a pair is rRNA when either mate matches; only pairs in which
#                 neither mate matches are kept
#
# Input : ${BASE}/02_qc/<sample>_R1.trim.fastq.gz, _R2.trim.fastq.gz   (01_quality.sh)
# Output: ${BASE}/03_rmRNA/<sample>.mRNA_fwd.fq.gz, <sample>.mRNA_rev.fq.gz  (non-rRNA pairs)
#         ${BASE}/03_rmRNA/<sample>.rRNA_fwd/rev.fq.gz                      (rRNA pairs)
# Run   : BASE=/path/to/project bash 02_remove_rRNA.sh     (or sbatch)
# Samples with complete, valid output are skipped.
# =============================================================================
set -eo pipefail          # -u is set after the environment is activated

BASE="${BASE:-/path/to/project}"                 # project folder (edit or export BASE)
ENV_NAME="${ENV_NAME:-rnaseq_sa_ca}"             # environment with SortMeRNA (envs/preprocessing.yml)
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-25}}"
export OMP_NUM_THREADS="${THREADS}" MKL_NUM_THREADS="${THREADS}"

RAW="${BASE}/02_qc"
DB="${DB:-${BASE}/db/sortmerna}"                 # folder with the four rRNA FASTA files
OUT="${BASE}/03_rmRNA"

REFS=(
  "${DB}/silva-bac-16s-id90.fasta"
  "${DB}/silva-bac-23s-id98.fasta"
  "${DB}/silva-euk-18s-id95.fasta"
  "${DB}/silva-euk-28s-id98.fasta"
)

eval "$(micromamba shell hook --shell bash)"
micromamba activate "${ENV_NAME}"
set -u
mkdir -p "${OUT}"
for r in "${REFS[@]}"; do [[ -s "${r}" ]] || { echo "ERROR: rRNA database not found: ${r}"; exit 1; }; done
echo "Job: ${SLURM_JOB_ID:-local}   threads: ${THREADS}"

shopt -s nullglob
for R1 in "${RAW}"/*_R1.trim.fastq.gz; do
    SAMPLE=$(basename "${R1}" _R1.trim.fastq.gz)
    R2="${RAW}/${SAMPLE}_R2.trim.fastq.gz"
    WORKDIR="${OUT}/${SAMPLE}_work"
    OUT1="${OUT}/${SAMPLE}.mRNA_fwd.fq.gz"
    OUT2="${OUT}/${SAMPLE}.mRNA_rev.fq.gz"

    echo ">>> ${SAMPLE}"
    if [[ -f "${OUT1}" && -f "${OUT2}" ]]; then
        if gzip -t "${OUT1}" && gzip -t "${OUT2}"; then
            echo ">>> SKIP ${SAMPLE}: already done"
            continue
        fi
        echo ">>> damaged output, processing again"
        rm -f "${OUT1}" "${OUT2}"
    fi

    rm -rf "${WORKDIR}"; mkdir -p "${WORKDIR}"     # SortMeRNA refuses a non-empty work directory
    sortmerna \
        --ref "${REFS[0]}" --ref "${REFS[1]}" --ref "${REFS[2]}" --ref "${REFS[3]}" \
        --reads "${R1}" --reads "${R2}" \
        --paired_out --fastx --out2 \
        --threads "${THREADS}" \
        --workdir "${WORKDIR}" \
        --other "${OUT}/${SAMPLE}.mRNA" \
        --aligned "${OUT}/${SAMPLE}.rRNA"
    rm -rf "${WORKDIR}"
    echo ">>> finished ${SAMPLE}"
done
