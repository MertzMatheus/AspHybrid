# ============================================================
# lib/common.sh - shared functions and standard paths (sourced by every step)
# ============================================================
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${HERE}/bin"

die()  { echo >&2; echo "ERROR: $*" >&2; exit 1; }
need() { local f; for f in "$@"; do [[ -s "${f}" ]] || die "missing input: ${f}
       -> run the earlier step that produces it"; done; }
banner() {
    echo "=================================================================="
    echo " $*"
    echo " $(date '+%F %T')   $(hostname)   threads: ${THREADS}"
    echo "=================================================================="
}
activate() { eval "$(micromamba shell hook --shell bash)"; micromamba activate "$1"; }

in_list() { local x="$1"; shift; [[ " $* " == *" ${x} "* ]]; }
species_arg() {     # species_arg <allowed list...> -- "$@"
    local allowed=()
    while [[ "$1" != "--" ]]; do allowed+=("$1"); shift; done; shift
    SP="${1:-}"
    [[ -n "${SP}" ]] || die "usage: bash $(basename "$0") <species>   (${allowed[*]})"
    in_list "${SP}" "${allowed[@]}" || die "'${SP}' not valid here. Use one of: ${allowed[*]}"
}
done_flag() { echo "$1/.done"; }
mark_done() { date '+%F %T' > "$(done_flag "$1")"; }
is_done()   { [[ -f "$(done_flag "$1")" ]]; }

# ---------------- standard output paths ----------------
D01="${OUT}/01_genome_qc"
D02="${OUT}/02_repeat_masking"
D03="${OUT}/03_rnaseq_mapping"
D04="${OUT}/04_gene_prediction"
D05="${OUT}/05_latus_subgenomes"
D06="${OUT}/06_reference_files"
D07="${OUT}/07_orthology"
D08="${OUT}/08_count_matrices"
D09="${OUT}/09_tables"
D10="${OUT}/10_functional_annotation"
D11="${OUT}/11_differential_expression"
FIG="${OUT}/figures"
LOGS="${OUT}/logs"
mkdir -p "${OUT}" "${FIG}" "${LOGS}"

genome_fa() { echo "${D01}/$1/$1.genome.fasta"; }          # clean genome (all later steps)
masked_fa() { echo "${D02}/$1/$1.softmasked.fasta"; }
ref_dir()   { echo "${D06}/$1"; }

# BUSCO lineage: unpack the .tar.gz once if needed
busco_lineage() {
    if [[ ! -s "${BUSCO_LINEAGE_DIR}/dataset.cfg" ]]; then
        need "${BUSCO_LINEAGE_TGZ}"
        echo "Unpacking $(basename "${BUSCO_LINEAGE_TGZ}") (once, a few minutes) ..."
        tar -xzf "${BUSCO_LINEAGE_TGZ}" -C "$(dirname "${BUSCO_LINEAGE_TGZ}")"
    fi
    need "${BUSCO_LINEAGE_DIR}/dataset.cfg"
}

# ---------------- RNA-seq sample sheet ----------------
# columns (header line starting with '#'): species sample R1 R2 strain condition replicate
# check_sheet <sp> : stops with a clear message if the sheet or the species' FASTQ files are missing
check_sheet() {
    [[ -s "${SAMPLES}" ]] || die "sample sheet not found: ${SAMPLES}
       put rnaseq_samples.tsv next to config.sh, or set SAMPLES= in config.sh"
    local n; n="$(sheet_samples "$1" | wc -l)"
    (( n > 0 )) || die "no $1 rows in ${SAMPLES} (first column must be exactly '$1', tab-separated)"
    local bad=0 s r1 r2 grp rep
    while IFS=$'\t' read -r s r1 r2 grp rep; do
        for f in "${r1}" "${r2}"; do
            [[ -s "${f}" ]] || { echo "  missing FASTQ for ${s}: ${f}" >&2; bad=1; }
        done
    done < <(sheet_samples "$1")
    (( bad == 0 )) || die "FASTQ files listed in ${SAMPLES} not found (see above)"
}
# sheet_samples <sp> : prints "sample<TAB>R1<TAB>R2<TAB>group<TAB>replicate", group = <strain>_<condition>
sheet_samples() {
    awk -F'\t' -v sp="$1" '
        /^#/ { if (!h) { for (i = 1; i <= NF; i++) { c = $i; gsub(/^[# ]+|[ ]+$/, "", c); col[c] = i }
                         if ("sample" in col) h = 1 }
               next }
        !h   { col["species"]=1; col["sample"]=2; col["R1"]=3; col["R2"]=4; col["condition"]=5; col["replicate"]=6; h = 1 }
        $col["species"] == sp {
            g = $col["condition"]
            if ("strain" in col && $col["strain"] != "") g = $col["strain"] "_" g
            print $col["sample"] "\t" $col["R1"] "\t" $col["R2"] "\t" g "\t" $col["replicate"] }' "${SAMPLES}"
}
