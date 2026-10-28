#!/usr/bin/env bash
# ============================================================
# 06_reference_files.sh <species>           GOAL 2 (part 3): .gff3 .gtf .faa + locus_tags
#
# Curates the BRAKER3 models and writes the reference files used for
# orthology (step 07) and RNA-seq quantification (step 08):
#  - removed: internal stop codon, protein < 50 aa, CDS > 50 % repeat-masked,
#    duplicated isoforms (identical CDS); partial genes kept and flagged
#  - locus_tags in genome order, step 10:  AFUM_000010, AFUM_000020 ...
#    (ASPN_ for A. spinulosporus, ALAT_ for A. latus); transcripts and
#    proteins <locus_tag>-T1, -T2 (T1 = longest protein)
#  - A. latus: subgenome=S|Q|ambiguous on every gene (needs step 05)
#  - checks: gffread re-translation of the GTF must reproduce >= 98 % of the
#    proteins; BUSCO on one protein per gene
#
# Output: final_results/06_reference_files/<sp>/
#   <sp>.genome.fasta  <sp>.gff3  <sp>.gtf  <sp>.proteins.faa  <sp>.proteins.longest.faa
#   <sp>.cds.fna  <sp>.transcripts.fna  <sp>.id_map.tsv  <sp>.annotation_summary.tsv
#   <sp>.removed_models.tsv  validation.tsv  busco_proteins/  md5sums.txt
# + annotation_summary_all.tsv and figures/Fig_annotation_qc.{png,pdf}
# Time: ~10-20 min. Safe to re-run (locus_tags are deterministic).
# ============================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"; source "$(dirname "$0")/lib/common.sh"
species_arg "${RNASEQ_SPECIES[@]}" -- "$@"
activate "${ENV_TOOLS}"
OFFICIAL=0; [[ "${SP}" == A_fumigatus && "${AF_ANNOTATION:-braker}" == official ]] && OFFICIAL=1
LT_PREFIX="${PREFIX[$SP]}"; [[ "${SP}" == A_fumigatus && ${OFFICIAL} == 0 ]] && LT_PREFIX=AFUM
banner "06 reference files: ${SP} (locus_tags ${LT_PREFIX}_...$( ((OFFICIAL)) && echo ', official Af293 RefSeq annotation'))"

need "$(genome_fa "${SP}")"
R="$(ref_dir "${SP}")"; mkdir -p "${R}"

if (( OFFICIAL )); then
    # ---------- A. fumigatus: official AFUA_ annotation (NCBI RefSeq Af293) ----------
    # keep the earlier de novo AFUM annotation (and everything made with it)
    OLD="${D06}/A_fumigatus_BRAKER_AFUM"
    if [[ -s "${R}/${SP}.gff3" ]] && grep -q -m1 "ID=AFUM_" "${R}/${SP}.gff3"; then
        echo "Keeping the de novo AFUM annotation in ${OLD}"
        rm -rf "${OLD}"; mv "${R}" "${OLD}"; mkdir -p "${R}"
    fi
    mkdir -p "${AF293_DIR}"
    for ext in genomic.gff.gz genomic.fna.gz; do
        f="${AF293_DIR}/${AF293_ACC}_${ext}"
        if [[ ! -s "${f}" ]] || ! gzip -t "${f}" 2>/dev/null; then
            echo "Downloading ${AF293_ACC}_${ext} (once) ..."
            wget -q -c -O "${f}" "${AF293_URL}/${AF293_ACC}_${ext}" || die "download failed: ${AF293_URL}/${AF293_ACC}_${ext}"
        fi
    done
    python3 "${BIN}/refseq_to_project.py" --gff "${AF293_DIR}/${AF293_ACC}_genomic.gff.gz" \
        --refseq-fasta "${AF293_DIR}/${AF293_ACC}_genomic.fna.gz" --genome "$(genome_fa "${SP}")" \
        --out "${R}/${SP}.curated" > "${R}/curation.log" 2>&1 || { cat "${R}/curation.log"; die "official annotation could not be placed"; }
else
    RAW="${D04}/${SP}/braker/braker.gff3"
    need "${RAW}" "$(masked_fa "${SP}")"
    SUB=()
    if [[ "${SP}" == A_latus ]]; then
        need "${D05}/A_latus.blocks.bed"          # run 05_latus_subgenomes.sh first
        SUB=(--subgenome-bed "${D05}/A_latus.blocks.bed")
    fi
    python3 "${BIN}/curate_annotation.py" --gff "${RAW}" --fasta "$(masked_fa "${SP}")" \
        --prefix "${LT_PREFIX}" --source BRAKER3 --min-aa 50 --max-masked 0.5 "${SUB[@]}" \
        --out "${R}/${SP}.curated" > "${R}/curation.log"
fi
cat "${R}/curation.log"
cp -f "$(genome_fa "${SP}")"            "${R}/${SP}.genome.fasta"
mv -f "${R}/${SP}.curated.gff3"           "${R}/${SP}.gff3"
mv -f "${R}/${SP}.curated.proteins.faa"   "${R}/${SP}.proteins.faa"
mv -f "${R}/${SP}.curated.proteins.longest.faa" "${R}/${SP}.proteins.longest.faa"
mv -f "${R}/${SP}.curated.cds.fna"        "${R}/${SP}.cds.fna"
mv -f "${R}/${SP}.curated.id_map.tsv"     "${R}/${SP}.id_map.tsv"
mv -f "${R}/${SP}.curated.removed.tsv"    "${R}/${SP}.removed_models.tsv"
mv -f "${R}/${SP}.curated.stats.tsv"      "${R}/${SP}.curation_stats.tsv"
[[ -f "${R}/${SP}.curated.seqid_map.tsv" ]] && mv -f "${R}/${SP}.curated.seqid_map.tsv" "${R}/${SP}.refseq_seqid_map.tsv"

python3 "${BIN}/gff3_to_gtf.py" "${R}/${SP}.gff3" "${R}/${SP}.gtf"
samtools faidx "${R}/${SP}.genome.fasta"
gffread "${R}/${SP}.gtf" -g "${R}/${SP}.genome.fasta" -w "${R}/${SP}.transcripts.fna" -y "${R}/check.faa"
python3 "${BIN}/validate_reference.py" --dir "${R}" --species "${SP}"
rm -f "${R}/check.faa"
grep -q $'\tOK\t' "${R}/validation.tsv" || die "validation FAILED - see ${R}/validation.tsv"

busco_lineage
STAMP="$(md5sum "${R}/${SP}.proteins.longest.faa" | cut -c1-32)"
if [[ "$(cat "${R}/busco_proteins/.input_md5" 2>/dev/null)" != "${STAMP}" ]]; then
    rm -rf "${R}/busco_proteins"
    busco -i "${R}/${SP}.proteins.longest.faa" -m proteins -l "${BUSCO_LINEAGE_DIR}" --offline \
          -c "${THREADS}" -o busco_proteins --out_path "${R}" --force > "${R}/busco.log" 2>&1 \
        || { tail -20 "${R}/busco.log"; die "BUSCO failed"; }
    echo "${STAMP}" > "${R}/busco_proteins/.input_md5"
fi

# ---------- A. fumigatus: link the earlier AFUM genes to the official AFUA genes ----------
if (( OFFICIAL )) && [[ -s "${D06}/A_fumigatus_BRAKER_AFUM/A_fumigatus.gff3" ]]; then
    echo "AFUM (BRAKER3) -> AFUA (official) correspondence:"
    TR=(); for f in "${D07}/Table_orthologs_Afumigatus.tsv" "${D07}/Table_pairs_Afumigatus_Alatus.tsv" \
                   "${D07}/Table_pairs_Afumigatus_Aspinulosporus.tsv"; do
              [[ -s "${f}" ]] && head -2 "${f}" | grep -q AFUM_ && TR+=("${f}"); done
    python3 "${BIN}/map_afum_afua.py" --afum "${D06}/A_fumigatus_BRAKER_AFUM/A_fumigatus.gff3" \
        --afua "${R}/${SP}.gff3" --out "${R}" --translate "${TR[@]}"
fi

# ---------- per-species summary ----------
python3 - "${R}" "${SP}" "${STRAIN[$SP]}" "${PREFIX[$SP]}" <<'PY'
import glob, re, statistics as st, sys
R, sp, strain, prefix = sys.argv[1:]
hdr, seq, name = {}, {}, None
for l in open(f"{R}/{sp}.proteins.faa"):
    l = l.rstrip("\n")
    if l.startswith(">"): name = l[1:].split()[0]; hdr[name] = l; seq[name] = ""
    else: seq[name] += l
t1 = [t for t in seq if t.endswith("-T1")]
cs = dict(l.rstrip("\n").split("\t") for l in open(f"{R}/{sp}.curation_stats.tsv"))
sub = {}
for t in t1:
    m = re.search(r"subgenome=(\S+)", hdr[t])
    if m: sub[m.group(1)] = sub.get(m.group(1), 0) + 1
b = open(glob.glob(f"{R}/busco_proteins/short_summary*.txt")[0]).read()
m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%", b)
L = [len(seq[t]) for t in t1]
prefix = min(t1).split("_")[0]
row = {"species": sp, "strain": strain, "locus_tag_prefix": prefix,
       "first_locus_tag": min(t1)[:-3], "last_locus_tag": max(t1)[:-3],
       "raw_models_genes": cs.get("input_genes", ""), "genes": len(t1), "transcripts": len(seq),
       "transcripts_per_gene": round(len(seq) / len(t1), 2),
       "partial_genes": sum(1 for t in t1 if "partial=" in hdr[t]),
       "mean_protein_aa": round(st.mean(L)), "median_protein_aa": round(st.median(L)),
       "BUSCO_C": m.group(1), "BUSCO_S": m.group(2), "BUSCO_D": m.group(3), "BUSCO_F": m.group(4), "BUSCO_M": m.group(5)}
row.update({k: v for k, v in cs.items() if k.startswith("removed_")})
row.update({f"subgenome_{k}": v for k, v in sorted(sub.items())})
open(f"{R}/{sp}.annotation_summary.tsv", "w").write("\t".join(row) + "\n" + "\t".join(map(str, row.values())) + "\n")
for k, v in row.items(): print(f"  {k:28s} {v}")
PY
( cd "${R}" && md5sum "${SP}".genome.fasta "${SP}".gff3 "${SP}".gtf "${SP}".proteins.faa "${SP}".proteins.longest.faa \
                      "${SP}".cds.fna "${SP}".transcripts.fna > md5sums.txt )
mark_done "${R}"
python3 "${BIN}/qc_figures.py" annotation --ref-dir "${D06}" --fig-dir "${FIG}" --species "${RNASEQ_SPECIES[@]}"
echo "Figure: ${FIG}/Fig_annotation_qc.png    Files: ${R}"
