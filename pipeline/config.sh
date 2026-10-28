# ============================================================
# config.sh  -  the only file you need to edit
# ============================================================

# ---------------- inputs ----------------
BASE="/path/to/project"                     # project folder on a large disk (edit this)
GENOMES="${BASE}/Genomes"                    # Genomes/<species>/<one .fasta>
# sample sheet: by default the one next to this config.sh (wherever the pipeline folder is)
SAMPLES="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/rnaseq_samples.tsv"
DB="${BASE}/db"

# RNA-seq design for differential expression (one row per sample; sample_old = name in rnaseq_samples.tsv)
# default: the design of the study, in config/ of the repository
METADATA="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/RNAseq_metadata.csv"

# ---------------- outputs ----------------
OUT="${BASE}/final_results"
# optional earlier run: its repeat masking is reused when the genome is identical (ignored if absent)
PREVIOUS_RESULTS="${BASE}/Results"

THREADS=28

# ---------------- species ----------------
ALL_SPECIES=( A_fumigatus A_spinulosporus A_latus A_quadrilineatus )
RNASEQ_SPECIES=( A_fumigatus A_spinulosporus A_latus )      # annotated + quantified
# A_quadrilineatus: genome QC only (reference for the A. latus subgenome phasing)

declare -A STRAIN=(
    [A_fumigatus]="Af293"
    [A_spinulosporus]="NRRL 2395"
    [A_latus]="NRRL 200"
    [A_quadrilineatus]="NRRL 201"
)
# locus_tag prefixes (NCBI may assign different ones; change here, re-run step 06)
declare -A PREFIX=(
    [A_fumigatus]="AFUA"
    [A_spinulosporus]="ASPN"
    [A_latus]="ALAT"
)
# A. fumigatus annotation:
#   official = NCBI RefSeq Af293 annotation, official locus_tags AFUA_1G00100 ...
#              (placed on the project genome; needs it to be the Af293 RefSeq assembly)
#   braker   = de novo BRAKER3 models like the other species (locus_tags AFUM_...)
AF_ANNOTATION="official"
AF293_ACC="GCF_000002655.1_ASM265v1"
AF293_DIR="${DB}/Af293_RefSeq"
AF293_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/${AF293_ACC}"
LATUS_PARENT_S="A_spinulosporus"
LATUS_PARENT_Q="A_quadrilineatus"

# ---------------- databases (already in db/) ----------------
BUSCO_LINEAGE_DIR="${DB}/busco/lineages/eurotiales_odb12.2"
BUSCO_LINEAGE_TGZ="${DB}/busco/lineages/eurotiales_odb12.2.tar.gz"
ORTHODB_FUNGI="${DB}/orthodb/Fungi.fa.gz"

EGGNOG_DATA_DIR="${DB}/eggnog"          # eggnog.db, eggnog_proteins.dmnd, eggnog.taxa.db ...
GO_OBO="${DB}/go/go-basic.obo"           # downloaded once if missing
KEGG_PATHWAY_NAMES="${DB}/kegg/pathway_names.tsv"   # downloaded once (rest.kegg.jp) if missing

# ---------------- functional annotation / DE settings ----------------
EGGNOG_TAX_SCOPE="auto"                  # eggNOG default: narrowest informative taxonomic level
EGGNOG_GO_EVIDENCE="non-electronic"      # eggNOG default; "all" also transfers IEA GO terms
EGGNOG_DBMEM=1                           # 1 = load eggnog.db in RAM (~45 GB of RAM): much faster; 0 on smaller machines
DE_PADJ=0.05                             # DEG: BH-adjusted p < DE_PADJ ...
DE_LFC=1                                 # ... and |shrunken log2 fold change| >= DE_LFC (2-fold)
DE_MIN_COUNT=10                          # pre-filter: >= DE_MIN_COUNT reads in >= 3 samples
ENRICH_PADJ=0.05                         # enrichment: BH-adjusted p < ENRICH_PADJ
ENRICH_MIN_SET=5; ENRICH_MAX_SET=500     # gene-set size limits (genes in the universe)

# ---------------- software ----------------
SOFT="${BASE}/software"
ENV_TOOLS="${SOFT}/envs/AnnotTools"           # minimap2 hisat2 samtools subread busco ...
ENV_ORTHO="Hibrid"                            # OrthoFinder
ENV_DE="${SOFT}/envs/DEtools"                 # R + DESeq2 + apeglm
IMG_BRAKER="teambraker/braker3@sha256:5f8b3c508a9fe1bbc2e9a74dcc013eeed82f91dd5945adca7823514d9c8aecf8"
IMG_TETOOLS="dfam/tetools@sha256:7fe6afce694d0759b56d65b623a9d684269e8f82667f0f2eb79fca6a63c51189"

# BRAKER3 RNA-seq evidence: one replicate per group is enough for gene
# structure (diversity of conditions matters, not depth); keeps the BAM small.
BRAKER_EVIDENCE_REPLICATE=1
# a single spaln alignment (ProtHint) running longer than this is killed (see 04)
SPALN_MAX_HOURS=2
# RNA-seq evidence for BRAKER: max read pairs per sample (random subsample; 0 = all reads)
EVIDENCE_MAX_PAIRS=15000000

export CONDA_PKGS_DIRS="${SOFT}/pkgs"
export TMPDIR="${BASE}/tmp"
export XDG_CACHE_HOME="${BASE}/.cache"
mkdir -p "${TMPDIR}" "${XDG_CACHE_HOME}" 2>/dev/null || true
