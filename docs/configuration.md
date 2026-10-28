# Configuration (`pipeline/config.sh`)

`config.sh` is the only file to edit. Every step reads it.

## Paths

| Variable | Meaning | Default |
|---|---|---|
| `BASE` | project folder on the large disk | `/path/to/project` |
| `GENOMES` | one sub-folder per species with one FASTA | `${BASE}/Genomes` |
| `SAMPLES` | sample sheet | `rnaseq_samples.tsv` next to `config.sh` |
| `METADATA` | design table for differential expression | `config/RNAseq_metadata.csv` of the repository |
| `DB` | databases | `${BASE}/db` |
| `OUT` | all results | `${BASE}/final_results` |
| `PREVIOUS_RESULTS` | optional earlier run; its repeat masking is reused when the genome is identical | `${BASE}/Results` |
| `THREADS` | threads per step | `28` |

## Species and identifiers

| Variable | Meaning |
|---|---|
| `ALL_SPECIES` | every genome (step 01) |
| `RNASEQ_SPECIES` | species annotated and quantified (steps 02–11) |
| `STRAIN` | assembly strain shown in tables (Af293, NRRL 2395, NRRL 200, NRRL 201) |
| `PREFIX` | locus_tag prefix per species (`AFUA`, `ASPN`, `ALAT`) |
| `AF_ANNOTATION` | `official` = RefSeq Af293 annotation (AFUA_ locus_tags); `braker` = de novo like the other species (AFUM_) |
| `LATUS_PARENT_S`, `LATUS_PARENT_Q` | parental genomes for the *A. latus* subgenome assignment |

## Databases and software

| Variable | Meaning |
|---|---|
| `BUSCO_LINEAGE_DIR`, `BUSCO_LINEAGE_TGZ` | eurotiales_odb12.2 (folder, or the .tar.gz unpacked on first use) |
| `ORTHODB_FUNGI` | OrthoDB v12 Fungi proteins (BRAKER3 evidence) |
| `EGGNOG_DATA_DIR` | eggNOG database folder |
| `GO_OBO`, `KEGG_PATHWAY_NAMES` | downloaded once if missing |
| `AF293_ACC`, `AF293_DIR`, `AF293_URL` | official Af293 annotation (downloaded once) |
| `ENV_TOOLS`, `ENV_DE`, `ENV_ORTHO` | micromamba environments |
| `IMG_BRAKER`, `IMG_TETOOLS` | Docker images pinned by digest |

## Analysis parameters

| Variable | Default | Effect |
|---|---|---|
| `BRAKER_EVIDENCE_REPLICATE` | `1` | replicate of each strain × condition used as BRAKER3 RNA-seq evidence |
| `EVIDENCE_MAX_PAIRS` | `15000000` | max read pairs per library in the BRAKER3 evidence (0 = all) |
| `SPALN_MAX_HOURS` | `2` | a single Spaln alignment running longer is stopped (ProtHint) |
| `EGGNOG_TAX_SCOPE` | `auto` | eggNOG taxonomic scope |
| `EGGNOG_GO_EVIDENCE` | `non-electronic` | `all` also transfers electronically inferred GO terms |
| `EGGNOG_DBMEM` | `1` | load eggnog.db in RAM (~45 GB); set 0 on smaller machines |
| `DE_PADJ` | `0.05` | DEG: adjusted p below this ... |
| `DE_LFC` | `1` | ... and absolute shrunken log2 fold change at least this (2-fold) |
| `DE_MIN_COUNT` | `10` | pre-filter: at least this many reads in at least 3 samples |
| `ENRICH_PADJ` | `0.05` | enriched term: adjusted p below this |
| `ENRICH_MIN_SET`, `ENRICH_MAX_SET` | `5`, `500` | gene-set size in the universe (max applied to GO only) |

Changing a threshold only requires re-running the steps that use it
(e.g. `DE_LFC` → step 11, then step 09).
