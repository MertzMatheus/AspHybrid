# Input tables

## `rnaseq_samples.example.tsv` → copy to `pipeline/rnaseq_samples.tsv`

One row per paired-end library, tab-separated, every column filled
(the header line starts with `#`):

| Column | Meaning | Example |
|---|---|---|
| species | must match a name in `ALL_SPECIES` (config.sh) | `A_latus` |
| sample | library name; becomes the BAM name and the count-matrix column | `R171001` |
| R1, R2 | full paths of the clean reads (output of `preprocessing/`) | `/data/03_rmRNA/A_latus/R171001.mRNA_fwd.fq.gz` |
| strain | strain name | `1710` |
| condition | `control` (fungus alone in DMEM) or `macrophage` | `macrophage` |
| replicate | biological replicate number | `1` |

Used by steps 03 (mapping), 04 (BRAKER3 evidence: replicate 1 of each strain × condition) and 08 (counting).

## `RNAseq_metadata.csv` → design table of the study (`METADATA` in config.sh points here)

Design table for differential expression (step 11), comma-separated. Required columns:

| Column | Meaning |
|---|---|
| sample_old | library name, identical to `sample` in the sample sheet |
| sample_id | readable name used in figures (e.g. `Alat_1710_Macrophage_01`) |
| species_id | `A_fumigatus`, `A_spinulosporus`, `A_latus` |
| strain | strain name; one contrast is made per strain |
| treatment | `DMEN` (control, fungus alone) or `Macrophage` |
| replicate | biological replicate |
| strain_color | colour of the strain in the PCA (hex) |

Other columns (ploidy, colours per group, slots) are carried along but not required.
