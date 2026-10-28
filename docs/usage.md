# Usage

Each step is a separate script, run by hand in the order below. A step checks
its inputs and tells you which earlier step to run if something is missing;
finished work is skipped when a step is run again, so re-running after an
error is safe. Long steps are best run in the background:

```bash
nohup bash 04_gene_prediction.sh A_latus > logs_04_A_latus.log 2>&1 &
tail -f logs_04_A_latus.log
```

`<sp>` = `A_fumigatus`, `A_spinulosporus` or `A_latus`.

## 1. Read pre-processing (`preprocessing/`)

```bash
export BASE=/path/to/project
bash preprocessing/01_quality.sh      # fastp: adapters, Q20, reads >= 30 nt  -> 02_qc/
bash preprocessing/02_remove_rRNA.sh  # SortMeRNA vs SILVA rRNA              -> 03_rmRNA/*.mRNA_fwd|rev.fq.gz
```

Both scripts carry SLURM headers (`sbatch` on a cluster) and also run with
`bash` on a workstation.

## 2. Fungal pipeline (`pipeline/`)

```bash
cd pipeline

# 2.1 genome quality (minutes each)
for sp in A_fumigatus A_spinulosporus A_latus A_quadrilineatus; do bash 01_genome_qc.sh $sp; done

# 2.2 read mapping and repeat masking (independent; can run at the same time)
for sp in A_fumigatus A_spinulosporus A_latus; do bash 03_map_rnaseq.sh $sp; done
for sp in A_spinulosporus A_latus; do nohup bash 02_repeat_masking.sh $sp > logs_02_$sp.log 2>&1; done &

# 2.3 A. latus subgenomes (needs 01 for A_latus, A_spinulosporus, A_quadrilineatus)
bash 05_latus_subgenomes.sh

# 2.4 gene prediction (A. spinulosporus and A. latus; 6-24 h each)
nohup bash 04_gene_prediction.sh A_spinulosporus > logs_04_A_spinulosporus.log 2>&1 &
nohup bash 04_gene_prediction.sh A_latus         > logs_04_A_latus.log 2>&1 &

# 2.5 reference files (A_fumigatus: official AFUA annotation, no step 04 needed)
for sp in A_fumigatus A_spinulosporus A_latus; do bash 06_reference_files.sh $sp; done

# 2.6 orthology, counts, functional annotation, differential expression
bash 07_orthology.sh
for sp in A_fumigatus A_spinulosporus A_latus; do bash 08_count_matrix.sh $sp; done
for sp in A_fumigatus A_spinulosporus A_latus; do bash 10_eggnog.sh $sp; done
for sp in A_fumigatus A_spinulosporus A_latus; do bash 11_differential_expression.sh $sp; done

# 2.7 all supplementary tables in one workbook
bash 09_summary_tables.sh
```

Two jobs that each use all threads (e.g. two BRAKER3 runs) slow each other
down; run them one after the other or lower `THREADS` in `config.sh`.

### Dependencies between steps

| Step | Needs |
|---|---|
| 02, 03 | 01 of the same species |
| 04 | 02 and 03 of the same species |
| 05 | 01 of *A. latus*, *A. spinulosporus*, *A. quadrilineatus* |
| 06 | 04 (*A. spinulosporus*, *A. latus*); 05 (*A. latus*); 01 only (*A. fumigatus*, official annotation) |
| 07 | 06 of the three species |
| 08 | 03 and 06 (07 optional: adds orthologue columns) |
| 10 | 06 (07 optional: adds orthologue columns) |
| 11 | 08 and 10 |
| 09 | anything; missing tables are skipped |

### What to check after each step

| After | Look at | Expected |
|---|---|---|
| 01 | `01_genome_qc/genome_qc_summary.tsv` | haploid genomes: BUSCO C > 95 %, D low; *A. latus*: D high, size ≈ sum of the parents |
| 02 | `02_repeat_masking/repeat_summary.tsv` | a few % of each genome masked |
| 03 | `03_rnaseq_mapping/<sp>/<sp>.mapping_rates.tsv` | overall alignment > 85–90 % |
| 04 | last lines of the log | ~10–12 k genes per haploid genome, about twice that in *A. latus* |
| 05 | `figures/Fig_latus_subgenomes.png` | two clear peaks of identity difference; each subgenome BUSCO C high, D low |
| 06 | `06_reference_files/<sp>/validation.tsv` | status `OK`, re-translation ≥ 98 % identical |
| 07 | `07_orthology/ortholog_statistics.tsv` | most AFUA genes 1:1 with ASPN and 1:2 with ALAT |
| 08 | `figures/Fig_RNAseq_QC_<sp>.png` | replicates cluster by strain and condition |
| 10 | `10_functional_annotation/functional_annotation_summary.tsv` | most genes with an eggNOG hit |
| 11 | `11_differential_expression/<sp>/DE_summary.tsv`, `figures/Fig_DE_<sp>.png` | DEGs per strain, PCA separates macrophage from DMEM |

## 3. Publication figures (`figures/`)

After curating the enrichment terms to show in the paper (one sheet with the
columns of `enrichment_significant.tsv`):

```bash
Rscript figures/enrichment_figures.R /path/to/Supplementary_Tables_3.xlsx
```

Figures are written next to the table, in `Enrichment_figures/` (PDF, PNG 300 dpi, TIFF 600 dpi).
