# AspHybrid pipeline (steps 01–11)

Genome quality, repeat masking, read mapping, gene prediction, *A. latus*
subgenome assignment, curated reference files, orthology, count matrices,
functional annotation, differential expression and supplementary tables for
*A. fumigatus*, *A. spinulosporus* and *A. latus*.

```
pipeline/
├── config.sh                      paths, species, databases, environments, thresholds (edit this)
├── rnaseq_samples.tsv             sample sheet (template: ../config/rnaseq_samples.example.tsv)
├── 01_genome_qc.sh <sp>           clean genome, assembly statistics, BUSCO
├── 02_repeat_masking.sh <sp>      RepeatModeler2 + RepeatMasker (Docker)
├── 03_map_rnaseq.sh <sp>          HISAT2, fastp report, mapping table
├── 04_gene_prediction.sh <sp>     BRAKER3 ETP (Docker)
├── 05_latus_subgenomes.sh         S / Q assignment of A. latus
├── 06_reference_files.sh <sp>     curated GFF3 / GTF / FASTA with locus_tags, validation, BUSCO
├── 07_orthology.sh                OrthoFinder, AFUA-anchored orthologue tables, gene sharing
├── 08_count_matrix.sh <sp>        featureCounts, count / TPM / CPM matrices, RNA-seq QC
├── 09_summary_tables.sh           Supplementary_Tables.xlsx and software versions (run last)
├── 10_eggnog.sh <sp>              eggNOG-mapper, gene annotation table, gene sets
├── 11_differential_expression.sh <sp>   DESeq2 per strain, enrichment, figures, Excel
├── lib/common.sh                  shared functions and output folders
└── bin/                           helpers called by the steps
    ├── clean_genome.py            clean FASTA + assembly statistics
    ├── qc_figures.py              genome and annotation QC tables and figures
    ├── phase_subgenomes.py        5-kb windows, identity to each parent, S/Q blocks
    ├── subgenome_figure.py        A. latus subgenome figure and summary
    ├── curate_annotation.py       BRAKER3 model curation and locus_tags
    ├── refseq_to_project.py       official Af293 annotation placed on the project genome
    ├── map_afum_afua.py           de novo AFUM genes ↔ official AFUA genes (benchmark)
    ├── gff3_to_gtf.py             GTF with gene_id = locus_tag
    ├── validate_reference.py      consistency checks of the reference files
    ├── build_ortholog_table.py    orthologue tables from OrthoFinder
    ├── ortholog_summary.py        AFUA-anchored tables, statistics and figure
    ├── gene_sharing.py            shared / species-specific genes, Venn
    ├── make_count_matrix.py       count, TPM, CPM matrices and sample metadata
    ├── rnaseq_figures.py          RNA-seq QC panel (read fate, PCA, correlation)
    ├── functional_db.py           eggNOG gene table and GO/KEGG/COG gene sets
    ├── functional_figure.py       annotation coverage and COG figure
    ├── deseq2_strains.R           DESeq2 model and contrasts per strain
    ├── de_analysis.py             DEG calls, enrichment, PCA, volcano, Venn, tables
    └── tables_to_excel.py         TSV tables → Excel workbook
```

Run order and expected results: [../docs/usage.md](../docs/usage.md).
