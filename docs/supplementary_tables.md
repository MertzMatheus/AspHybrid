# Supplementary tables

`pipeline/09_summary_tables.sh` collects every statistics table into
`09_tables/Supplementary_Tables.xlsx`, one sheet per table. Run it last; tables
of steps not yet run are skipped. Each species also gets its own complete
workbook from step 11 (`Supplementary_DE_<sp>.xlsx`).

| Sheet | Content | Source file | Step |
|---|---|---|---|
| S1 Genome QC | assembly size, N50, L50, N90, GC, N, BUSCO | `01_genome_qc/genome_qc_summary.tsv` | 01 |
| S2 Repeats | masked fraction by repeat class | `02_repeat_masking/*/*.repeat_summary.tsv` | 02 |
| S3 A latus subgenomes | genome fraction S / Q / ambiguous, BUSCO per subgenome | `05_latus_subgenomes/subgenome_summary.tsv` | 05 |
| S4 A latus contigs | majority subgenome and switches per contig | `05_latus_subgenomes/A_latus.contigs.tsv` | 05 |
| S5 Annotation summary | genes, isoforms, partial genes, protein length, BUSCO, curation filters | `06_reference_files/annotation_summary_all.tsv` | 06 |
| S6 Validation | consistency checks of the reference files | `06_reference_files/*/validation.tsv` | 06 |
| S7 Orthology statistics | % of AFUA genes per relation (1:1, 1:2 …) | `07_orthology/ortholog_statistics.tsv` | 07 |
| S8 Orthologs A fumigatus | AFUA → ASPN / ALAT orthologues | `07_orthology/Table_orthologs_Afumigatus.tsv` | 07 |
| S8b Pairs Afum-Alat | one row per AFUA–ALAT pair, with subgenome | `07_orthology/Table_pairs_Afumigatus_Alatus.tsv` | 07 |
| S9 Gene sharing | shared and species-specific genes per species | `07_orthology/sharing/sharing_by_species.tsv` | 07 |
| S10 Sharing patterns | genes and orthogroups per presence/absence pattern | `07_orthology/sharing/sharing_patterns.tsv` | 07 |
| S11 Read mapping | reads, Q30, GC, duplication, alignment rate per library | `03_rnaseq_mapping/*/*.mapping_rates.tsv` | 03 |
| S12 Counting stats | reads assigned to genes, multi-mapping, no feature | `08_count_matrices/*/*.mapping_stats.tsv` | 08 |
| S13 Functional annotation | % of genes with eggNOG, GO, KEGG, COG, CAZy, Pfam | `10_functional_annotation/functional_annotation_summary.tsv` | 10 |
| S14 COG categories | genes per COG category and species | `10_functional_annotation/COG_summary.tsv` | 10 |
| S15 DE summary | up / down DEGs per strain | `11_differential_expression/*/DE_summary.tsv` | 11 |
| S16 DEG sharing strains | DEGs shared between strains of a species | `11_differential_expression/*/DEG_sharing.tsv` | 11 |
| S17 Enrichment summary | significant terms per strain, direction and collection | `11_differential_expression/*/enrichment_summary.tsv` | 11 |
| S18 Enrichment significant | every significant term with its genes | `11_differential_expression/*/enrichment_significant.tsv` | 11 |
| S19 Software | version of every tool and database | `09_tables/software_versions.tsv` | 09 |

Excel cells hold at most 32,767 characters; longer cells (very long gene lists)
are truncated in the workbook only, and the TSV files keep the full content.
