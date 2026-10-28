# Outputs

Everything is written to `OUT` (`final_results/`), one folder per step.

```
final_results/
├── 01_genome_qc/
│   ├── <sp>/<sp>.genome.fasta            clean genome used by every later step
│   ├── <sp>/<sp>.header_map.tsv           new ID ↔ original FASTA header
│   ├── <sp>/busco/                        BUSCO genome mode
│   └── genome_qc_summary.tsv              size, N50, L50, N90, GC, N, BUSCO (Table S1)
├── 02_repeat_masking/
│   └── <sp>/<sp>.softmasked.fasta, <sp>-families.fa, <sp>.repeat_summary.tsv (Table S2)
├── 03_rnaseq_mapping/
│   └── <sp>/bam/<sample>.bam(.bai), qc/ (fastp, HISAT2, flagstat), <sp>.mapping_rates.tsv (Table S11)
├── 04_gene_prediction/
│   └── <sp>/braker/braker.gff3            raw BRAKER3 models
│       <sp>/evidence_subsampling.tsv      read pairs used as RNA-seq evidence
│       <sp>/spaln_watchdog.log            Spaln alignments stopped after SPALN_MAX_HOURS
├── 05_latus_subgenomes/
│   ├── A_latus.windows.tsv                identity to each parent per 5-kb window
│   ├── A_latus.blocks.bed                 S / Q / ambiguous blocks
│   ├── A_latus.contigs.tsv                majority subgenome per contig (Table S4)
│   └── subgenome_summary.tsv              genome fraction and BUSCO per subgenome (Table S3)
├── 06_reference_files/
│   ├── <sp>/<sp>.gff3                     curated gene models; gene ID = locus_tag
│   ├── <sp>/<sp>.gtf                      same models for read counting; gene_id = locus_tag
│   ├── <sp>/<sp>.proteins.faa             all isoforms (<locus_tag>-T1, -T2 ...)
│   ├── <sp>/<sp>.proteins.longest.faa     one protein per gene (orthology, eggNOG, BUSCO)
│   ├── <sp>/<sp>.cds.fna, <sp>.transcripts.fna
│   ├── <sp>/<sp>.id_map.tsv               locus_tag ↔ BRAKER / RefSeq IDs, coordinates, subgenome
│   ├── <sp>/<sp>.removed_models.tsv       models removed by curation and why
│   ├── <sp>/validation.tsv, md5sums.txt, busco_proteins/
│   ├── A_fumigatus/AFUM_to_AFUA.tsv       de novo BRAKER3 gene → official AFUA gene
│   ├── A_fumigatus/AFUM_vs_AFUA_summary.tsv   % of official genes recovered (benchmark)
│   ├── A_fumigatus_BRAKER_AFUM/           de novo A. fumigatus annotation (kept)
│   └── annotation_summary_all.tsv         genes, isoforms, protein length, BUSCO (Table S5)
├── 07_orthology/
│   ├── Table_orthologs_Afumigatus.tsv     one row per AFUA gene → ASPN and ALAT orthologues (Table S8)
│   ├── Table_pairs_Afumigatus_Alatus.tsv  one row per AFUA–ALAT pair, with subgenome
│   ├── Table_pairs_Afumigatus_Aspinulosporus.tsv
│   ├── ortholog_statistics.tsv            % of AFUA genes per relation type (Table S7)
│   ├── tables/correspondence_<sp>.tsv     orthologues anchored on each species
│   ├── sharing/                           shared and species-specific genes (Tables S9–S10)
│   └── OrthoFinder/Results_*/             raw OrthoFinder output
├── 08_count_matrices/
│   └── <sp>/<sp>.counts.tsv               raw counts, genes × samples (DESeq2 input)
│       <sp>/<sp>.TPM.tsv, <sp>.CPM.tsv, <sp>.counts.annotated.tsv
│       <sp>/<sp>.sample_metadata.tsv, <sp>.mapping_stats.tsv (Table S12), strandedness.txt
├── 09_tables/
│   ├── Supplementary_Tables.xlsx          all tables, one sheet each
│   └── software_versions.tsv              versions of every tool (Table S19)
├── 10_functional_annotation/
│   ├── <sp>/<sp>.emapper.annotations      raw eggNOG-mapper output
│   ├── <sp>/<sp>.gene_annotation.tsv      one row per gene: product, eggNOG, COG, GO, KEGG, CAZy, Pfam, orthologues
│   ├── <sp>/genesets/<sp>.<collection>.TERM2GENE.tsv / TERM2NAME.tsv   gene sets (GO_BP, GO_MF, GO_CC, KEGG_pathway, COG, KO, Pfam, CAZy)
│   ├── functional_annotation_summary.tsv  annotation coverage (Table S13)
│   └── COG_summary.tsv                    COG categories per species (Table S14)
├── 11_differential_expression/
│   └── <sp>/deseq2/DE_<strain>.tsv        DESeq2 results per strain; vst.tsv, normalized_counts.tsv
│       <sp>/DE_<strain>.annotated.tsv     all tested genes + call (up/down/ns) + annotation
│       <sp>/DEGs_<strain>.tsv             significant genes only
│       <sp>/DE_summary.tsv                up / down per strain (Table S15)
│       <sp>/DEG_sharing.tsv, DEG_strain_membership.tsv   DEGs shared between strains (Table S16)
│       <sp>/enrichment_significant.tsv, enrichment_all_tests.tsv, enrichment_summary.tsv (Tables S17–S18)
│       <sp>/Supplementary_DE_<sp>.xlsx    everything for one species in one workbook
└── figures/
    ├── Fig_genome_qc, Fig_latus_subgenomes, Fig_annotation_qc
    ├── Fig_orthologs, Fig_gene_sharing
    ├── Fig_RNAseq_QC_<sp>, Fig_functional_annotation
    └── Fig_DE_<sp>, Fig_enrichment_<sp>     (PNG 300 dpi + PDF)
```

## Identifiers

| Species | Genes | Isoforms / proteins | Source |
|---|---|---|---|
| *A. fumigatus* | `AFUA_1G00100` | `AFUA_1G00100-T1` | official RefSeq Af293 locus_tags |
| *A. spinulosporus* | `ASPN_000010` | `ASPN_000010-T1` | this work, genome order, steps of 10 |
| *A. latus* | `ALAT_000010` | `ALAT_000010-T1` | this work, genome order, steps of 10 |

`-T1` is always the longest protein. The same locus_tag is used in the GFF3,
GTF, FASTA files, count matrices, orthology, eggNOG and DE tables, so any two
tables can be joined on it.

## Main columns

**`DE_<strain>.annotated.tsv`**: `locus_tag`, `baseMean`, `log2FC` (apeglm-shrunken),
`log2FC_MLE` (unshrunken), `lfcSE`, `stat`, `pvalue`, `padj`, `mean_norm_DMEN`,
`mean_norm_Macrophage`, `DE` (up / down / ns / not_tested), then the gene annotation
and orthologue columns.

**`enrichment_significant.tsv`**: `species`, `strain`, `direction`, `collection`,
`term`, `name`, `DEGs_in_term` (k), `DEGs_annotated` (n), `term_size` (K),
`universe` (N), `gene_ratio` (k/n), `fold_enrichment` ((k/n)/(K/N)), `pvalue`, `padj`,
`significant`, `genes`, `gene_names`.

**`Table_orthologs_Afumigatus.tsv`**: `A_fumigatus`, `orthogroup`, `A_spinulosporus`,
`A_spinulosporus_relation`, `A_latus`, `A_latus_n`, `A_latus_relation`,
`A_latus_subgenome`, `A_latus_S`, `A_latus_Q`, `one_to_one_all`,
`A_fumigatus_1to1_As_and_1to2_Al`.
