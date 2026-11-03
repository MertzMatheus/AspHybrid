# Publication figures

| Script | Figures | Input |
|---|---|---|
| `paper_figures.R` | Fig. 1–5 of the article | Supplementary Tables 2 and 3 |
| `enrichment_figures.R` | Fig. 6 (enrichment) | curated enrichment sheet of Supplementary Table 3 |

## paper_figures.R

```bash
Rscript paper_figures.R Supplementary_Tables_2.xlsx Supplementary_Tables_3.xlsx [output folder]
```

| Figure | Panels |
|---|---|
| `Fig1_genomes` | a genome size and gene number; b BUSCO of genomes and proteomes; c *A. latus* S/Q subgenomes; d orthology of *A. fumigatus* genes with the other species; e orthogroups shared (UpSet) |
| `Fig2_rnaseq` | a alignment and gene assignment per library; b up/down DEGs per strain |
| `Fig3_cross_species` | a DEGs shared by the three species or lineage-specific; b Jaccard similarity of the responses (same direction); c DEGs shared between species (UpSet, *A. fumigatus* orthologues) |
| `Fig4_genes` | log2FC heatmap: conserved core (same direction in all six strains), *A. fumigatus*-only and *Nidulantes*-only genes |
| `Fig5_hybrid` | a subgenome of *A. latus* DEGs vs genome; b response of S–Q homeolog pairs |

DEGs of *A. spinulosporus* and *A. latus* are compared through their *A. fumigatus*
orthologues (column `A_fumigatus_orthologs` of the DEG sheets), on the 8,223
*A. fumigatus* genes that have orthologues in both other species. The folder also
receives one CSV per panel (the data as plotted) and `figure_numbers.txt`, with
every number quoted in the Results. Needs R ≥ 4.1 with readxl, ggplot2, patchwork
and scales.

## enrichment_figures.R

`enrichment_figures.R` draws the enrichment figures of the article from the
curated enrichment panel (a sheet with the columns of
`enrichment_significant.tsv`: species, strain, direction, collection, term,
name, DEGs_in_term, DEGs_annotated, term_size, universe, gene_ratio,
fold_enrichment, pvalue, padj, significant, genes).

```bash
Rscript enrichment_figures.R <table.xlsx> [sheet] [output folder]
```

| Figure | Content |
|---|---|
| `Fig_enrichment_dotplot_all` | terms (rows, grouped by collection) × species › strain › up/down; size = gene ratio, colour = −log10 adjusted p |
| `Fig_enrichment_dotplot_<species>` | the same, one species per figure |
| `Fig_enrichment_bars_<species>` | diverging bars per strain (up right, down left), number of DEGs per term |
| `Fig_enrichment_counts` | enriched terms per strain, direction and collection |

Output: PDF (vector), PNG 300 dpi and TIFF 600 dpi (LZW), plus
`enrichment_panel_plotted.csv`. Decimal commas ("2,826") are read correctly.
Needs R ≥ 4.1 with readxl, ggplot2 and scales (`envs/detools.yml`).
