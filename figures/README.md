# Publication figures

| Script | Figures | Input |
|---|---|---|
| `paper_figures.R` | every panel of Fig. 1–4 and Supplementary Figs 1–2, one file per panel | Supplementary Tables 2 and 3; optionally the pipeline's `final_results/11_differential_expression` |
| `enrichment_figures.R` | alternative enrichment dot plots of the curated panel | curated enrichment sheet of Supplementary Table 3 |

## paper_figures.R

```bash
Rscript paper_figures.R Supplementary_Tables_2.xlsx Supplementary_Tables_3.xlsx \
        /path/to/project/final_results/11_differential_expression  [output folder]
```

Use `NA` as the third argument to skip the PCA and volcano plots (they need
`<species>/deseq2/vst.tsv`, `samples_used.tsv` and `DE_<strain>.tsv` from step 11).

| Figure | Panels (one file each, in `Paper_figures/FigN/`) |
|---|---|
| Fig. 1 – annotation and gene sharing | a genes and locus-tag prefixes; b BUSCO of genomes and proteomes; c functional annotation coverage; d COG categories; e orthology of *A. fumigatus* genes; f orthogroups shared (UpSet); g gene sharing per proteome |
| Fig. 2 – RNA-seq | a PCA per species; b DEGs per strain; c volcano per strain; d DEGs shared between strains of a species (UpSet) |
| Fig. 3 – enrichment | a top terms per strain (up and down); b terms recurring across strains |
| Fig. 4 – comparison through *A. fumigatus* orthologues | a DEGs in shared or lineage-specific genes; b DEGs of each strain against *A. fumigatus*; c Jaccard similarity; d DEGs shared between species (UpSet); e conserved and lineage-specific genes (heatmap) |
| Supplementary | S1 read alignment and assignment; S2 *A. latus* subgenomes (genome, DEGs, homeolog pairs) |

DEGs of *A. spinulosporus* and *A. latus* are compared through their *A. fumigatus*
orthologues (column `A_fumigatus_orthologs` of the DEG sheets), on the 8,223
*A. fumigatus* genes that have orthologues in both other species. Each panel is
saved as PDF and PNG (add `"tiff"` to `FORMATS` for 600 dpi TIFF) with a CSV of the
plotted data; `figure_numbers.txt` lists every number quoted in the Results.
Needs R ≥ 4.1 with readxl, ggplot2, patchwork and scales.

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
