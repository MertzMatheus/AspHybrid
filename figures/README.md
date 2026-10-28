# Publication figures

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
