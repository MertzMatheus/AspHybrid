# Changelog – AspHybrid

All notable changes are listed here. Versions follow [Semantic Versioning](https://semver.org).

## [1.0.0] – [release date]

Version used for the article [Article title].

### Workflow
- Read pre-processing: fastp trimming and SortMeRNA rRNA removal (`preprocessing/`).
- Fungal pipeline, steps 01–11 (`pipeline/`): genome QC, repeat masking, read
  mapping, BRAKER3 gene prediction, *A. latus* subgenome assignment, curated
  reference files, orthology, count matrices, eggNOG annotation, differential
  expression and enrichment, supplementary tables.
- *A. fumigatus* uses the official RefSeq Af293 annotation (AFUA_ locus_tags);
  the de novo BRAKER3 models are kept as a benchmark and linked to AFUA by CDS overlap.
- Publication figures from the curated enrichment panel (`figures/`).

### Robustness (found while running on the real data)
- BRAKER3 RNA-seq evidence limited to mapped primary alignments and at most
  15 M read pairs per library, avoiding a failure of BRAKER's internal `samtools sort`.
- Watchdog that stops any single Spaln alignment running for more than 2 h
  (ProtHint), so BRAKER3 cannot hang on one very long protein.
- One-run-per-species lock in the mapping step; temporary files and output
  folders recreated per sample.
- Python helpers compatible with Python 3.10+.
