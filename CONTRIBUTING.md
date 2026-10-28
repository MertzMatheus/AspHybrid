# Contributing

Thank you for your interest in this workflow.

## Reporting a problem
Open a [GitHub issue](../../issues) with:
1. the step and command you ran (e.g. `bash 04_gene_prediction.sh A_latus`);
2. the last 30 lines of the terminal output or of the step's log file;
3. your operating system, number of threads and RAM;
4. the tool versions (`pipeline/09_summary_tables.sh` writes them to `09_tables/software_versions.tsv`).

Please check [docs/troubleshooting.md](docs/troubleshooting.md) first: most
known problems and their fixes are listed there.

## Suggesting a change
1. Fork the repository and create a branch (`git checkout -b fix-short-name`).
2. Keep one step per script and all user settings in `pipeline/config.sh`.
3. Test the changed step on at least one species and describe the test in the pull request.
4. Update the documentation in `docs/` and `CHANGELOG.md` when behaviour or outputs change.

## Using the workflow for other species
The pipeline is written for this study but is general: add the species to
`ALL_SPECIES` / `RNASEQ_SPECIES`, give it a locus_tag prefix in `PREFIX`, and
provide its genome and sample sheet rows. Subgenome assignment (step 05) is
specific to allodiploid hybrids with two available parental genomes.
