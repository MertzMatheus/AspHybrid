# Software environments

| File | Environment | Used by |
|---|---|---|
| `annottools.yml` | `AnnotTools` (path in `ENV_TOOLS`) | pipeline steps 01, 03, 05, 06, 08, 09, 10 and all Python helpers |
| `detools.yml` | `DEtools` (path in `ENV_DE`) | step 11 (DESeq2) and `figures/enrichment_figures.R` |
| `orthofinder.yml` | `Hibrid` (name in `ENV_ORTHO`) | step 07 |
| `preprocessing.yml` | `rnaseq_sa_ca` | `preprocessing/` |

Containers (pinned by digest in `pipeline/config.sh`):

| Image | Tools | Steps |
|---|---|---|
| `teambraker/braker3@sha256:5f8b3c50…aecf8` | BRAKER 3.0.8, GeneMark-ETP, AUGUSTUS, ProtHint, Spaln | 04 |
| `dfam/tetools@sha256:7fe6afce…51189` | RepeatModeler 2.0.9, RepeatMasker 4.2.6 | 02 |

The YAML files list packages without exact versions (except BUSCO and
OrthoFinder) so they solve on new systems. For full reproducibility, export the
environments you actually used and commit the lock files:

```bash
micromamba env export -p /path/to/software/envs/AnnotTools > envs/annottools.lock.yml
micromamba env export -p /path/to/software/envs/DEtools    > envs/detools.lock.yml
micromamba env export -n Hibrid                            > envs/orthofinder.lock.yml
micromamba env export -n rnaseq_sa_ca                      > envs/preprocessing.lock.yml
```

`pipeline/09_summary_tables.sh` also writes the version of every tool to
`09_tables/software_versions.tsv` (Supplementary Table S19).
