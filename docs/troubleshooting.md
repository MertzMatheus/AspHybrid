# Troubleshooting

Problems met while running the workflow on the study data, and their fixes.

| Message or symptom | Cause | Fix |
|---|---|---|
| `No space left on device` while creating environments | `/home` partition full | keep environments, package cache and `TMPDIR` on the data disk (`CONDA_PKGS_DIRS`, `TMPDIR`, `XDG_CACHE_HOME` in `config.sh`) |
| `ERROR: missing input: …/rnaseq_samples.tsv` | sample sheet not next to `config.sh` | put `rnaseq_samples.tsv` in `pipeline/` or set `SAMPLES` |
| `missing FASTQ for <sample>` | wrong R1/R2 path in the sample sheet | use the full path of the clean reads (`*.mRNA_fwd.fq.gz`, `*.mRNA_rev.fq.gz`) |
| `SyntaxError: f-string expression part cannot include a backslash` | Python < 3.12 with an older helper | update to the current scripts (compatible with Python ≥ 3.10) |
| `samtools sort: failed writing to …/bam/<sample>.bam.tmp.bam: No such file or directory` | output folder removed, or two runs of the same species at once | the mapping step now refuses a second run of the same species and recreates its folders; run it again, finished samples are skipped |
| BRAKER3 stops at `samtools sort … merging from N files` with no other message | RNA-seq evidence BAM too large for BRAKER's internal sort (disk or memory) | evidence is now limited to mapped primary reads and `EVIDENCE_MAX_PAIRS` per library; run step 04 again |
| BRAKER3 apparently frozen for hours, one CPU at 100 %; `docker top` shows a single old `spaln` | ProtHint aligning one very long protein (e.g. a PKS/NRPS) | step 04 kills any Spaln alignment older than `SPALN_MAX_HOURS` (logged in `spaln_watchdog.log`); in a run already going, `kill <PID>` of that `spaln` only |
| BRAKER3 GeneMark-ETP step takes many hours | 24 million OrthoDB Fungi proteins | expected (6–20 h per genome); check progress with `docker stats` and `docker top braker_<sp>` |
| `D10: unbound variable` | old `lib/common.sh` | update `lib/common.sh` (current steps 09–11 also work with the old one) |
| `ZeroDivisionError` in `gene_sharing.py` | OrthoFinder run with different proteome names than expected | update step 07; it reuses OrthoFinder only when run on exactly the current proteins and otherwise archives the old run in `07_orthology/previous_<date>/` |
| `RefSeq sequences with genes are not identical …` in step 06 (*A. fumigatus*) | project genome is not the Af293 RefSeq assembly | use the Af293 assembly, or set `AF_ANNOTATION="braker"` |
| BUSCO tries to download the lineage and is slow | online mode | the pipeline runs BUSCO `--offline` with the local `eurotiales_odb12.2` |
| `docker.io/interpro/interproscan:latest not found` | the image has no `latest` tag | not needed by this workflow; if used, pin a version tag |
| eggNOG-mapper very slow | database read from disk | `EGGNOG_DBMEM=1` (needs ~45 GB RAM) |

## Monitoring long steps

```bash
tail -f logs_04_A_latus.log                       # progress messages
docker stats --no-stream                          # BRAKER3 / RepeatModeler CPU and memory
docker top braker_A_latus -eo pid,etimes,comm     # which program is running and for how long
df -h /path/to/project                            # free disk
```
