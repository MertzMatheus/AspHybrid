# Read pre-processing

Trimming and rRNA removal of the 36 fungal RNA-seq libraries.

| Script | Tool | Settings | Output |
|---|---|---|---|
| `01_quality.sh` | fastp | `--detect_adapter_for_pe --qualified_quality_phred 20 --length_required 30` | `02_qc/<sample>_R1/R2.trim.fastq.gz`, HTML/JSON reports |
| `02_remove_rRNA.sh` | SortMeRNA | SILVA bacterial 16S (id90), 23S (id98), eukaryotic 18S (id95), 28S (id98); `--paired_out` | `03_rmRNA/<sample>.mRNA_fwd.fq.gz`, `.mRNA_rev.fq.gz` (non-rRNA pairs) |

Set the project folder with `BASE` (exported or edited at the top of each script):

```bash
BASE=/path/to/project bash preprocessing/01_quality.sh
BASE=/path/to/project bash preprocessing/02_remove_rRNA.sh
```

Raw reads are expected as `${BASE}/01.RawData/<sample>/<sample>_1.fq.gz` and
`_2.fq.gz`; the SortMeRNA databases (`silva-bac-16s-id90.fasta`,
`silva-bac-23s-id98.fasta`, `silva-euk-18s-id95.fasta`, `silva-euk-28s-id98.fasta`,
distributed with SortMeRNA) in `${BASE}/db/sortmerna/` (or set `DB`). Both scripts
carry SLURM headers (`sbatch`) and also run with `bash`; finished samples are
skipped. Environment: `envs/preprocessing.yml` (`ENV_NAME`, default `rnaseq_sa_ca`).

The clean files (`*.mRNA_fwd.fq.gz`, `*.mRNA_rev.fq.gz`) are the R1/R2 inputs of
`pipeline/rnaseq_samples.tsv`.
