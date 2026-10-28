# Installation

## 1. System

| Resource | Used here | Minimum advised |
|---|---|---|
| OS | Ubuntu Linux | any Linux with Bash ≥ 4 |
| CPU threads | 28 | 16 |
| RAM | 123 GB | 64 GB (set `EGGNOG_DBMEM=0` below 64 GB) |
| Free disk | ~1 TB | 500 GB for databases, BAMs and intermediates |
| Software | micromamba, Docker | conda/mamba also work |

Put everything (environments, package cache, temporary files) on the large
data disk. A full `/home` partition stops installations with
"No space left on device"; `pipeline/config.sh` already points `TMPDIR`,
`CONDA_PKGS_DIRS` and `XDG_CACHE_HOME` to `BASE`.

## 2. Get the code

```bash
cd /path/to/project            # e.g. /mnt/bigdata/users/<you>/RNAseq_Diploid
git clone https://github.com/MertzMatheus/AspHybrid.git
cd AspHybrid
chmod +x pipeline/*.sh pipeline/bin/* preprocessing/*.sh
```

## 3. Environments

```bash
export CONDA_PKGS_DIRS=/path/to/project/software/pkgs
micromamba create -y -f envs/annottools.yml  -p /path/to/project/software/envs/AnnotTools
micromamba create -y -f envs/detools.yml     -p /path/to/project/software/envs/DEtools
micromamba create -y -f envs/orthofinder.yml          # named environment "Hibrid"
micromamba create -y -f envs/preprocessing.yml        # named environment "rnaseq_sa_ca"
```

Set the same paths in `pipeline/config.sh` (`SOFT`, `ENV_TOOLS`, `ENV_DE`,
`ENV_ORTHO`). Details per environment: [envs/README.md](../envs/README.md).

## 4. Containers

```bash
docker pull teambraker/braker3@sha256:5f8b3c508a9fe1bbc2e9a74dcc013eeed82f91dd5945adca7823514d9c8aecf8
docker pull dfam/tetools@sha256:7fe6afce694d0759b56d65b623a9d684269e8f82667f0f2eb79fca6a63c51189
```

The digests pin BRAKER 3.0.8 and RepeatModeler 2.0.9 / RepeatMasker 4.2.6.
Your user must be allowed to run Docker (`docker run hello-world`).

## 5. Databases

Expected layout (path `DB` in `config.sh`):

```
db/
├── busco/lineages/eurotiales_odb12.2.tar.gz   (unpacked automatically on first use)
├── orthodb/Fungi.fa.gz                         OrthoDB v12 Fungi proteins
├── eggnog/                                     eggnog.db, eggnog_proteins.dmnd, eggnog.taxa.db ...
├── go/go-basic.obo                             downloaded automatically by step 10
├── kegg/pathway_names.tsv                      downloaded automatically by step 10
└── Af293_RefSeq/                               downloaded automatically by step 06
```

| Database | How to get it | Size |
|---|---|---|
| BUSCO `eurotiales_odb12.2` | `https://busco-data.ezlab.org/v5/data/lineages/` (file `eurotiales_odb12.2*.tar.gz`) | ~0.5 GB |
| OrthoDB v12 Fungi | `wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Fungi.fa.gz` | ~3 GB |
| eggNOG 5 | in the AnnotTools env: `download_eggnog_data.py --data_dir db/eggnog -y` | ~50 GB |

All BUSCO runs are offline; nothing else is downloaded during a run except the
three automatic files above (once).

## 6. Input data

```
Genomes/
├── A_fumigatus/A_fumigatus.fasta           Af293 assembly (must be identical to RefSeq GCF_000002655.1)
├── A_spinulosporus/A_spinulosporus.fasta
├── A_latus/A_latus.fasta
└── A_quadrilineatus/A_quadrilineatus.fasta
```

One FASTA per folder. Raw reads are processed with `preprocessing/` first;
the clean read paths go into `pipeline/rnaseq_samples.tsv`
(see [config/README.md](../config/README.md)).

## 7. Check

```bash
cd pipeline
bash 01_genome_qc.sh A_fumigatus     # minutes; confirms paths, environment and BUSCO lineage
```
