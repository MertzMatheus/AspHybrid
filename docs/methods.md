# Methods

Detailed methods of the AspHybrid analyses, as described in the
article. Items in [brackets] are filled in from the article and from
`09_tables/software_versions.tsv` (Supplementary Table S19).

## Experimental design and RNA sequencing

Six strains of three *Aspergillus* section *Fumigati* species were studied:
*A. fumigatus* A1160 (reference), the haploid *A. spinulosporus* strains 3 and
5845, and the allodiploid hybrid *A. latus* strains 1, 1710 and 3622
(Steenwyk et al., 2020). Each strain was grown alone in DMEM (fungus-only
control) and in co-culture with murine macrophages [cell type], in three
independent biological replicates per condition (6 strains × 2 conditions × 3
replicates = 36 fungal libraries).

| Species | Ploidy | Strains | Fungal libraries |
|---|---|---|---|
| *A. fumigatus* | haploid | A1160 | 6 |
| *A. spinulosporus* | haploid | 3, 5845 | 12 |
| *A. latus* | allodiploid | 1, 1710, 3622 | 18 |

[Culture and infection conditions: spore preparation, medium, temperature and
CO₂, multiplicity of infection, co-culture time, macrophage-only controls.]
Total RNA was extracted [method and kit] and its integrity checked [instrument,
RIN threshold]. [Strand-specific / unstranded] poly(A)-enriched libraries were
prepared [kit] and sequenced on an Illumina [platform] as [2 × 150 bp]
paired-end reads [provider]. The fungal reads of each library are analysed here; the host (macrophage)
transcriptome is analysed separately.

## Read pre-processing

Raw reads were trimmed and filtered with fastp [version] (Chen et al., 2018)
with automatic adapter detection for paired-end data
(`--detect_adapter_for_pe`). A base was considered qualified at Phred ≥ 20
(`--qualified_quality_phred 20`); reads with more than 40 % unqualified bases,
more than 5 undetermined bases (fastp defaults) or shorter than 30 nt after
trimming (`--length_required 30`) were discarded, together with their mate.

Ribosomal RNA reads were removed with SortMeRNA [version] (Kopylova et al.,
2012) against SILVA-derived databases (Quast et al., 2013) of bacterial 16S
(90 % identity clustering) and 23S (98 %) and eukaryotic 18S (95 %) and 28S
(98 %) rRNA. Pairs were kept together (`--paired_out`): only pairs in which
neither mate matched an rRNA sequence were retained.

## Genome assemblies and quality assessment

Four assemblies were analysed: *A. fumigatus* Af293 (reference assembly, 8
chromosomes, 29.38 Mb; Nierman et al., 2005), *A. spinulosporus* NRRL 2395,
*A. latus* NRRL 200 and *A. quadrilineatus* NRRL 201, the last used only as the
reference for the second parental subgenome of *A. latus*. Each assembly was
converted to a clean copy (alphanumeric sequence identifiers, upper-case
sequence, original headers kept in a correspondence table), used by all later
steps so that genome, annotation and alignments share coordinates. Assembly
statistics were computed with a custom script, and gene-space completeness with
BUSCO v6.1.0 (Manni et al., 2021) in genome mode against eurotiales_odb12.2,
offline.

## Repeat identification and masking

Repeats were identified de novo with RepeatModeler v2.0.9 (Flynn et al., 2020),
including the LTR structural module (`-LTRStruct`), and each genome was
soft-masked with its own library using RepeatMasker v4.2.6 (Smit et al.,
2013–2015; `-xsmall`), both from the Dfam TE Tools container (Storer et al.,
2021) pinned by digest. Soft-masked genomes were used only for gene prediction.

## Alignment of fungal reads

Clean reads were aligned to the genome of their own species with HISAT2
[version] (Kim et al., 2019) with `--dta` and a maximum intron length of
5,000 bp, and alignments were sorted, indexed and checked with SAMtools
[version] (Danecek et al., 2021). fastp was re-run in report-only mode to record
the input of each alignment. *A. latus* reads were aligned to the complete
hybrid genome containing both subgenomes, so that each read could be placed on
its own homeolog. *A. fumigatus* A1160 reads were aligned to the Af293 reference
genome; genes specific to the A1160 lineage are therefore not represented.

## Gene annotation

**A. fumigatus.** The official NCBI RefSeq annotation of Af293
(GCF_000002655.1, ASM265v1; O'Leary et al., 2016) was used with its locus_tags
(AFUA_1G00100 …). Every RefSeq sequence carrying genes was confirmed to be
identical to a sequence of the project genome (MD5 of the upper-case sequence),
so the annotation was placed without liftover. Only protein-coding genes were
kept; isoforms were named `<locus_tag>-T1`, `-T2` … by decreasing protein
length, keeping RefSeq accessions and product names as attributes.

**A. spinulosporus and A. latus.** Genes were predicted with BRAKER v3.0.8 in
ETP mode (Gabriel et al., 2024), which integrates GeneMark-ETP (Brůna et al.,
2024) and AUGUSTUS (Stanke et al., 2008), using the soft-masked genome, the
fungal branch-point model (`--fungus`) and two types of evidence:

- proteins: the OrthoDB v12 Fungi partition (Tegenfeldt et al., 2025),
  deduplicated by sequence (24,091,311 sequences);
- RNA-seq: alignments of one biological replicate of every strain × condition
  group, merged (4 libraries for *A. spinulosporus*, 6 for *A. latus*); only
  mapped primary alignments were kept, and for *A. latus* each library was
  randomly subsampled to at most 15 million read pairs (pair-consistent, fixed
  seed).

During protein-hint generation (ProtHint, Brůna et al., 2021; Spaln, Iwata &
Gotoh, 2012), any single alignment running for more than 2 h was stopped, losing
only the hints of that protein (one alignment in *A. spinulosporus*).

**Curation and locus_tags.** Models with internal stop codons, proteins shorter
than 50 aa, models whose CDS was more than 50 % soft-masked and duplicated
isoforms with identical CDS were removed; partial models were kept and flagged.
Genes received locus_tags in genome order in steps of 10 (ASPN_000010,
ALAT_000010). The longest isoform (-T1) represents each gene downstream.

**Reference files and validation.** GFF3, GTF (gene_id = locus_tag), protein,
CDS and transcript FASTA files were produced and validated: features within the
genome, identical gene sets across files, no internal stop codons, and ≥ 98 % of
proteins identical after independent re-translation of the GTF with gffread
[version] (Pertea & Pertea, 2020). Proteome completeness was assessed with BUSCO
in protein mode.

**Benchmark.** BRAKER3 was also run on Af293 with the same settings, and its
models were matched to the official AFUA genes by CDS overlap on the same strand
(same gene: ≥ 50 % reciprocal overlap) to estimate the accuracy of the de novo
annotation.

## Subgenome assignment in A. latus

The *A. latus* assembly was split into non-overlapping 5-kb windows (windows
< 1 kb discarded), and each window was aligned to the *A. spinulosporus* and
*A. quadrilineatus* genomes with minimap2 [version] (Li, 2018; `-x asm20`, no
secondary alignments). Per parent, the best alignment covering ≥ 50 % of the
window was kept; the window was assigned to S or Q when its identity to one
parent exceeded that to the other by ≥ 2 percentage points, and otherwise
labelled ambiguous or unassigned. Single ambiguous or unassigned windows flanked
by the same subgenome took that subgenome. Windows were merged into blocks,
every gene received the label of its block (`subgenome=S|Q|ambiguous|unassigned`),
and each subgenome was assessed with BUSCO.

## Orthology and gene sharing

Orthologues were inferred with OrthoFinder v3.1.5 (Emms & Kelly, 2019), using
DIAMOND (Buchfink et al., 2021) and gene-tree-based orthologue inference, from
one protein per gene of *A. fumigatus* (AFUA genes), *A. spinulosporus* and
*A. latus*; *A. latus* was one proteome containing both subgenomes. For every
AFUA gene, orthologues were classified as 1:1, 1:2 (one gene per subgenome, the
pattern expected in the hybrid), 1:>2, many:1, many:many or none, and each
*A. latus* orthologue was tagged S or Q. Genes were considered shared with
another species when their orthogroup contained a gene of that species.

## Read counting and expression matrices

Read pairs were assigned to genes with featureCounts (Subread [version]; Liao
et al., 2014), counting fragments (`-p --countReadPairs`), requiring both mates
aligned (`-B`) to the same sequence (`-C`), and using primary alignments with
MAPQ ≥ 10 (`--primary -Q 10`), summarised per gene from exons. Strandedness was
detected automatically from the first library. The MAPQ filter excludes reads
that align equally well to both *A. latus* homeologs. Raw count, TPM and CPM
matrices (rows = locus_tags) were produced per species.

## Functional annotation

The representative protein of every gene was annotated with eggNOG-mapper
[version] (Cantalapiedra et al., 2021) against a local eggNOG 5 database
(Huerta-Cepas et al., 2019), DIAMOND mode, automatic taxonomic scope and
non-electronic GO evidence, with identical settings for all species. GO terms
were mapped to primary identifiers, obsolete terms removed, and each gene
assigned all ancestors of its terms (is_a and part_of; go-basic.obo [release];
Gene Ontology Consortium, 2023); root terms were excluded. KEGG pathway names
were taken from the KEGG REST service (Kanehisa et al., 2023).

## Differential expression and functional enrichment

Differential expression was analysed with DESeq2 [version] (Love et al., 2014)
in R [version], separately for each species. Genes with ≥ 10 reads in ≥ 3
samples were retained. One model per species included all its samples
(`~ group`, group = strain × treatment), and for every strain the macrophage
co-culture was compared with the same strain alone in DMEM (Wald test).
P-values were adjusted with the Benjamini–Hochberg method (Benjamini &
Hochberg, 1995) after independent filtering (Bourgon et al., 2010), and log2
fold changes were shrunk with apeglm [version] (Zhu et al., 2019). Genes with
adjusted p < 0.05 and absolute shrunken log2 fold change ≥ 1 were called DEGs.
Variance-stabilised counts (`vst`, blind = FALSE) of the 500 most variable genes
were used for PCA.

Enrichment of up- and down-regulated DEGs was tested per strain by
over-representation analysis against GO BP, MF and CC, KEGG pathways and COG
categories. For each strain × direction × collection, the universe was the set
of tested genes annotated in that collection, and only annotated DEGs formed the
test set. Gene sets with ≥ 5 genes in the universe (GO sets also ≤ 500) were
tested with a one-sided hypergeometric test, P(X ≥ k), with BH correction within
each strain × direction × collection; adjusted p < 0.05 was considered enriched.
Fold enrichment = (k/n)/(K/N).

## Statistics and reproducibility

All analyses used three independent biological replicates per condition. Wald
tests are two-sided; enrichment tests are one-sided hypergeometric tests. All
P-values were adjusted with the Benjamini–Hochberg false discovery rate.
Environments were managed with micromamba, BRAKER3 and RepeatModeler2 ran in
Docker images pinned by digest, and all tool versions are exported to
Supplementary Table S19. Figures were produced with matplotlib (Hunter, 2007)
and ggplot2.

## References

- Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B* 57:289–300.
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. *PNAS* 107:9546–9551.
- Brůna T, Hoff KJ, Lomsadze A, Stanke M, Borodovsky M. 2021. BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database. *NAR Genom Bioinform* 3:lqaa108.
- Brůna T, Lomsadze A, Borodovsky M. 2024. GeneMark-ETP significantly improves the accuracy of automatic annotation of large eukaryotic genomes. *Genome Res* 34:757–768.
- Buchfink B, Reuter K, Drost HG. 2021. Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat Methods* 18:366–368.
- Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J. 2021. eggNOG-mapper v2. *Mol Biol Evol* 38:5825–5829.
- Chen S, Zhou Y, Chen Y, Gu J. 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34:i884–i890.
- Danecek P, Bonfield JK, Liddle J, et al. 2021. Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008.
- Emms DM, Kelly S. 2019. OrthoFinder: phylogenetic orthology inference for comparative genomics. *Genome Biol* 20:238.
- Flynn JM, Hubley R, Goubert C, et al. 2020. RepeatModeler2 for automated genomic discovery of transposable element families. *PNAS* 117:9451–9457.
- Gabriel L, Brůna T, Hoff KJ, et al. 2024. BRAKER3: fully automated genome annotation using RNA-seq and protein evidence with GeneMark-ETP, AUGUSTUS, and TSEBRA. *Genome Res* 34:769–777.
- Gene Ontology Consortium. 2023. The Gene Ontology knowledgebase in 2023. *Genetics* 224:iyad031.
- Huerta-Cepas J, Szklarczyk D, Heller D, et al. 2019. eggNOG 5.0. *Nucleic Acids Res* 47:D309–D314.
- Hunter JD. 2007. Matplotlib: a 2D graphics environment. *Comput Sci Eng* 9:90–95.
- Iwata H, Gotoh O. 2012. Benchmarking spliced alignment programs including Spaln2. *Nucleic Acids Res* 40:e161.
- Kanehisa M, Furumichi M, Sato Y, Kawashima M, Ishiguro-Watanabe M. 2023. KEGG for taxonomy-based analysis of pathways and genomes. *Nucleic Acids Res* 51:D587–D592.
- Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. 2019. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nat Biotechnol* 37:907–915.
- Kopylova E, Noé L, Touzet H. 2012. SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. *Bioinformatics* 28:3211–3217.
- Li H. 2018. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34:3094–3100.
- Liao Y, Smyth GK, Shi W. 2014. featureCounts. *Bioinformatics* 30:923–930.
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550.
- Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. 2021. BUSCO update. *Mol Biol Evol* 38:4647–4654.
- Nierman WC, Pain A, Anderson MJ, et al. 2005. Genomic sequence of the pathogenic and allergenic filamentous fungus *Aspergillus fumigatus*. *Nature* 438:1151–1156.
- O'Leary NA, Wright MW, Brister JR, et al. 2016. Reference sequence (RefSeq) database at NCBI. *Nucleic Acids Res* 44:D733–D745.
- Pertea G, Pertea M. 2020. GFF utilities: GffRead and GffCompare. *F1000Research* 9:304.
- Quast C, Pruesse E, Yilmaz P, et al. 2013. The SILVA ribosomal RNA gene database project. *Nucleic Acids Res* 41:D590–D596.
- Smit AFA, Hubley R, Green P. 2013–2015. RepeatMasker Open-4.0. http://www.repeatmasker.org.
- Stanke M, Diekhans M, Baertsch R, Haussler D. 2008. Using native and syntenically mapped cDNA alignments to improve de novo gene finding. *Bioinformatics* 24:637–644.
- Steenwyk JL, Lind AL, Ries LNA, et al. 2020. Pathogenic allodiploid hybrids of *Aspergillus* fungi. *Curr Biol* 30:2495–2507.
- Storer J, Hubley R, Rosen J, Wheeler TJ, Smit AF. 2021. The Dfam community resource of transposable element families. *Mob DNA* 12:2.
- Tegenfeldt F, Kuznetsov D, Manni M, et al. 2025. OrthoDB and BUSCO update. *Nucleic Acids Res* 53:D516–D522.
- Zhu A, Ibrahim JG, Love MI. 2019. Heavy-tailed prior distributions for sequence count data. *Bioinformatics* 35:2084–2092.
