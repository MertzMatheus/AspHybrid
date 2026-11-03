#!/usr/bin/env Rscript
# =============================================================================
# paper_figures.R  -  every figure panel of the article as a separate file
#
# Story: (1) one standardized annotation of A. fumigatus, A. spinulosporus and
# A. latus, and the genes they share; (2) the transcriptional response of each
# strain to macrophages (PCA, DEGs, enrichment); (3) the responses compared through
# A. fumigatus orthologues.
#
# Input
#   Supplementary_Tables_2.xlsx  sheets S1, S3, S5, S7, S8, S9, S10, S12, S13, S14,
#                                S15, S16 and S18 (full enrichment table)
#   Supplementary_Tables_3.xlsx  sheets Samples, DEGs_<strain>, Enrichment significant
#   DE_DIR (optional)            final_results/11_differential_expression of the pipeline:
#                                <species>/deseq2/vst.tsv, samples_used.tsv, DE_<strain>.tsv
#                                -> PCA and volcano plots (skipped when DE_DIR is NA or absent)
#
# Run
#   Rscript paper_figures.R
#   Rscript paper_figures.R <Supplementary_Tables_2.xlsx> <Supplementary_Tables_3.xlsx> \
#                           [DE_DIR or NA] [output folder]
#
# Output (<output folder>/Fig1 ... Fig4, Supplementary): one file per panel, in the
# formats of FORMATS, plus one CSV per panel (the data as plotted) and
# figure_numbers.txt (every number quoted in the Results).
#
#   Fig1  a genes and locus-tag prefix       b BUSCO (genome, proteome)
#         c functional annotation coverage   d COG categories
#         e orthology relative to A. fumigatus
#         f orthogroups shared (UpSet)       g gene sharing per proteome
#   Fig2  a PCA per species                  b DEGs per strain
#         c volcano per strain               d DEGs shared between strains of a species
#   Fig3  a enrichment per strain (top terms, up and down)
#         b terms recurring across strains
#   Fig4  a DEGs in shared or lineage-specific genes
#         b DEGs of each strain compared with A. fumigatus
#         c similarity of the responses (Jaccard)
#         d DEGs shared between species (UpSet)
#         e conserved and lineage-specific genes (heatmap)
#   Supplementary  S1 read alignment and assignment
#                  S2 A. latus subgenomes: a genome, b DEGs, c homeolog pairs
#
# Needs: R >= 4.1 with readxl, ggplot2, patchwork, scales
# =============================================================================

# ---------------- settings (edit here) ----------------
TABLE2  <- "Supplementary_Tables_2.xlsx"
TABLE3  <- "Supplementary_Tables_3.xlsx"
DE_DIR  <- NA          # e.g. "/path/to/project/final_results/11_differential_expression"
OUTDIR  <- NA          # NA = folder "Paper_figures" next to TABLE3
FORMATS <- c("pdf", "png")          # add "tiff" for 600 dpi LZW TIFF
BASE_SIZE <- 7                      # Nature: 5-7 pt
COL1 <- 89; COL2 <- 183             # single / double column width (mm)
PCA_NTOP  <- 500                    # most variable genes in the PCA (as DESeq2::plotPCA)
N_SPECIFIC <- 12                    # Fig. 4e: genes shown per lineage-specific group
MIN_RECUR  <- 3                     # Fig. 3b: term significant in >= this many strain x direction
LFC_CAP    <- 6                     # colour cap in heatmaps

SPECIES  <- c("A_fumigatus", "A_spinulosporus", "A_latus")
SP_LABEL <- c(A_fumigatus = "A. fumigatus", A_spinulosporus = "A. spinulosporus", A_latus = "A. latus")
SP_PREFIX <- c(A_fumigatus = "AFUA_", A_spinulosporus = "ASPN_", A_latus = "ALAT_")
SP_COL   <- c(A_fumigatus = "#E69F00", A_spinulosporus = "#0072B2", A_latus = "#009E73")
STRAINS  <- data.frame(
  species = c("A_fumigatus", "A_spinulosporus", "A_spinulosporus", "A_latus", "A_latus", "A_latus"),
  strain  = c("A1160", "3", "5845", "1", "1710", "3622"), stringsAsFactors = FALSE)
STRAINS$label <- paste(SP_LABEL[STRAINS$species], STRAINS$strain)
STRAIN_COL <- c(A1160 = "#E69F00", `3` = "#56B4E9", `5845` = "#0072B2",
                `1` = "#7FD3B5", `1710` = "#009E73", `3622` = "#00563F")
DIR_COL  <- c(Up = "#B2182B", Down = "#2166AC")
COND_COL <- c(DMEM = "grey60", Macrophage = "#D55E00")
SUB_COL  <- c(S = "#0072B2", Q = "#CC79A7", Ambiguous = "#999999", Unassigned = "#DDDDDD")

# ---------------- command line ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) TABLE2 <- args[1]
if (length(args) >= 2) TABLE3 <- args[2]
if (length(args) >= 3 && args[3] != "NA") DE_DIR <- args[3]
if (length(args) >= 4) OUTDIR <- args[4]
if (is.na(OUTDIR)) OUTDIR <- file.path(dirname(TABLE3), "Paper_figures")

# ---------------- packages ----------------
need <- c("readxl", "ggplot2", "patchwork", "scales")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss))
  stop("missing R packages: ", paste(miss, collapse = ", "),
       "\n  install with: install.packages(c(", paste0('"', miss, '"', collapse = ", "), "))",
       "\n  or: micromamba install -c conda-forge ", paste0("r-", miss, collapse = " "), call. = FALSE)
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
for (f in c(TABLE2, TABLE3)) if (!file.exists(f)) stop("file not found: ", f, call. = FALSE)
HAVE_DE <- !is.na(DE_DIR) && dir.exists(DE_DIR)
for (d in c("Fig1", "Fig2", "Fig3", "Fig4", "Supplementary"))
  dir.create(file.path(OUTDIR, d), showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Helpers
# =============================================================================
num <- function(x) {                      # numbers, "2,826", "5,83E-40"
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x)); x[x %in% c("", "NA", "None", "-")] <- NA
  as.numeric(gsub(",", ".", x, fixed = TRUE))
}
sheet <- function(file, pattern, required = TRUE) {
  sh <- readxl::excel_sheets(file); hit <- sh[grepl(pattern, sh)]
  if (!length(hit)) { if (required) stop("no sheet '", pattern, "' in ", file, call. = FALSE); return(NULL) }
  as.data.frame(readxl::read_excel(file, sheet = hit[1], guess_max = 1e5))
}
split_ids <- function(x) {
  x <- as.character(x); x[is.na(x) | x %in% c("None", "NA", "")] <- ""
  lapply(strsplit(x, ",", fixed = TRUE), trimws)
}
it <- function(x) parse(text = paste0("italic('", x, "')"))
it_strain <- function(lab) {
  sp <- sub(" [^ ]+$", "", lab); st <- sub("^.* ", "", lab)
  parse(text = sprintf("italic('%s')~'%s'", sp, st))
}
theme_paper <- function() {
  theme_classic(base_size = BASE_SIZE) +
    theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          strip.background = element_blank(), strip.text = element_text(face = "bold"),
          legend.key.size = unit(3, "mm"), legend.title = element_text(size = BASE_SIZE),
          legend.text = element_text(hjust = 0),
          panel.grid = element_blank(), plot.title = element_text(size = BASE_SIZE + 1, face = "bold"))
}
save_panel <- function(p, fig, name, w, h) {
  base <- file.path(OUTDIR, fig, name)
  for (fmt in FORMATS) {
    f <- paste0(base, ".", fmt)
    if (fmt == "pdf") ggsave(f, p, width = w, height = h, units = "mm",
                             device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf)
    if (fmt == "png") ggsave(f, p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
    if (fmt == "tiff") ggsave(f, p, width = w, height = h, units = "mm", dpi = 600, compression = "lzw", bg = "white")
  }
  message(sprintf("  %-14s %s", fig, name))
}
write_tab <- function(d, fig, name) write.csv(d, file.path(OUTDIR, fig, paste0(name, ".csv")), row.names = FALSE)
NUMBERS <- character()
note <- function(...) { s <- sprintf(...); NUMBERS <<- c(NUMBERS, s); invisible(s) }
lab_num <- function(v) ifelse(v < 100, scales::comma(v, accuracy = 0.1), scales::comma(v, accuracy = 1))

# UpSet plot: bars of intersection sizes over a dot matrix
upset_plot <- function(d, sets, set_labels, ylab, fill = NULL, fill_values = NULL, italic_sets = TRUE) {
  pats <- unique(d$pattern)
  ord  <- pats[order(-vapply(pats, function(p) sum(d$n[d$pattern == p]), 0))]
  d$pattern <- factor(d$pattern, levels = ord)
  bars <- ggplot(d, aes(pattern, n)) +
    { if (is.null(fill)) geom_col(width = 0.7, fill = "grey30")
      else geom_col(aes(fill = .data[[fill]]), width = 0.75, position = position_dodge(0.8)) } +
    { if (is.null(fill)) geom_text(aes(label = scales::comma(n)), vjust = -0.3, size = BASE_SIZE / 3.2)
      else geom_text(aes(label = scales::comma(n), group = .data[[fill]]), position = position_dodge(0.8),
                     vjust = -0.3, size = BASE_SIZE / 3.6) } +
    { if (!is.null(fill_values)) scale_fill_manual(values = fill_values, name = NULL) } +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::comma) +
    labs(x = NULL, y = ylab) + theme_paper() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          legend.position = "top")
  m <- expand.grid(pattern = ord, set = sets, stringsAsFactors = FALSE)
  m$on <- mapply(function(p, s) s %in% strsplit(p, "+", fixed = TRUE)[[1]], m$pattern, m$set)
  m$pattern <- factor(m$pattern, levels = ord); m$set <- factor(m$set, levels = rev(sets))
  seg <- do.call(rbind, lapply(ord, function(p) {
    on <- m[m$pattern == p & m$on, ]
    if (nrow(on) < 2) return(NULL)
    data.frame(pattern = p, y1 = min(as.integer(on$set)), y2 = max(as.integer(on$set)))
  }))
  if (!is.null(seg)) seg$pattern <- factor(seg$pattern, levels = ord)
  dots <- ggplot(m, aes(pattern, set)) +
    geom_point(aes(colour = on), size = 1.8, show.legend = FALSE) +
    { if (!is.null(seg)) geom_segment(data = seg, aes(x = pattern, xend = pattern, y = y1, yend = y2),
                                      inherit.aes = FALSE, linewidth = 0.5) } +
    scale_colour_manual(values = c(`TRUE` = "grey15", `FALSE` = "grey85")) +
    scale_y_discrete(labels = function(x) if (italic_sets) it_or_strain(set_labels[x]) else set_labels[x]) +
    labs(x = NULL, y = NULL) + theme_paper() +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
  bars / dots + plot_layout(heights = c(3, 1.1))
}
it_or_strain <- function(x) {      # "A. latus" -> italic; "A. latus 1710" -> italic name + strain
  out <- vapply(x, function(s) if (grepl(" [0-9A-Z]+[0-9]+$", s) || grepl(" [0-9]+$", s))
    sprintf("italic('%s')~'%s'", sub(" [^ ]+$", "", s), sub("^.* ", "", s)) else sprintf("italic('%s')", s), "")
  parse(text = out)
}

# =============================================================================
# Data
# =============================================================================
message("Reading tables")
qc    <- sheet(TABLE2, "^S1 ");  qc  <- qc[qc$species %in% SPECIES, ]
ann   <- sheet(TABLE2, "^S5 ");  ann <- ann[ann$species %in% SPECIES, ]
subg  <- sheet(TABLE2, "^S3 ")
ost   <- sheet(TABLE2, "^S7 ")
orth  <- sheet(TABLE2, "^S8 ")
gshr  <- sheet(TABLE2, "^S9 ")
ogshr <- sheet(TABLE2, "^S10 ")
cnt   <- sheet(TABLE2, "^S12 ")
fann  <- sheet(TABLE2, "^S13 ")
cog   <- sheet(TABLE2, "^S14 ")
desm  <- sheet(TABLE2, "^S15 ")
dshr  <- sheet(TABLE2, "^S16 ")
enr_full <- sheet(TABLE2, "^S18 ", required = FALSE)
smp   <- sheet(TABLE3, "^Samples")
enr_cur  <- sheet(TABLE3, "^Enrichment")
if (is.null(enr_full)) enr_full <- enr_cur
degs <- lapply(STRAINS$strain, function(s) sheet(TABLE3, paste0("^DEGs_", s, "$")))
names(degs) <- STRAINS$strain
for (s in names(degs)) {
  degs[[s]]$log2FC <- num(degs[[s]]$log2FC); degs[[s]]$DE <- tolower(trimws(degs[[s]]$DE))
  degs[[s]] <- degs[[s]][degs[[s]]$DE %in% c("up", "down"), ]
}
smp$species_key <- c(`Aspergillus fumigatus` = "A_fumigatus", `Aspergillus spinulosporus` = "A_spinulosporus",
                     `Aspergillus latus` = "A_latus")[trimws(smp$species)]

# orthology universe: A. fumigatus genes with orthologues in both other species
orth$As_rel <- as.character(orth$A_spinulosporus_relation)
orth$Al_rel <- as.character(orth$A_latus_relation)
U3 <- orth$A_fumigatus[!(orth$As_rel %in% c("none", NA)) & !(orth$Al_rel %in% c("none", NA))]
none_both <- orth$A_fumigatus[orth$As_rel %in% "none" & orth$Al_rel %in% "none"]
note("A. fumigatus genes with orthologues in A. spinulosporus and A. latus: %d of %d (%.1f %%)",
     length(U3), nrow(orth), 100 * length(U3) / nrow(orth))

# per strain: A. fumigatus gene -> log2FC (largest |log2FC| when several genes map to it)
afua_fc <- function(st) {
  d <- degs[[st]]
  ids <- if (st == "A1160") as.list(d$locus_tag) else split_ids(d$A_fumigatus_orthologs)
  x <- data.frame(afua = unlist(ids), fc = d$log2FC[rep(seq_len(nrow(d)), lengths(ids))])
  x <- x[x$afua != "", ]; x <- x[order(-abs(x$fc)), ]; x <- x[!duplicated(x$afua), ]
  setNames(x$fc, x$afua)
}
FC <- lapply(STRAINS$strain, afua_fc); names(FC) <- STRAINS$strain
gene_set <- function(st, dirn, universe = U3) {
  f <- FC[[st]]; intersect(names(f)[if (dirn == "up") f > 0 else f < 0], universe)
}

# =============================================================================
# Figure 1 - standardized annotation and gene sharing
# =============================================================================
message("Figure 1")
spx <- function(x) factor(x, levels = SPECIES)

# 1a genes, with the locus-tag prefix of each species
d <- data.frame(species = spx(ann$species), genes = num(ann$genes), transcripts = num(ann$transcripts),
                prefix = SP_PREFIX[ann$species], first = ann$first_locus_tag)
p <- ggplot(d, aes(species, genes, fill = species)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(genes)), vjust = -0.4, size = BASE_SIZE / 3) +
  scale_fill_manual(values = SP_COL) +
  scale_x_discrete(labels = function(x) parse(text = sprintf("atop(italic('%s'), bold('%s'))", SP_LABEL[x],
                                                              sub("_$", "_*", SP_PREFIX[x])))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = scales::comma) +
  labs(x = NULL, y = "Protein-coding genes") + theme_paper()
save_panel(p, "Fig1", "Fig1a_genes_locus_tags", COL1, 65); write_tab(d, "Fig1", "Fig1a_genes_locus_tags")
for (i in seq_len(nrow(d))) note("%s: %s genes, %s transcripts, locus tags %s (first %s)", d$species[i],
                                 d$genes[i], d$transcripts[i], d$prefix[i], d$first[i])

# 1b BUSCO of genome and proteome
busco <- function(x, what) do.call(rbind, lapply(c("S", "D", "F", "M"), function(k)
  data.frame(species = x$species, set = what, class = k, value = num(x[[paste0("BUSCO_", k)]]))))
d <- rbind(busco(qc, "Genome"), busco(ann, "Proteome"))
d$species <- factor(d$species, levels = rev(SPECIES))
d$class <- factor(d$class, levels = c("M", "F", "D", "S"),
                  labels = c("Missing", "Fragmented", "Duplicated", "Single-copy"))
p <- ggplot(d, aes(value, species, fill = class)) +
  geom_col(width = 0.7) + facet_wrap(~set, ncol = 1) +
  scale_fill_manual(values = c(`Single-copy` = "#56B4E9", Duplicated = "#0072B2", Fragmented = "#F0E442",
                               Missing = "#D55E00"), name = NULL, guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(labels = function(x) it(SP_LABEL[x])) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, 100.5)) +
  labs(x = "BUSCO genes (%)", y = NULL) + theme_paper() + theme(legend.position = "bottom")
save_panel(p, "Fig1", "Fig1b_BUSCO", COL1, 65); write_tab(d, "Fig1", "Fig1b_BUSCO")

# 1c functional annotation coverage
keep_ann <- c("eggNOG hit", "named product (not hypothetical)", "Pfam domain", "COG category", "KEGG KO",
              "GO (any)", "KEGG pathway", "EC number", "CAZy")
lab_ann <- c("eggNOG hit", "Named product", "Pfam domain", "COG category", "KEGG orthology",
             "GO term", "KEGG pathway", "EC number", "CAZy family")
f <- fann[fann$annotation %in% keep_ann, ]
d <- do.call(rbind, lapply(SPECIES, function(sp)
  data.frame(species = sp, annotation = f$annotation, pct = num(f[[paste0(sp, "_percent")]]),
             genes = num(f[[paste0(sp, "_genes")]]))))
d$annotation <- factor(lab_ann[match(d$annotation, keep_ann)], levels = rev(lab_ann))
d$species <- factor(d$species, levels = rev(SPECIES))
p <- ggplot(d, aes(pct, annotation, fill = species)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = SP_COL, labels = function(x) it(SP_LABEL[x]), name = NULL,
                    breaks = SPECIES) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, 100)) +
  labs(x = "Proteins annotated (%)", y = NULL) + theme_paper() + theme(legend.position = "bottom")
save_panel(p, "Fig1", "Fig1c_functional_annotation", COL1, 75); write_tab(d, "Fig1", "Fig1c_functional_annotation")

# 1d COG categories
d <- do.call(rbind, lapply(SPECIES, function(sp)
  data.frame(species = sp, COG = cog$COG, description = cog$description, pct = num(cog[[paste0(sp, "_percent")]]))))
d <- d[!is.na(d$pct) & d$pct >= 0.2, ]
d$label <- paste0(d$COG, "  ", ifelse(nchar(d$description) > 42, paste0(substr(d$description, 1, 40), "..."),
                                      d$description))
ordc <- unique(d$label[order(d$COG, decreasing = TRUE)])
d$label <- factor(d$label, levels = ordc); d$species <- factor(d$species, levels = rev(SPECIES))
p <- ggplot(d, aes(pct, label, colour = species)) +
  geom_line(aes(group = label), colour = "grey80", linewidth = 0.4) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = SP_COL, labels = function(x) it(SP_LABEL[x]), name = NULL, breaks = SPECIES) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Proteins in COG category (%)", y = NULL) + theme_paper() +
  theme(legend.position = c(0.78, 0.2), legend.background = element_blank())
save_panel(p, "Fig1", "Fig1d_COG_categories", COL1, 95); write_tab(d, "Fig1", "Fig1d_COG_categories")
qrow <- d[d$COG == "Q", ]
note("COG Q (secondary metabolites): %s", paste(sprintf("%s %.2f %%", qrow$species, qrow$pct), collapse = "; "))

# 1e orthology of A. fumigatus genes
rel <- c("1:1", "1:2", "1:>2", "many:1", "many:many", "none")
rel_lab <- c("1:1", "1:2", "1:>2", "many:1", "many:many", "no orthologue")
d <- ost[ost$comparison %in% c("A_fumigatus -> A_spinulosporus", "A_fumigatus -> A_latus") & ost$category %in% rel, ]
d$target <- factor(ifelse(grepl("spinul", d$comparison), "A_spinulosporus", "A_latus"),
                   levels = c("A_latus", "A_spinulosporus"))
d$category <- factor(rel_lab[match(d$category, rel)], levels = rev(rel_lab))
d$percent <- num(d$percent); d$genes <- num(d$A_fumigatus_genes)
p <- ggplot(d, aes(percent, target, fill = category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(percent >= 8, sprintf("%.0f%%", percent), "")),
            position = position_stack(vjust = 0.5), size = BASE_SIZE / 3.2, colour = "white") +
  scale_fill_manual(values = setNames(c("#1B9E77", "#66A61E", "#A6D854", "#7570B3", "#BEAED4", "#BBBBBB"), rel_lab),
                    name = NULL, guide = guide_legend(reverse = TRUE, nrow = 2)) +
  scale_y_discrete(labels = function(x) it(SP_LABEL[x])) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(paste(italic("A. fumigatus"), " genes (%)")), y = NULL) +
  theme_paper() + theme(legend.position = "bottom")
save_panel(p, "Fig1", "Fig1e_orthology_vs_A_fumigatus", COL1, 50); write_tab(d, "Fig1", "Fig1e_orthology_vs_A_fumigatus")
note("A. fumigatus -> A. spinulosporus 1:1 %s genes (%.1f %%); -> A. latus 1:2 %s genes (%.1f %%)",
     d$genes[d$target == "A_spinulosporus" & d$category == "1:1"], d$percent[d$target == "A_spinulosporus" & d$category == "1:1"],
     d$genes[d$target == "A_latus" & d$category == "1:2"], d$percent[d$target == "A_latus" & d$category == "1:2"])

# 1f orthogroups shared (UpSet)
ab <- c(Af = "A_fumigatus", As = "A_spinulosporus", Al = "A_latus")
d <- data.frame(pattern = vapply(strsplit(ogshr$pattern, "+", fixed = TRUE), function(v) paste(ab[v], collapse = "+"), ""),
                n = num(ogshr$orthogroups))
p <- upset_plot(d, SPECIES, SP_LABEL, "Orthogroups")
save_panel(p, "Fig1", "Fig1f_orthogroups_upset", COL1, 75); write_tab(d, "Fig1", "Fig1f_orthogroups_upset")

# 1g gene sharing per proteome
g <- gshr[gshr$category != "TOTAL", ]
g$class <- ifelse(grepl("^core", g$category), "Core (three species)",
           ifelse(grepl("specific", g$category), "Species-specific", "Shared with one species"))
d <- aggregate(cbind(genes = num(g$genes), pct = num(g$percent_of_genes)) ~ species + class, g, sum)
cls <- c("Core (three species)", "Shared with one species", "Species-specific")
d$class <- factor(d$class, levels = rev(cls)); d$species <- factor(d$species, levels = rev(SPECIES))
p <- ggplot(d, aes(pct, species, fill = class)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(pct >= 6, sprintf("%.0f%%", pct), "")), position = position_stack(vjust = 0.5),
            size = BASE_SIZE / 3.2, colour = "white") +
  scale_fill_manual(values = setNames(c("#332288", "#88CCEE", "#AA4499"), cls), name = NULL,
                    guide = guide_legend(reverse = TRUE, nrow = 2)) +
  scale_y_discrete(labels = function(x) it(SP_LABEL[x])) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Genes (%)", y = NULL) + theme_paper() + theme(legend.position = "bottom")
save_panel(p, "Fig1", "Fig1g_gene_sharing", COL1, 50); write_tab(d, "Fig1", "Fig1g_gene_sharing")
for (i in seq_len(nrow(d))) note("Gene sharing %s - %s: %s genes (%.1f %%)", d$species[i], d$class[i], d$genes[i], d$pct[i])

# =============================================================================
# Figure 2 - RNA-seq: PCA, DEGs, strain overlap
# =============================================================================
message("Figure 2")
# 2a PCA (vst of the pipeline)
if (HAVE_DE) {
  for (sp in SPECIES) {
    vf <- file.path(DE_DIR, sp, "deseq2", "vst.tsv"); mf <- file.path(DE_DIR, sp, "deseq2", "samples_used.tsv")
    if (!file.exists(vf)) { message("  no ", vf, " - PCA skipped"); next }
    v <- as.matrix(read.delim(vf, row.names = 1, check.names = FALSE))
    meta <- if (file.exists(mf)) read.delim(mf, colClasses = "character")
            else smp[smp$species_key == sp, c("sample_old", "strain", "treatment")]
    meta <- meta[match(colnames(v), meta$sample_old), ]
    top <- head(order(apply(v, 1, var), decreasing = TRUE), PCA_NTOP)
    pc <- prcomp(t(v[top, ]))
    ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
    d <- data.frame(sample = colnames(v), PC1 = pc$x[, 1], PC2 = pc$x[, 2], strain = as.character(meta$strain),
                    condition = ifelse(meta$treatment == "DMEN", "DMEM", "Macrophage"))
    p <- ggplot(d, aes(PC1, PC2, colour = strain, shape = condition)) +
      geom_point(size = 2) +
      scale_colour_manual(values = STRAIN_COL, name = "Strain") +
      scale_shape_manual(values = c(DMEM = 1, Macrophage = 16), name = NULL) +
      labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), title = it(SP_LABEL[sp])) +
      theme_paper() + theme(panel.border = element_rect(fill = NA, linewidth = 0.3), axis.line = element_blank())
    save_panel(p, "Fig2", paste0("Fig2a_PCA_", sp), 62, 52); write_tab(d, "Fig2", paste0("Fig2a_PCA_", sp))
    note("PCA %s (top %d genes): PC1 %.1f %%, PC2 %.1f %%", sp, length(top), ve[1], ve[2])
  }
} else message("  DE_DIR not given - PCA (Fig2a) and volcano plots (Fig2c) skipped")

# 2b DEGs per strain
desm$label <- factor(STRAINS$label[match(paste(desm$species, desm$strain), paste(STRAINS$species, STRAINS$strain))],
                     levels = rev(STRAINS$label))
d <- rbind(data.frame(label = desm$label, direction = "Up", n = num(desm$up)),
           data.frame(label = desm$label, direction = "Down", n = -num(desm$down)))
d$direction <- factor(d$direction, levels = c("Up", "Down"))
mx <- max(abs(d$n))
pct <- data.frame(label = desm$label, txt = sprintf("%.1f%%", num(desm$percent_DEGs_of_tested)))
p <- ggplot(d, aes(n, label, fill = direction)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = scales::comma(abs(n)), hjust = ifelse(n > 0, -0.15, 1.15)), size = BASE_SIZE / 3.2) +
  geom_text(data = pct, aes(x = mx * 1.6, y = label, label = txt), inherit.aes = FALSE,
            size = BASE_SIZE / 3.2, colour = "grey30") +
  annotate("text", x = mx * 1.6, y = nrow(desm) + 0.75, label = "% tested", size = BASE_SIZE / 3.2,
           colour = "grey30", fontface = "bold") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = DIR_COL, name = NULL) + scale_y_discrete(labels = it_strain) +
  scale_x_continuous(labels = function(x) scales::comma(abs(x)), limits = c(-mx * 1.35, mx * 1.8)) +
  coord_cartesian(clip = "off") +
  labs(x = "Differentially expressed genes", y = NULL) + theme_paper() + theme(legend.position = "top")
save_panel(p, "Fig2", "Fig2b_DEG_numbers", COL1, 60); write_tab(d, "Fig2", "Fig2b_DEG_numbers")
for (i in seq_len(nrow(desm))) note("DEGs %s %s: %s up, %s down, %s total (%.1f %% of %s tested)", desm$species[i],
     desm$strain[i], desm$up[i], desm$down[i], desm$DEGs[i], num(desm$percent_DEGs_of_tested[i]), desm$genes_tested[i])

# 2c volcano plots (all tested genes from the pipeline)
if (HAVE_DE) {
  for (i in seq_len(nrow(STRAINS))) {
    sp <- STRAINS$species[i]; st <- STRAINS$strain[i]
    f <- file.path(DE_DIR, sp, "deseq2", paste0("DE_", st, ".tsv"))
    if (!file.exists(f)) { message("  no ", f, " - volcano skipped"); next }
    r <- read.delim(f); r <- r[!is.na(r$padj), ]
    r$class <- ifelse(r$padj < 0.05 & r$log2FC >= 1, "Up", ifelse(r$padj < 0.05 & r$log2FC <= -1, "Down", "NS"))
    r$y <- pmin(-log10(pmax(r$padj, 1e-300)), 100)
    lab <- head(r[r$class != "NS", ][order(-abs(r$log2FC[r$class != "NS"]) * r$y[r$class != "NS"]), ], 6)
    p <- ggplot(r, aes(log2FC, y, colour = class)) +
      geom_point(size = 0.35, alpha = 0.7) +
      geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.25, colour = "grey40") +
      geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25, colour = "grey40") +
      geom_text(data = lab, aes(label = locus_tag), size = BASE_SIZE / 3.8, colour = "black",
                hjust = "inward", vjust = "inward", check_overlap = TRUE) +
      scale_colour_manual(values = c(DIR_COL, NS = "grey80"), guide = "none") +
      labs(x = expression(log[2]~"fold change"), y = expression(-log[10]~"adjusted"~italic(P)),
           title = it_strain(STRAINS$label[i])) + theme_paper()
    save_panel(p, "Fig2", paste0("Fig2c_volcano_", st), 60, 55)
  }
}

# 2d DEGs shared between strains of the same species (UpSet)
for (sp in c("A_spinulosporus", "A_latus")) {
  x <- dshr[dshr$species == sp, ]
  sts <- STRAINS$strain[STRAINS$species == sp]
  d <- data.frame(pattern = as.character(x$strains), n = num(x$genes),
                  direction = factor(ifelse(x$direction == "up", "Up", "Down"), levels = c("Up", "Down")))
  p <- upset_plot(d, sts, setNames(paste(SP_LABEL[sp], sts), sts), "DEGs", fill = "direction", fill_values = DIR_COL)
  save_panel(p, "Fig2", paste0("Fig2d_strain_overlap_", sp), COL1, 70)
  write_tab(d, "Fig2", paste0("Fig2d_strain_overlap_", sp))
  allp <- paste(sts, collapse = "+")
  note("%s: DEGs in all strains %d up, %d down; strain-specific up %s; down %s", sp,
       sum(d$n[d$pattern == allp & d$direction == "Up"]), sum(d$n[d$pattern == allp & d$direction == "Down"]),
       paste(sprintf("%s=%d", sts, vapply(sts, function(s) sum(d$n[d$pattern == s & d$direction == "Up"]), 0)), collapse = " "),
       paste(sprintf("%s=%d", sts, vapply(sts, function(s) sum(d$n[d$pattern == s & d$direction == "Down"]), 0)), collapse = " "))
}

# =============================================================================
# Figure 3 - functional enrichment
# =============================================================================
message("Figure 3")
prep_enr <- function(e) {
  e$strain <- as.character(e$strain); e$padj <- num(e$padj); e$fold_enrichment <- num(e$fold_enrichment)
  e$k <- num(e$DEGs_in_term)
  nm <- if ("name" %in% names(e)) e$name else e$Name
  e$name <- ifelse(is.na(nm) | nm == "None", e$term, nm)
  cap <- grepl("^[a-z][a-z]", e$name)                  # keep "rRNA", "tRNA" as they are
  e$name[cap] <- paste0(toupper(substr(e$name[cap], 1, 1)), substr(e$name[cap], 2, nchar(e$name[cap])))
  e$direction <- ifelse(tolower(e$direction) == "up", "Up", "Down")
  e$coll <- c(GO_BP = "BP", GO_MF = "MF", GO_CC = "CC", KEGG_pathway = "KEGG", COG = "COG")[e$collection]
  e
}
enr_cur <- prep_enr(enr_cur); enr_full <- prep_enr(enr_full)
if ("significant" %in% names(enr_full)) enr_full <- enr_full[enr_full$significant %in% c("yes", TRUE) | enr_full$padj < 0.05, ]
wrap <- function(x, n = 45) vapply(strwrap(x, n, simplify = FALSE), paste, "", collapse = "\n")

# 3a one panel per strain: top terms, induced to the right, repressed to the left
for (i in seq_len(nrow(STRAINS))) {
  st <- STRAINS$strain[i]
  d <- enr_cur[enr_cur$strain == st, ]
  if (!nrow(d)) next
  d$x <- ifelse(d$direction == "Up", 1, -1) * -log10(d$padj)
  d$lab <- paste0(d$name, " (", d$coll, ")")
  d <- d[order(d$direction == "Up", ifelse(d$direction == "Up", d$x, -d$x)), ]
  d$lab <- factor(make.unique(d$lab, sep = " "), levels = make.unique(d$lab, sep = " "))
  p <- ggplot(d, aes(x, lab, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = k, hjust = ifelse(x > 0, -0.2, 1.2)), size = BASE_SIZE / 3.4) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = DIR_COL, name = NULL, breaks = c("Up", "Down")) +
    scale_x_continuous(labels = function(v) abs(v), expand = expansion(mult = 0.12)) +
    scale_y_discrete(labels = function(v) wrap(sub(" [0-9]+$", "", v), 50)) +
    labs(x = expression(-log[10]~"adjusted"~italic(P)), y = NULL, title = it_strain(STRAINS$label[i])) +
    theme_paper() + theme(legend.position = "top", axis.text.y = element_text(size = BASE_SIZE - 1))
  save_panel(p, "Fig3", paste0("Fig3a_enrichment_", st), COL1, 20 + 4.2 * nrow(d))
  write_tab(d[, c("strain", "direction", "collection", "term", "name", "k", "fold_enrichment", "padj")],
            "Fig3", paste0("Fig3a_enrichment_", st))
}

# 3b curated terms that recur across strains (values from the full enrichment table)
key <- function(e) paste(e$collection, e$term)
cur_keys <- unique(key(enr_cur))
recur <- table(key(enr_full)[key(enr_full) %in% cur_keys])
sel <- names(recur)[recur >= MIN_RECUR]
d <- enr_full[key(enr_full) %in% sel, ]
d$col <- factor(paste(STRAINS$label[match(d$strain, STRAINS$strain)], d$direction),
                levels = as.vector(t(outer(STRAINS$label, c("Up", "Down"), paste))))
d$y <- paste0(d$name, " (", d$coll, ")")
pattern <- tapply(paste(d$strain, d$direction), d$y, function(v) length(v))
mean_dir <- tapply(ifelse(d$direction == "Up", 1, -1), d$y, mean)
d$y <- factor(d$y, levels = names(sort(mean_dir + pattern / 100)))
d$mlp <- pmin(-log10(d$padj), 20)
p <- ggplot(d, aes(col, y)) +
  geom_point(aes(size = mlp, colour = direction)) +
  geom_vline(xintercept = seq(2.5, 10.5, 2), colour = "grey90", linewidth = 0.3) +
  scale_colour_manual(values = DIR_COL, name = NULL) +
  scale_size_continuous(range = c(0.6, 3.2), name = expression(-log[10]~italic(P)[adj]), limits = c(1.3, 20)) +
  scale_x_discrete(drop = FALSE, labels = function(v) {
    s <- sub(" (Up|Down)$", "", v); dd <- sub("^.* ", "", v)
    parse(text = sprintf("italic('%s')~'%s %s'", sub(" [^ ]+$", "", s), sub("^.* ", "", s), dd)) }) +
  labs(x = NULL, y = NULL) + theme_paper() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1), axis.text.y = element_text(size = BASE_SIZE - 1),
        panel.border = element_rect(fill = NA, linewidth = 0.3), axis.line = element_blank())
save_panel(p, "Fig3", "Fig3b_recurrent_terms", COL2, 30 + 3.4 * length(sel))
write_tab(d[, c("y", "strain", "direction", "collection", "term", "k", "fold_enrichment", "padj")], "Fig3", "Fig3b_recurrent_terms")
note("Recurrent enrichment terms (curated, significant in >= %d strain x direction): %d", MIN_RECUR, length(sel))

# =============================================================================
# Figure 4 - the responses compared through A. fumigatus orthologues
# =============================================================================
message("Figure 4")
cls_lev <- c("Shared by the three species", "Shared with one other species", "Species- or lineage-specific")
d <- do.call(rbind, lapply(seq_len(nrow(STRAINS)), function(i) {
  st <- STRAINS$strain[i]; x <- degs[[st]]
  if (st == "A1160") cls <- ifelse(x$locus_tag %in% U3, cls_lev[1], ifelse(x$locus_tag %in% none_both, cls_lev[3], cls_lev[2]))
  else cls <- vapply(split_ids(x$A_fumigatus_orthologs), function(v) {
    v <- v[v != ""]; if (!length(v)) cls_lev[3] else if (any(v %in% U3)) cls_lev[1] else cls_lev[2] }, "")
  tb <- table(factor(cls, levels = cls_lev))
  data.frame(label = STRAINS$label[i], class = names(tb), n = as.integer(tb))
}))
d$pct <- ave(d$n, d$label, FUN = function(v) 100 * v / sum(v))
d$label <- factor(d$label, levels = rev(STRAINS$label)); d$class <- factor(d$class, levels = rev(cls_lev))
p <- ggplot(d, aes(pct, label, fill = class)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ifelse(pct >= 7, sprintf("%.0f%%", pct), "")), position = position_stack(vjust = 0.5),
            size = BASE_SIZE / 3.2, colour = "white") +
  scale_fill_manual(values = setNames(c("#332288", "#88CCEE", "#AA4499"), cls_lev), name = NULL,
                    guide = guide_legend(reverse = TRUE, ncol = 1)) +
  scale_y_discrete(labels = it_strain) + scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "DEGs (%)", y = NULL) + theme_paper() + theme(legend.position = "bottom")
save_panel(p, "Fig4", "Fig4a_DEGs_shared_or_specific", COL1, 65); write_tab(d, "Fig4", "Fig4a_DEGs_shared_or_specific")
for (l in STRAINS$label) { x <- d[d$label == l, ]
  note("%s: %d DEGs; shared by the three species %.1f %%; species/lineage-specific %.1f %% (%d)", l, sum(x$n),
       x$pct[x$class == cls_lev[1]], x$pct[x$class == cls_lev[3]], x$n[x$class == cls_lev[3]]) }

# 4b each Nidulantes strain against A. fumigatus (genes in the comparison universe)
af_fc <- FC[["A1160"]]
d <- do.call(rbind, lapply(2:nrow(STRAINS), function(i) {
  f <- FC[[STRAINS$strain[i]]]; f <- f[names(f) %in% U3]
  a <- af_fc[names(f)]
  cls <- ifelse(is.na(a), "Not DE in A. fumigatus", ifelse(sign(a) == sign(f), "Same direction in A. fumigatus",
                                                             "Opposite direction in A. fumigatus"))
  tb <- table(factor(cls, levels = c("Same direction in A. fumigatus", "Opposite direction in A. fumigatus",
                                     "Not DE in A. fumigatus")))
  data.frame(label = STRAINS$label[i], class = names(tb), n = as.integer(tb))
}))
d$label <- factor(d$label, levels = rev(STRAINS$label[-1]))
d$class <- factor(d$class, levels = rev(c("Same direction in A. fumigatus", "Opposite direction in A. fumigatus",
                                          "Not DE in A. fumigatus")))
d$pct <- ave(d$n, d$label, FUN = function(v) 100 * v / sum(v))
tot <- aggregate(n ~ label, d, sum)
p <- ggplot(d, aes(pct, label, fill = class)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ifelse(pct >= 6, sprintf("%.0f%%", pct), "")), position = position_stack(vjust = 0.5),
            size = BASE_SIZE / 3.2, colour = "white") +
  geom_text(data = tot, aes(x = 101, y = label, label = paste0("n = ", scales::comma(n))), inherit.aes = FALSE,
            hjust = 0, size = BASE_SIZE / 3.2, colour = "grey30") +
  coord_cartesian(xlim = c(0, 100), clip = "off") +
  scale_fill_manual(values = c(`Same direction in A. fumigatus` = "#117733", `Opposite direction in A. fumigatus` = "#882255",
                               `Not DE in A. fumigatus` = "#BBBBBB"),
                    labels = function(x) parse(text = sprintf("'%s'~italic('A. fumigatus')", sub(" A. fumigatus$", "", x))),
                    name = NULL, guide = guide_legend(reverse = TRUE, ncol = 1)) +
  scale_y_discrete(labels = it_strain) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = expression(paste("DEGs with an ", italic("A. fumigatus"), " orthologue (%)")), y = NULL) +
  theme_paper() + theme(legend.position = "bottom", legend.text = element_text(hjust = 0),
                        plot.margin = margin(5, 28, 5, 5))
save_panel(p, "Fig4", "Fig4b_DEGs_vs_A_fumigatus", COL1, 60); write_tab(d, "Fig4", "Fig4b_DEGs_vs_A_fumigatus")
for (l in levels(d$label)) { x <- d[d$label == l, ]
  note("%s vs A. fumigatus: %d orthologue DEGs; same direction %d; opposite %d; not DE in A. fumigatus %d (%.1f %%)", l,
       sum(x$n), x$n[x$class == "Same direction in A. fumigatus"], x$n[x$class == "Opposite direction in A. fumigatus"],
       x$n[x$class == "Not DE in A. fumigatus"], 100 * x$n[x$class == "Not DE in A. fumigatus"] / sum(x$n)) }

# 4c Jaccard similarity (same direction, comparison universe)
SG <- lapply(STRAINS$strain, function(st) { f <- FC[[st]]; f <- f[names(f) %in% U3]; paste(names(f), sign(f)) })
names(SG) <- STRAINS$label
jac <- expand.grid(a = STRAINS$label, b = STRAINS$label, stringsAsFactors = FALSE)
jac$J <- mapply(function(a, b) length(intersect(SG[[a]], SG[[b]])) / length(union(SG[[a]], SG[[b]])), jac$a, jac$b)
jac$a <- factor(jac$a, levels = STRAINS$label); jac$b <- factor(jac$b, levels = rev(STRAINS$label))
low <- jac[as.integer(jac$a) < (nrow(STRAINS) + 1 - as.integer(jac$b)), ]
p <- ggplot(low, aes(a, b, fill = J)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", J), colour = J > 0.22), size = BASE_SIZE / 3, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08306B", name = "Jaccard", limits = c(0, NA)) +
  scale_x_discrete(labels = it_strain) + scale_y_discrete(labels = it_strain) +
  labs(x = NULL, y = NULL) + coord_equal() + theme_paper() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
save_panel(p, "Fig4", "Fig4c_jaccard", COL1, 80); write_tab(jac, "Fig4", "Fig4c_jaccard")
J <- function(a, b) jac$J[as.character(jac$a) == a & as.character(jac$b) == b]
L <- STRAINS$label
note("Jaccard: A. fumigatus vs others %.2f-%.2f; As 3-As 5845 %.2f; As 5845-Al 3622 %.2f; As 3-Al 1 %.2f; Al 1-Al 1710 %.2f; Al 1710-Al 3622 %.2f",
     min(sapply(L[-1], J, a = L[1])), max(sapply(L[-1], J, a = L[1])), J(L[2], L[3]), J(L[3], L[6]), J(L[2], L[4]),
     J(L[4], L[5]), J(L[5], L[6]))

# 4d DEGs shared between species (>= 1 strain per species)
sp_set <- function(sp, dirn, how) {
  s <- lapply(STRAINS$strain[STRAINS$species == sp], gene_set, dirn = dirn)
  if (how == "any") Reduce(union, s) else Reduce(intersect, s)
}
upset_counts <- function(how) do.call(rbind, lapply(c("up", "down"), function(dirn) {
  S <- lapply(SPECIES, sp_set, dirn = dirn, how = how); names(S) <- SPECIES
  allg <- Reduce(union, S)
  memb <- vapply(allg, function(g) paste(SPECIES[vapply(S, function(s) g %in% s, TRUE)], collapse = "+"), "")
  tb <- table(memb)
  data.frame(pattern = names(tb), n = as.integer(tb), direction = ifelse(dirn == "up", "Up", "Down"), how = how)
}))
up_any <- upset_counts("any"); up_all <- upset_counts("all")
up_any$direction <- factor(up_any$direction, levels = c("Up", "Down"))
p <- upset_plot(up_any, SPECIES, SP_LABEL, expression(paste("DEGs (", italic("A. fumigatus"), " orthologues)")),
                fill = "direction", fill_values = DIR_COL)
save_panel(p, "Fig4", "Fig4d_species_overlap_upset", COL2 * 0.62, 75)
write_tab(rbind(up_any, up_all), "Fig4", "Fig4d_species_overlap_upset")
for (h in list(up_any, up_all)) for (dd in c("Up", "Down")) {
  x <- h[h$direction == dd, ]; s <- function(p) sum(x$n[x$pattern == p])
  note("Shared DEGs (%s strain(s) per species), %s: all three %d; As+Al %d; Af+As %d; Af+Al %d; Af only %d; As only %d; Al only %d",
       ifelse(h$how[1] == "any", ">= 1", "all"), dd, s(paste(SPECIES, collapse = "+")), s("A_spinulosporus+A_latus"),
       s("A_fumigatus+A_spinulosporus"), s("A_fumigatus+A_latus"), s("A_fumigatus"), s("A_spinulosporus"), s("A_latus"))
}

# 4e conserved and lineage-specific genes
af <- degs[["A1160"]]
short_name <- function(prod, pref) {
  p <- ifelse(is.na(pref) | pref %in% c("None", "-", ""), "", as.character(pref))
  d <- ifelse(is.na(prod) | prod %in% c("None", ""), "uncharacterized protein", as.character(prod))
  d <- sub("^(putative|probable) ", "", d)
  d <- ifelse(nchar(d) > 38, paste0(substr(d, 1, 36), "..."), d)
  has <- mapply(function(a, b) a != "" && !grepl(tolower(a), tolower(b), fixed = TRUE), p, d)
  ifelse(has, paste0(d, " (", p, ")"), d)
}
af_name <- setNames(short_name(af$product, af$eggNOG_preferred_name), af$locus_tag)
core <- lapply(c("up", "down"), function(dirn) Reduce(intersect, lapply(STRAINS$strain, gene_set, dirn = dirn)))
af_spec <- af[af$locus_tag %in% none_both, ]
af_spec <- head(af_spec[order(-abs(af_spec$log2FC)), ], N_SPECIFIC)
al_fc <- lapply(c("1", "1710", "3622"), function(st) setNames(degs[[st]]$log2FC, degs[[st]]$locus_tag))
names(al_fc) <- c("1", "1710", "3622")
as_fc <- lapply(c("3", "5845"), function(st) setNames(degs[[st]]$log2FC, degs[[st]]$locus_tag)); names(as_fc) <- c("3", "5845")
nid <- unique(do.call(rbind, lapply(c("3", "5845"), function(st) {
  x <- degs[[st]]; x[vapply(split_ids(x$A_fumigatus_orthologs), function(v) all(v == ""), TRUE),
                     c("locus_tag", "product", "eggNOG_preferred_name", "A_latus_orthologs")] })))
nid <- nid[!duplicated(nid$locus_tag), ]
nid_fc <- t(vapply(seq_len(nrow(nid)), function(i) {
  g <- nid$locus_tag[i]; alo <- split_ids(nid$A_latus_orthologs[i])[[1]]
  c(A1160 = NA, `3` = unname(as_fc[["3"]][g]), `5845` = unname(as_fc[["5845"]][g]),
    vapply(c("1", "1710", "3622"), function(st) { v <- al_fc[[st]][intersect(alo, names(al_fc[[st]]))]
      if (length(v)) v[which.max(abs(v))] else NA_real_ }, 0))
}, numeric(6)))
colnames(nid_fc) <- STRAINS$strain
ok <- apply(nid_fc[, -1, drop = FALSE], 1, function(v) { v <- v[!is.na(v)]; length(v) >= 4 && (all(v > 0) || all(v < 0)) })
note("Nidulantes-specific DEGs (no A. fumigatus orthologue; same direction in >= 4 of 5 Nidulantes strains): %d", sum(ok))
nid_sel <- which(ok)
nid_name <- short_name(nid$product, nid$eggNOG_preferred_name)
named <- !grepl("^(uncharacterized|hypothetical) protein$", nid_name)
nid_sel <- nid_sel[order(!named[nid_sel], -rowMeans(abs(nid_fc[nid_sel, -1, drop = FALSE]), na.rm = TRUE))]
nid_sel <- head(nid_sel, N_SPECIFIC)
rows <- function(id, nm, grp, fcs) data.frame(id = id, gene = paste0(nm, "  ", id), group = grp,
                                             label = STRAINS$label, fc = fcs)
d <- rbind(
  do.call(rbind, lapply(unlist(core), function(g) rows(g, af_name[g], "core", vapply(STRAINS$strain, function(st) unname(FC[[st]][g]), 0)))),
  do.call(rbind, lapply(af_spec$locus_tag, function(g) rows(g, af_name[g], "af", c(af_fc[g], rep(NA, 5))))),
  do.call(rbind, lapply(nid_sel, function(i) rows(nid$locus_tag[i], nid_name[i], "nid", nid_fc[i, ]))))
d$absent <- (d$group == "af" & d$label != L[1]) | (d$group == "nid" & d$label == L[1])
d$group <- factor(d$group, levels = c("core", "af", "nid"),
                  labels = c("'Conserved in all strains'", "italic('A. fumigatus')~'only'", "italic('Nidulantes')~'only'"))
ordg <- unique(d$gene[order(d$group, -ave(d$fc, d$id, FUN = function(v) mean(v, na.rm = TRUE)))])
d$gene <- factor(d$gene, levels = rev(ordg)); d$label <- factor(d$label, levels = L)
d$fc_plot <- pmax(pmin(d$fc, LFC_CAP), -LFC_CAP)
p <- ggplot(d, aes(label, gene)) +
  geom_tile(aes(fill = fc_plot), colour = "white", linewidth = 0.3) +
  geom_point(data = d[d$absent, ], shape = 4, size = 0.9, colour = "grey45") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", labeller = label_parsed) +
  scale_fill_gradient2(low = DIR_COL[["Down"]], mid = "white", high = DIR_COL[["Up"]], midpoint = 0,
                       limits = c(-LFC_CAP, LFC_CAP), na.value = "grey92", name = expression(log[2]~"FC")) +
  scale_x_discrete(labels = it_strain, position = "top") +
  labs(x = NULL, y = NULL) + theme_paper() +
  theme(axis.text.x.top = element_text(angle = 40, hjust = 0), axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0), axis.text.y = element_text(size = BASE_SIZE - 1.5))
save_panel(p, "Fig4", "Fig4e_conserved_and_specific_genes", COL2 * 0.75, 40 + 3.1 * length(unique(d$gene)))
write_tab(d, "Fig4", "Fig4e_conserved_and_specific_genes")
note("Conserved core (same direction in all six strains): %d up, %d down", length(core[[1]]), length(core[[2]]))
note("A. fumigatus DEGs without orthologue in either other species: %d (%d up, %d down)",
     sum(af$locus_tag %in% none_both), sum(af$DE == "up" & af$locus_tag %in% none_both),
     sum(af$DE == "down" & af$locus_tag %in% none_both))

# =============================================================================
# Supplementary - mapping and A. latus subgenomes
# =============================================================================
message("Supplementary")
smp$key <- paste(smp$species_key, smp$sample_old); cnt$key <- paste(cnt$species, cnt$sample)
cnt$strain <- as.character(smp$strain[match(cnt$key, smp$key)])
cnt$condition <- ifelse(tolower(cnt$condition) == "control", "DMEM", "Macrophage")
cnt$label <- factor(L[match(paste(cnt$species, cnt$strain), paste(STRAINS$species, STRAINS$strain))], levels = L)
d <- rbind(data.frame(cnt[, c("sample", "label", "condition")], metric = "Aligned to genome", value = num(cnt$overall_alignment_pct)),
           data.frame(cnt[, c("sample", "label", "condition")], metric = "Assigned to genes", value = num(cnt$pct_assigned_of_clean)))
p <- ggplot(d, aes(label, value, colour = condition)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1), size = 1.1) +
  geom_text(data = d[d$value < 80 & d$metric == "Aligned to genome", ], aes(label = sample), colour = "black",
            size = BASE_SIZE / 3.2, hjust = -0.25, show.legend = FALSE) +
  facet_wrap(~metric) + scale_colour_manual(values = COND_COL, name = NULL) +
  scale_x_discrete(labels = it_strain) + scale_y_continuous(limits = c(50, 100)) +
  labs(x = NULL, y = "Read pairs (%)") + theme_paper() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
save_panel(p, "Supplementary", "FigS1_read_mapping", COL2 * 0.7, 70); write_tab(d, "Supplementary", "FigS1_read_mapping")

sv <- setNames(num(subg$value), subg$metric); al_ann <- ann[ann$species == "A_latus", ]
d <- rbind(data.frame(level = "Sequence (Mb)", class = c("S", "Q", "Ambiguous", "Unassigned"),
                      value = sv[c("S_Mb", "Q_Mb", "ambiguous_Mb", "unassigned_Mb")]),
           data.frame(level = "Genes", class = c("S", "Q", "Ambiguous", "Unassigned"),
                      value = num(unlist(al_ann[c("subgenome_S", "subgenome_Q", "subgenome_ambiguous", "subgenome_unassigned")]))))
d$pct <- ave(d$value, d$level, FUN = function(v) 100 * v / sum(v))
d$class <- factor(d$class, levels = rev(names(SUB_COL))); d$level <- factor(d$level, levels = c("Genes", "Sequence (Mb)"))
p <- ggplot(d, aes(pct, level, fill = class)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(class %in% c("S", "Q"), paste0(class, "\n", lab_num(value)), "")),
            position = position_stack(vjust = 0.5), size = BASE_SIZE / 3.2, colour = "white", lineheight = 0.85) +
  scale_fill_manual(values = SUB_COL, name = NULL, guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(paste("Share of the ", italic("A. latus"), " genome (%)")), y = NULL) +
  theme_paper() + theme(legend.position = "bottom")
save_panel(p, "Supplementary", "FigS2a_A_latus_subgenomes", COL1, 45); write_tab(d, "Supplementary", "FigS2a_A_latus_subgenomes")

bg <- num(unlist(al_ann[c("subgenome_S", "subgenome_Q")])); bg_other <- num(al_ann$genes) - sum(bg)
d <- rbind(data.frame(x = "Genome", sub = c("S", "Q", "Other"), n = c(bg, bg_other)),
           do.call(rbind, lapply(c("1", "1710", "3622"), function(st) {
             x <- degs[[st]]; s <- ifelse(x$subgenome %in% c("S", "Q"), x$subgenome, "Other")
             tb <- as.data.frame(table(paste(st, ifelse(x$DE == "up", "Up", "Down")), s), stringsAsFactors = FALSE)
             setNames(tb, c("x", "sub", "n")) })))
d$pct <- ave(d$n, d$x, FUN = function(v) 100 * v / sum(v))
d$x <- factor(d$x, levels = c("Genome", as.vector(t(outer(c("1", "1710", "3622"), c("Up", "Down"), paste)))))
d$sub <- factor(d$sub, levels = c("Other", "Q", "S"))
p <- ggplot(d, aes(x, pct, fill = sub)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 100 * bg[1] / (sum(bg) + bg_other), linetype = 2, linewidth = 0.3) +
  geom_text(aes(label = ifelse(sub != "Other", sprintf("%.0f", pct), "")), position = position_stack(vjust = 0.5),
            size = BASE_SIZE / 3.2, colour = "white") +
  scale_fill_manual(values = c(S = SUB_COL[["S"]], Q = SUB_COL[["Q"]], Other = "grey80"), name = "Subgenome",
                    breaks = c("S", "Q", "Other")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression(paste(italic("A. latus"), " strain and direction")), y = "Genes (%)") +
  theme_paper() + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_panel(p, "Supplementary", "FigS2b_subgenome_of_DEGs", COL1, 60); write_tab(d, "Supplementary", "FigS2b_subgenome_of_DEGs")

pairs <- orth[orth$Al_rel %in% "1:2", c("A_fumigatus", "A_latus_S", "A_latus_Q")]
pairs <- pairs[!is.na(pairs$A_latus_S) & !is.na(pairs$A_latus_Q) & !grepl(",", pairs$A_latus_S) &
               !grepl(",", pairs$A_latus_Q) & !(pairs$A_latus_S %in% "None") & !(pairs$A_latus_Q %in% "None"), ]
pl <- c("Both, same direction", "S homeolog only", "Q homeolog only", "Opposite directions")
d <- do.call(rbind, lapply(c("1", "1710", "3622"), function(st) {
  f <- al_fc[[st]]; s <- f[pairs$A_latus_S]; q <- f[pairs$A_latus_Q]
  cls <- ifelse(!is.na(s) & !is.na(q), ifelse(sign(s) == sign(q), pl[1], pl[4]),
         ifelse(!is.na(s), pl[2], ifelse(!is.na(q), pl[3], NA)))
  tb <- table(factor(cls, levels = pl))
  data.frame(label = L[STRAINS$strain == st], class = names(tb), n = as.integer(tb))
}))
d$class <- factor(d$class, levels = rev(pl)); d$label <- factor(d$label, levels = L[4:6])
p <- ggplot(d, aes(label, n, fill = class)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ifelse(n >= 60, scales::comma(n), "")), position = position_stack(vjust = 0.5),
            size = BASE_SIZE / 3.2, colour = "white") +
  scale_fill_manual(values = setNames(c("#117733", SUB_COL[["S"]], SUB_COL[["Q"]], "black"), pl), name = NULL,
                    guide = guide_legend(reverse = TRUE)) +
  scale_x_discrete(labels = it_strain) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::comma) +
  labs(x = NULL, y = "S-Q homeolog pairs with a DE copy") + theme_paper()
save_panel(p, "Supplementary", "FigS2c_homeolog_pairs", COL1, 60); write_tab(d, "Supplementary", "FigS2c_homeolog_pairs")
note("One-to-one S-Q homeolog pairs: %d", nrow(pairs))

writeLines(c(paste("paper_figures.R -", format(Sys.time(), "%Y-%m-%d %H:%M")),
             paste("Table 2:", normalizePath(TABLE2)), paste("Table 3:", normalizePath(TABLE3)),
             paste("DE_DIR :", if (HAVE_DE) normalizePath(DE_DIR) else "not given (no PCA / volcano)"), "", NUMBERS),
           file.path(OUTDIR, "figure_numbers.txt"))
message("Done: ", normalizePath(OUTDIR))
