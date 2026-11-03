#!/usr/bin/env Rscript
# =============================================================================
# paper_figures.R  -  main figures of the article (A. fumigatus, A. spinulosporus,
#                     A. latus; fungal RNA-seq after co-culture with macrophages)
#
# Input : Supplementary_Tables_2.xlsx  genomes, annotation, orthology, mapping
#                                      (sheets S1, S3, S5, S7, S8, S10, S12, S15)
#         Supplementary_Tables_3.xlsx  RNA-seq: Samples, DEGs_<strain>
#
# Run   : Rscript paper_figures.R
#         Rscript paper_figures.R <Supplementary_Tables_2.xlsx> <Supplementary_Tables_3.xlsx> <output folder>
#         or in RStudio: source("paper_figures.R")
#
# Figures (PDF vector + PNG 300 dpi + TIFF 600 dpi LZW), 183 mm wide (Nature double column):
#   Fig1_genomes          a genome size and gene number   b BUSCO, genome and proteome
#                         c A. latus subgenomes           d orthology relative to A. fumigatus
#                         e orthogroup sharing (UpSet)
#   Fig2_rnaseq           a read alignment and gene assignment per library
#                         b DEGs per strain (up / down)
#   Fig3_cross_species    a DEGs with an A. fumigatus orthologue
#                         b similarity of the responses (Jaccard, same direction)
#                         c DEGs shared between species (UpSet, A. fumigatus orthologues)
#   Fig4_genes            log2 fold change of the conserved core, the A. fumigatus-specific
#                         and the Nidulantes-specific genes
#   Fig5_hybrid           a subgenome of A. latus DEGs   b response of S-Q homeolog pairs
#   The enrichment figure (Fig. 6) is drawn by enrichment_figures.R.
#
# Tables (CSV, the data exactly as plotted) and figure_numbers.txt (every number
# quoted in the Results) are written to the same folder.
#
# Comparison across species: every DEG of A. spinulosporus and A. latus is converted
# to its A. fumigatus orthologue(s) (columns A_fumigatus_orthologs of the DEG sheets).
# Overlaps are counted on the 8,223 A. fumigatus genes that have orthologues in both
# other species, so that a gene absent from one genome cannot be counted as "not
# shared". For A. latus, an A. fumigatus gene takes the homeolog with the largest
# |log2FC| when both homeologs are DE.
#
# Needs: R >= 4.1 with readxl, ggplot2, patchwork, scales
# =============================================================================

# ---------------- settings (edit here) ----------------
TABLE2 <- "Supplementary_Tables_2.xlsx"
TABLE3 <- "Supplementary_Tables_3.xlsx"
OUTDIR <- NA                      # NA = folder "Paper_figures" next to TABLE3
BASE_SIZE <- 7                    # Nature: 5-7 pt text
WIDTH_MM  <- 183                  # double column; single column = 89
N_SPECIFIC <- 12                  # Fig. 4: genes shown per lineage-specific group

SPECIES <- c("A_fumigatus", "A_spinulosporus", "A_latus")
SP_LABEL <- c(A_fumigatus = "A. fumigatus", A_spinulosporus = "A. spinulosporus",
              A_latus = "A. latus")
SP_COL   <- c(A_fumigatus = "#E69F00", A_spinulosporus = "#0072B2", A_latus = "#009E73")
STRAINS  <- data.frame(
  species = c("A_fumigatus", "A_spinulosporus", "A_spinulosporus", "A_latus", "A_latus", "A_latus"),
  strain  = c("A1160", "3", "5845", "1", "1710", "3622"))
STRAINS$sheet <- paste0("DEGs_", STRAINS$strain)
STRAINS$label <- paste(SP_LABEL[STRAINS$species], STRAINS$strain)
DIR_COL  <- c(Up = "#B2182B", Down = "#2166AC")
SUB_COL  <- c(S = "#0072B2", Q = "#CC79A7", Ambiguous = "#999999", Unassigned = "#DDDDDD")

# ---------------- command line ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) TABLE2 <- args[1]
if (length(args) >= 2) TABLE3 <- args[2]
if (length(args) >= 3) OUTDIR <- args[3]
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
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---------------- helpers ----------------
num <- function(x) {                      # accepts numbers, "2,826" and "5,83E-40"
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x)); x[x %in% c("", "NA", "None", "-")] <- NA
  as.numeric(gsub(",", ".", x, fixed = TRUE))
}
sheet <- function(file, pattern) {         # sheet by name prefix, e.g. "S1 " or "DEGs_1$"
  sh <- readxl::excel_sheets(file)
  hit <- sh[grepl(pattern, sh)]
  if (!length(hit)) stop("no sheet matching '", pattern, "' in ", file, call. = FALSE)
  as.data.frame(readxl::read_excel(file, sheet = hit[1], guess_max = 1e5))
}
split_ids <- function(x) {
  x <- as.character(x); x[is.na(x) | x %in% c("None", "NA", "")] <- ""
  lapply(strsplit(x, ",", fixed = TRUE), trimws)
}
it <- function(x) parse(text = paste0("italic('", x, "')"))          # italic axis labels
it_strain <- function(lab) {                                           # "A. latus 1710" -> italic name + strain
  sp <- sub(" [^ ]+$", "", lab); st <- sub("^.* ", "", lab)
  parse(text = sprintf("italic('%s')~'%s'", sp, st))
}
theme_paper <- function() {
  theme_classic(base_size = BASE_SIZE) +
    theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          strip.background = element_blank(), strip.text = element_text(face = "bold"),
          legend.key.size = unit(3, "mm"), legend.title = element_text(size = BASE_SIZE),
          plot.tag = element_text(face = "bold", size = BASE_SIZE + 1),
          panel.grid = element_blank())
}
save_fig <- function(p, name, height_mm, width_mm = WIDTH_MM) {
  base <- file.path(OUTDIR, name)
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(paste0(base, ".pdf"), p, width = width_mm, height = height_mm, units = "mm", device = pdf_dev)
  ggsave(paste0(base, ".png"), p, width = width_mm, height = height_mm, units = "mm", dpi = 300, bg = "white")
  ggsave(paste0(base, ".tiff"), p, width = width_mm, height = height_mm, units = "mm", dpi = 600,
         compression = "lzw", bg = "white")
  message("  ", name, " (.pdf .png .tiff)")
}
write_tab <- function(d, name) write.csv(d, file.path(OUTDIR, paste0(name, ".csv")), row.names = FALSE)
NUMBERS <- character()
note <- function(...) { s <- sprintf(...); NUMBERS <<- c(NUMBERS, s); invisible(s) }

# UpSet-style plot: bars of intersection sizes over a dot matrix
upset_plot <- function(d, sets, ylab, fill = NULL, fill_values = NULL) {
  # d: data.frame with column 'pattern' (sets joined by "+"), 'n', optional fill column
  pats <- unique(d$pattern)
  ord  <- pats[order(-vapply(pats, function(p) sum(d$n[d$pattern == p]), 0))]
  d$pattern <- factor(d$pattern, levels = ord)
  bars <- ggplot(d, aes(pattern, n)) +
    { if (is.null(fill)) geom_col(width = 0.7, fill = "grey30")
      else geom_col(aes(fill = .data[[fill]]), width = 0.7, position = position_dodge(0.75)) } +
    { if (is.null(fill)) geom_text(aes(label = scales::comma(n)), vjust = -0.3, size = BASE_SIZE / 3)
      else geom_text(aes(label = scales::comma(n), group = .data[[fill]]),
                     position = position_dodge(0.75), vjust = -0.3, size = BASE_SIZE / 3.4) } +
    { if (!is.null(fill_values)) scale_fill_manual(values = fill_values, name = NULL) } +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::comma) +
    labs(x = NULL, y = ylab) + theme_paper() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          legend.position = "top")
  m <- expand.grid(pattern = factor(ord, levels = ord), set = sets, stringsAsFactors = FALSE)
  m$on <- mapply(function(p, s) s %in% strsplit(as.character(p), "+", fixed = TRUE)[[1]], m$pattern, m$set)
  m$set <- factor(m$set, levels = rev(sets))
  seg <- do.call(rbind, lapply(ord, function(p) {
    on <- m[m$pattern == p & m$on, ]
    if (nrow(on) < 2) return(NULL)
    data.frame(pattern = p, y1 = min(as.integer(on$set)), y2 = max(as.integer(on$set)))
  }))
  dots <- ggplot(m, aes(pattern, set)) +
    geom_point(aes(colour = on), size = 1.8, show.legend = FALSE) +
    { if (!is.null(seg)) geom_segment(data = seg, aes(x = pattern, xend = pattern, y = y1, yend = y2),
                                      inherit.aes = FALSE, linewidth = 0.5) } +
    scale_colour_manual(values = c(`TRUE` = "grey15", `FALSE` = "grey85")) +
    scale_y_discrete(labels = function(x) it(SP_LABEL[x])) +
    labs(x = NULL, y = NULL) + theme_paper() +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
  bars / dots + plot_layout(heights = c(3, 1.2))
}

message("Reading tables")
# =============================================================================
# Data
# =============================================================================
qc   <- sheet(TABLE2, "^S1 ");  qc   <- qc[qc$species %in% SPECIES, ]
ann  <- sheet(TABLE2, "^S5 ");  ann  <- ann[ann$species %in% SPECIES, ]
subg <- sheet(TABLE2, "^S3 ")
ost  <- sheet(TABLE2, "^S7 ")
orth <- sheet(TABLE2, "^S8 ")
shar <- sheet(TABLE2, "^S10 ")
cnt  <- sheet(TABLE2, "^S12 ")
desm <- sheet(TABLE2, "^S15 ")
smp  <- sheet(TABLE3, "^Samples")
degs <- lapply(STRAINS$sheet, function(s) sheet(TABLE3, paste0("^", s, "$")))
names(degs) <- STRAINS$strain
for (s in names(degs)) {
  degs[[s]]$log2FC <- num(degs[[s]]$log2FC)
  degs[[s]]$DE <- tolower(trimws(degs[[s]]$DE))
}

# A. fumigatus genes with orthologues in both other species (comparison universe)
orth$As_rel <- as.character(orth$A_spinulosporus_relation)
orth$Al_rel <- as.character(orth$A_latus_relation)
U3 <- orth$A_fumigatus[!(orth$As_rel %in% c("none", NA)) & !(orth$Al_rel %in% c("none", NA))]
note("A. fumigatus genes with orthologues in A. spinulosporus and A. latus: %d of %d (%.1f %%)",
     length(U3), nrow(orth), 100 * length(U3) / nrow(orth))

# per strain: A. fumigatus gene -> log2FC (DE genes only)
afua_fc <- function(st) {
  d <- degs[[st]]; d <- d[d$DE %in% c("up", "down"), ]
  ids <- if (st == "A1160") as.list(d$locus_tag) else split_ids(d$A_fumigatus_orthologs)
  rows <- rep(seq_len(nrow(d)), lengths(ids))
  x <- data.frame(afua = unlist(ids), fc = d$log2FC[rows], stringsAsFactors = FALSE)
  x <- x[x$afua != "", ]
  x <- x[order(-abs(x$fc)), ]
  x <- x[!duplicated(x$afua), ]             # one value per A. fumigatus gene (largest |log2FC|)
  setNames(x$fc, x$afua)
}
FC <- lapply(STRAINS$strain, afua_fc); names(FC) <- STRAINS$strain
none_both <- orth$A_fumigatus[orth$As_rel %in% "none" & orth$Al_rel %in% "none"]   # no orthologue in either
gene_set <- function(st, dirn, universe = U3) {
  f <- FC[[st]]; g <- names(f)[if (dirn == "up") f > 0 else f < 0]
  intersect(g, universe)
}

# =============================================================================
# Figure 1 - genomes, annotation, orthology
# =============================================================================
message("Figure 1")
sp_f <- function(x) factor(x, levels = SPECIES)
f1a_d <- rbind(data.frame(species = qc$species, metric = "Genome size (Mb)", value = num(qc$size_Mb)),
               data.frame(species = ann$species, metric = "Protein-coding genes", value = num(ann$genes)))
f1a_d$species <- sp_f(f1a_d$species)
f1a <- ggplot(f1a_d, aes(species, value, fill = species)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(value, accuracy = ifelse(value < 100, 0.1, 1))),
            vjust = -0.3, size = BASE_SIZE / 3) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = SP_COL) + scale_x_discrete(labels = function(x) it(SP_LABEL[x])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::comma) +
  labs(x = NULL, y = NULL) + theme_paper() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

busco <- function(d, what) {
  do.call(rbind, lapply(c("S", "D", "F", "M"), function(k)
    data.frame(species = d$species, set = what, class = k, value = num(d[[paste0("BUSCO_", k)]]))))
}
f1b_d <- rbind(busco(qc, "Genome"), busco(ann, "Proteome"))
f1b_d$species <- factor(f1b_d$species, levels = rev(SPECIES))
f1b_d$class <- factor(f1b_d$class, levels = c("M", "F", "D", "S"),
                      labels = c("Missing", "Fragmented", "Duplicated", "Single-copy"))
f1b <- ggplot(f1b_d, aes(value, species, fill = class)) +
  geom_col(width = 0.7) + facet_wrap(~set, ncol = 1) +
  scale_fill_manual(values = c(`Single-copy` = "#56B4E9", Duplicated = "#0072B2",
                               Fragmented = "#F0E442", Missing = "#D55E00"), name = "BUSCO",
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(labels = function(x) it(SP_LABEL[x])) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
  labs(x = "BUSCO genes (%)", y = NULL) + theme_paper() + theme(legend.position = "bottom")

sv <- setNames(num(subg$value), subg$metric)
al_ann <- ann[ann$species == "A_latus", ]
f1c_d <- rbind(
  data.frame(level = "Sequence (Mb)", class = c("S", "Q", "Ambiguous", "Unassigned"),
             value = sv[c("S_Mb", "Q_Mb", "ambiguous_Mb", "unassigned_Mb")]),
  data.frame(level = "Genes", class = c("S", "Q", "Ambiguous", "Unassigned"),
             value = num(unlist(al_ann[c("subgenome_S", "subgenome_Q", "subgenome_ambiguous",
                                         "subgenome_unassigned")]))))
f1c_d$pct <- ave(f1c_d$value, f1c_d$level, FUN = function(v) 100 * v / sum(v))
f1c_d$class <- factor(f1c_d$class, levels = rev(names(SUB_COL)))
f1c_d$level <- factor(f1c_d$level, levels = c("Genes", "Sequence (Mb)"))
f1c <- ggplot(f1c_d, aes(pct, level, fill = class)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(class %in% c("S", "Q"),
                               sprintf("%s\n%s", class, scales::comma(value, accuracy = ifelse(value < 100, 0.1, 1))), "")),
            position = position_stack(vjust = 0.5), size = BASE_SIZE / 3.2, colour = "white", lineheight = 0.85) +
  scale_fill_manual(values = SUB_COL, name = NULL, guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(paste("Share of ", italic("A. latus"), " genome (%)")), y = NULL) +
  theme_paper() + theme(legend.position = "bottom")
note("A. latus subgenomes: S %.2f Mb (%.1f %%), Q %.2f Mb (%.1f %%); genes S %d, Q %d",
     sv["S_Mb"], sv["S_percent"], sv["Q_Mb"], sv["Q_percent"],
     as.integer(al_ann$subgenome_S), as.integer(al_ann$subgenome_Q))

rel_levels <- c("1:1", "1:2", "1:>2", "many:1", "many:many", "none")
f1d_d <- ost[ost$comparison %in% c("A_fumigatus -> A_spinulosporus", "A_fumigatus -> A_latus") &
             ost$category %in% rel_levels, ]
f1d_d$target <- factor(ifelse(grepl("spinul", f1d_d$comparison), "A_spinulosporus", "A_latus"),
                       levels = c("A_latus", "A_spinulosporus"))
f1d_d$category <- factor(f1d_d$category, levels = rev(rel_levels),
                         labels = rev(c("1:1", "1:2", "1:>2", "many:1", "many:many", "no orthologue")))
f1d_d$percent <- num(f1d_d$percent)
f1d <- ggplot(f1d_d, aes(percent, target, fill = category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(percent >= 8, sprintf("%.0f%%", percent), "")),
            position = position_stack(vjust = 0.5), size = BASE_SIZE / 3.2, colour = "white") +
  scale_fill_manual(values = setNames(c("#1B9E77", "#66A61E", "#A6D854", "#7570B3", "#BEAED4", "#BBBBBB"),
                                      c("1:1", "1:2", "1:>2", "many:1", "many:many", "no orthologue")),
                    name = expression(paste(italic("A. fumigatus"), " : other")),
                    guide = guide_legend(reverse = TRUE, nrow = 2)) +
  scale_y_discrete(labels = function(x) it(SP_LABEL[x])) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(paste(italic("A. fumigatus"), " genes (%)")), y = NULL) +
  theme_paper() + theme(legend.position = "bottom")

ab <- c(Af = "A_fumigatus", As = "A_spinulosporus", Al = "A_latus")
f1e_d <- data.frame(pattern = vapply(strsplit(shar$pattern, "+", fixed = TRUE),
                                     function(v) paste(ab[v], collapse = "+"), ""),
                    n = num(shar$orthogroups))
f1e <- upset_plot(f1e_d, SPECIES, "Orthogroups")
write_tab(rbind(f1a_d[, c("species", "metric", "value")]), "Fig1a_genome_size_genes")
write_tab(f1b_d, "Fig1b_BUSCO"); write_tab(f1c_d, "Fig1c_subgenomes")
write_tab(f1d_d, "Fig1d_orthology"); write_tab(f1e_d, "Fig1e_orthogroups")

fig1 <- (f1a + f1b + plot_layout(widths = c(1.1, 1))) /
        (f1c + f1d) /
        wrap_elements(full = f1e) +
        plot_layout(heights = c(1.1, 0.75, 0.95)) + plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 1))
save_fig(fig1, "Fig1_genomes", 200)

# =============================================================================
# Figure 2 - RNA-seq overview
# =============================================================================
message("Figure 2")
smp$species_key <- c(`Aspergillus fumigatus` = "A_fumigatus", `Aspergillus spinulosporus` = "A_spinulosporus",
                     `Aspergillus latus` = "A_latus")[trimws(smp$species)]
smp$key <- paste(smp$species_key, smp$sample_old)
cnt$key <- paste(cnt$species, cnt$sample)
cnt$strain <- as.character(smp$strain[match(cnt$key, smp$key)])
cnt$condition <- factor(ifelse(tolower(cnt$condition) == "control", "DMEM", "Macrophage"),
                        levels = c("DMEM", "Macrophage"))
cnt$label <- factor(STRAINS$label[match(paste(cnt$species, cnt$strain), paste(STRAINS$species, STRAINS$strain))],
                    levels = STRAINS$label)
f2a_d <- rbind(data.frame(cnt[, c("sample", "label", "condition")], metric = "Aligned to genome",
                          value = num(cnt$overall_alignment_pct)),
               data.frame(cnt[, c("sample", "label", "condition")], metric = "Assigned to genes",
                          value = num(cnt$pct_assigned_of_clean)))
low <- f2a_d$metric == "Aligned to genome" & f2a_d$value < 80
f2a <- ggplot(f2a_d, aes(label, value, colour = condition)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1), size = 1.1) +
  geom_text(data = f2a_d[low, ], aes(label = sample), colour = "black", size = BASE_SIZE / 3.2,
            hjust = -0.25, show.legend = FALSE) +
  facet_wrap(~metric) +
  scale_colour_manual(values = c(DMEM = "grey55", Macrophage = "#D55E00"), name = NULL) +
  scale_x_discrete(labels = it_strain) + scale_y_continuous(limits = c(50, 100)) +
  labs(x = NULL, y = "Read pairs (%)") + theme_paper() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")

desm$label <- factor(STRAINS$label[match(paste(desm$species, desm$strain), paste(STRAINS$species, STRAINS$strain))],
                     levels = rev(STRAINS$label))
f2b_d <- rbind(data.frame(label = desm$label, direction = "Up", n = num(desm$up)),
               data.frame(label = desm$label, direction = "Down", n = -num(desm$down)))
f2b_d$direction <- factor(f2b_d$direction, levels = c("Down", "Up"))
pct <- data.frame(label = desm$label, n = num(desm$up), txt = sprintf("%.1f%%", num(desm$percent_DEGs_of_tested)))
f2b <- ggplot(f2b_d, aes(n, label, fill = direction)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = scales::comma(abs(n)), hjust = ifelse(n > 0, -0.15, 1.15)), size = BASE_SIZE / 3) +
  geom_text(data = pct, aes(x = max(abs(f2b_d$n)) * 1.55, y = label, label = txt), inherit.aes = FALSE,
            size = BASE_SIZE / 3, colour = "grey30") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = DIR_COL, name = NULL, breaks = c("Up", "Down")) +
  scale_y_discrete(labels = it_strain) +
  scale_x_continuous(labels = function(x) scales::comma(abs(x)),
                     limits = c(-max(abs(f2b_d$n)) * 1.3, max(abs(f2b_d$n)) * 1.7)) +
  annotate("text", x = max(abs(f2b_d$n)) * 1.55, y = nrow(desm) + 0.75, label = "% tested",
           size = BASE_SIZE / 3, colour = "grey30", fontface = "bold") +
  coord_cartesian(clip = "off") +
  labs(x = "Differentially expressed genes", y = NULL) +
  theme_paper() + theme(legend.position = "top")
write_tab(f2a_d, "Fig2a_mapping"); write_tab(f2b_d, "Fig2b_DEG_numbers")
for (i in seq_len(nrow(desm)))
  note("DEGs %s %s: %s up, %s down, %s total (%.1f %% of tested)", desm$species[i], desm$strain[i],
       desm$up[i], desm$down[i], desm$DEGs[i], num(desm$percent_DEGs_of_tested[i]))
fig2 <- f2a + f2b + plot_layout(widths = c(1.4, 1)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 1))
save_fig(fig2, "Fig2_rnaseq", 75)

# =============================================================================
# Figure 3 - comparison of the responses through A. fumigatus orthologues
# =============================================================================
message("Figure 3")
f3a_d <- do.call(rbind, lapply(seq_len(nrow(STRAINS)), function(i) {
  st <- STRAINS$strain[i]; d <- degs[[st]]; d <- d[d$DE %in% c("up", "down"), ]
  ids <- if (st == "A1160") as.list(d$locus_tag) else split_ids(d$A_fumigatus_orthologs)
  cls_lev <- c("Shared by the three species", "Shared with one other species", "Species- or lineage-specific")
  if (st == "A1160") {
    none <- d$locus_tag %in% none_both
    cls <- ifelse(d$locus_tag %in% U3, cls_lev[1], ifelse(none, cls_lev[3], cls_lev[2]))
  } else {
    cls <- vapply(ids, function(v) {
      v <- v[v != ""]
      if (!length(v)) cls_lev[3] else if (any(v %in% U3)) cls_lev[1] else cls_lev[2]
    }, "")
  }
  tb <- table(cls)
  data.frame(label = STRAINS$label[i], class = names(tb), n = as.integer(tb))
}))
f3a_d$pct <- ave(f3a_d$n, f3a_d$label, FUN = function(v) 100 * v / sum(v))
f3a_d$label <- factor(f3a_d$label, levels = rev(STRAINS$label))
cls_lev <- c("Shared by the three species", "Shared with one other species", "Species- or lineage-specific")
f3a_d$class <- factor(f3a_d$class, levels = rev(cls_lev))
f3a <- ggplot(f3a_d, aes(pct, label, fill = class)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ifelse(pct >= 7, sprintf("%.0f%%", pct), "")),
            position = position_stack(vjust = 0.5), size = BASE_SIZE / 3.2, colour = "white") +
  scale_fill_manual(values = setNames(c("#332288", "#88CCEE", "#AA4499"), cls_lev), name = NULL,
                    guide = guide_legend(reverse = TRUE, ncol = 1)) +
  scale_y_discrete(labels = it_strain, expand = expansion(add = c(0.6, 0.6))) + scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "DEGs (%)", y = NULL) + theme_paper() +
  theme(legend.position = "top", legend.justification = "left")
for (l in STRAINS$label) {
  x <- f3a_d[f3a_d$label == l, ]
  note("%s: %d DEGs; shared by the three species %.1f %%; with one other species %.1f %%; species/lineage-specific %.1f %% (%d)",
       l, sum(x$n), sum(x$pct[x$class == cls_lev[1]]), sum(x$pct[x$class == cls_lev[2]]),
       sum(x$pct[x$class == cls_lev[3]]), sum(x$n[x$class == cls_lev[3]]))
}

signed <- function(st) {
  f <- FC[[st]]; f <- f[names(f) %in% U3]
  paste(names(f), ifelse(f > 0, "u", "d"))
}
SG <- lapply(STRAINS$strain, signed); names(SG) <- STRAINS$label
jac <- expand.grid(a = STRAINS$label, b = STRAINS$label, stringsAsFactors = FALSE)
jac$J <- mapply(function(a, b) length(intersect(SG[[a]], SG[[b]])) / length(union(SG[[a]], SG[[b]])), jac$a, jac$b)
jac$a <- factor(jac$a, levels = STRAINS$label); jac$b <- factor(jac$b, levels = rev(STRAINS$label))
jac_low <- jac[as.integer(jac$a) < (nrow(STRAINS) + 1 - as.integer(jac$b)), ]
f3b <- ggplot(jac_low, aes(a, b, fill = J)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", J), colour = J > 0.22), size = BASE_SIZE / 3, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08306B", name = "Jaccard", limits = c(0, NA)) +
  scale_x_discrete(labels = it_strain) + scale_y_discrete(labels = it_strain) +
  labs(x = NULL, y = NULL) + coord_equal() + theme_paper() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(),
        axis.ticks = element_blank(), legend.position = "right")
J <- function(a, b) jac$J[jac$a == a & as.character(jac$b) == b]
note("Jaccard (same direction, shared orthologues): A. fumigatus vs others %.2f-%.2f; As 3 vs As 5845 %.2f; As 5845 vs Al 3622 %.2f; As 3 vs Al 1 %.2f; Al 1710 vs Al 3622 %.2f",
     min(jac$J[jac$a == STRAINS$label[1] & as.character(jac$b) != STRAINS$label[1]]),
     max(jac$J[jac$a == STRAINS$label[1] & as.character(jac$b) != STRAINS$label[1]]),
     J(STRAINS$label[2], STRAINS$label[3]), J(STRAINS$label[3], STRAINS$label[6]),
     J(STRAINS$label[2], STRAINS$label[4]), J(STRAINS$label[5], STRAINS$label[6]))

sp_set <- function(sp, dirn, how = c("any", "all")) {
  how <- match.arg(how); sts <- STRAINS$strain[STRAINS$species == sp]
  s <- lapply(sts, gene_set, dirn = dirn)
  if (how == "any") Reduce(union, s) else Reduce(intersect, s)
}
upset_counts <- function(how) {
  do.call(rbind, lapply(c("up", "down"), function(dirn) {
    S <- lapply(SPECIES, sp_set, dirn = dirn, how = how); names(S) <- SPECIES
    allg <- Reduce(union, S)
    memb <- vapply(allg, function(g) paste(SPECIES[vapply(S, function(s) g %in% s, TRUE)], collapse = "+"), "")
    tb <- table(memb)
    data.frame(pattern = names(tb), n = as.integer(tb), direction = ifelse(dirn == "up", "Up", "Down"), how = how)
  }))
}
up_any <- upset_counts("any"); up_all <- upset_counts("all")
up_any$direction <- factor(up_any$direction, levels = c("Up", "Down"))
f3c <- upset_plot(up_any, SPECIES, expression(paste("DEGs (", italic("A. fumigatus"), " orthologues)")), fill = "direction", fill_values = DIR_COL)
for (h in list(up_any, up_all)) for (dd in c("Up", "Down")) {
  x <- h[h$direction == dd, ]
  note("Shared DEGs (%s strain(s) of each species), %s: all three %d; As+Al only %d; Af+As only %d; Af+Al only %d; Af only %d; As only %d; Al only %d",
       ifelse(h$how[1] == "any", ">= 1", "all"), dd,
       sum(x$n[x$pattern == paste(SPECIES, collapse = "+")]),
       sum(x$n[x$pattern == "A_spinulosporus+A_latus"]), sum(x$n[x$pattern == "A_fumigatus+A_spinulosporus"]),
       sum(x$n[x$pattern == "A_fumigatus+A_latus"]), sum(x$n[x$pattern == "A_fumigatus"]),
       sum(x$n[x$pattern == "A_spinulosporus"]), sum(x$n[x$pattern == "A_latus"]))
}
write_tab(f3a_d, "Fig3a_DEG_orthology"); write_tab(jac, "Fig3b_jaccard")
write_tab(rbind(up_any, up_all), "Fig3c_shared_DEGs")
fig3 <- (f3a + f3b + plot_layout(widths = c(1, 1.15))) / wrap_elements(full = f3c) +
  plot_layout(heights = c(1, 1.1)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 1))
save_fig(fig3, "Fig3_cross_species", 170)

# =============================================================================
# Figure 4 - conserved and lineage-specific genes
# =============================================================================
message("Figure 4")
af <- degs[["A1160"]]
short_name <- function(prod, pref) {
  p <- ifelse(is.na(pref) | pref %in% c("None", "-", ""), "", as.character(pref))
  d <- ifelse(is.na(prod) | prod %in% c("None", ""), "uncharacterized protein", as.character(prod))
  d <- sub("^(putative|probable) ", "", d)
  d <- ifelse(nchar(d) > 38, paste0(substr(d, 1, 36), "..."), d)
  has <- mapply(function(a, b) a != "" && !grepl(a, b, ignore.case = TRUE, fixed = FALSE), p, d)
  ifelse(has, paste0(d, " (", p, ")"), d)
}
af_name <- setNames(short_name(af$product, af$eggNOG_preferred_name), af$locus_tag)
anno_U3 <- setNames(as.character(orth$A_fumigatus), orth$A_fumigatus)

# (i) conserved core: same direction in every strain
core <- lapply(c("up", "down"), function(dirn) Reduce(intersect, lapply(STRAINS$strain, gene_set, dirn = dirn)))
core_genes <- unlist(core)
# (ii) A. fumigatus-specific: DE in A. fumigatus, no orthologue in either other species
af_spec <- af[af$DE %in% c("up", "down") & af$locus_tag %in% none_both, ]
af_spec <- af_spec[order(-abs(af_spec$log2FC)), ][seq_len(min(N_SPECIFIC, nrow(af_spec))), ]
# (iii) Nidulantes-specific: A. spinulosporus DEG without A. fumigatus orthologue whose
#       A. latus orthologue(s) change in the same direction in every A. latus strain
al_fc <- lapply(c("1", "1710", "3622"), function(st) {
  d <- degs[[st]]; d <- d[d$DE %in% c("up", "down"), ]; setNames(d$log2FC, d$locus_tag) })
names(al_fc) <- c("1", "1710", "3622")
nid <- do.call(rbind, lapply(c("3", "5845"), function(st) {
  d <- degs[[st]]; d <- d[d$DE %in% c("up", "down") & vapply(split_ids(d$A_fumigatus_orthologs),
                                                              function(v) all(v == ""), TRUE), ]
  d$strainAs <- st; d
}))
nid_ids <- unique(nid$locus_tag)
nid_tab <- do.call(rbind, lapply(nid_ids, function(g) {
  r <- nid[nid$locus_tag == g, ][1, ]
  alo <- split_ids(r$A_latus_orthologs)[[1]]; alo <- alo[alo != ""]
  fcs <- c(A1160 = NA,
           `3` = { x <- degs[["3"]]; v <- x$log2FC[x$locus_tag == g & x$DE %in% c("up", "down")]; if (length(v)) v[1] else NA },
           `5845` = { x <- degs[["5845"]]; v <- x$log2FC[x$locus_tag == g & x$DE %in% c("up", "down")]; if (length(v)) v[1] else NA })
  for (st in c("1", "1710", "3622")) {
    v <- al_fc[[st]][intersect(alo, names(al_fc[[st]]))]
    fcs[st] <- if (length(v)) v[which.max(abs(v))] else NA
  }
  data.frame(id = g, name = short_name(r$product, r$eggNOG_preferred_name), t(fcs), check.names = FALSE)
}))
dir_ok <- apply(nid_tab[, c("3", "5845", "1", "1710", "3622")], 1, function(v) {
  v <- v[!is.na(v)]; length(v) >= 4 && (all(v > 0) || all(v < 0)) })
nid_all <- nid_tab[dir_ok, ]
nid_all <- nid_all[!(nid_all$name %in% c("uncharacterized protein", "hypothetical protein")) |
                     nrow(nid_all) <= N_SPECIFIC, ]
nid_all <- nid_all[order(-rowMeans(abs(nid_all[, c("3", "5845", "1", "1710", "3622")]), na.rm = TRUE)), ]
nid_show <- head(nid_all, N_SPECIFIC)
note("Conserved core (same direction in all six strains): %d up, %d down", length(core[[1]]), length(core[[2]]))
note("A. fumigatus DEGs without orthologue in either other species: %d (%d up, %d down)",
     sum(af$DE %in% c("up", "down") & af$locus_tag %in% none_both),
     sum(af$DE == "up" & af$locus_tag %in% none_both), sum(af$DE == "down" & af$locus_tag %in% none_both))
note("Nidulantes-specific DEGs (no A. fumigatus orthologue, same direction in >= 4 of 5 Nidulantes strains): %d",
     nrow(nid_tab[dir_ok, ]))

long_rows <- function(ids, names, group, get) {
  do.call(rbind, lapply(seq_along(ids), function(i) {
    data.frame(id = ids[i], gene = paste0(names[i], "  ", ids[i]), group = group,
               label = STRAINS$label, fc = vapply(STRAINS$strain, function(st) get(ids[i], st), 0))
  }))
}
get_afua <- function(g, st) { v <- FC[[st]][g]; if (is.na(v)) NA_real_ else unname(v) }
rows_core <- long_rows(core_genes, af_name[core_genes], "Conserved in all strains", get_afua)
rows_af   <- long_rows(af_spec$locus_tag, af_name[af_spec$locus_tag], "A. fumigatus only", get_afua)
rows_nid  <- do.call(rbind, lapply(seq_len(nrow(nid_show)), function(i)
  data.frame(id = nid_show$id[i], gene = paste0(nid_show$name[i], "  ", nid_show$id[i]),
             group = "Nidulantes only", label = STRAINS$label,
             fc = unlist(nid_show[i, STRAINS$strain]))))
f4_d <- rbind(rows_core, rows_af, rows_nid)
f4_d$absent <- (f4_d$group == "A. fumigatus only" & f4_d$label != STRAINS$label[1]) |
               (f4_d$group == "Nidulantes only" & f4_d$label == STRAINS$label[1])
f4_d$group <- factor(f4_d$group, levels = c("Conserved in all strains", "A. fumigatus only", "Nidulantes only"),
                     labels = c("'Conserved in all strains'", "italic('A. fumigatus')~'only'", "italic('Nidulantes')~'only'"))
ord <- unique(f4_d[order(f4_d$group, -ave(f4_d$fc, f4_d$id, FUN = function(v) mean(v, na.rm = TRUE))), "gene"])
f4_d$gene <- factor(f4_d$gene, levels = rev(ord))
f4_d$label <- factor(f4_d$label, levels = STRAINS$label)
lim <- 6
f4_d$fc_plot <- pmax(pmin(f4_d$fc, lim), -lim)
f4 <- ggplot(f4_d, aes(label, gene)) +
  geom_tile(aes(fill = fc_plot), colour = "white", linewidth = 0.3) +
  geom_point(data = f4_d[f4_d$absent, ], shape = 4, size = 0.9, colour = "grey45") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", labeller = label_parsed) +
  scale_fill_gradient2(low = DIR_COL[["Down"]], mid = "white", high = DIR_COL[["Up"]], midpoint = 0,
                       limits = c(-lim, lim), na.value = "grey92",
                       name = expression(log[2]~"fold change")) +
  scale_x_discrete(labels = it_strain, position = "top") +
  labs(x = NULL, y = NULL, caption = "Grey: not differentially expressed; x: no orthologue in that genome. Colour capped at |6|.") +
  theme_paper() +
  theme(axis.text.x.top = element_text(angle = 40, hjust = 0), axis.line = element_blank(),
        axis.ticks = element_blank(), strip.text.y = element_text(angle = 0, hjust = 0),
        axis.text.y = element_text(size = BASE_SIZE - 1.5), plot.caption = element_text(hjust = 0, size = BASE_SIZE - 1))
write_tab(f4_d, "Fig4_genes")
save_fig(f4, "Fig4_genes", max(120, 3.1 * length(unique(f4_d$gene)) + 40))

# =============================================================================
# Figure 5 - the two A. latus subgenomes in the response
# =============================================================================
message("Figure 5")
bg <- num(unlist(al_ann[c("subgenome_S", "subgenome_Q")])); bg_other <- num(al_ann$genes) - sum(bg)
sub_lab <- function(x) ifelse(x %in% c("S", "Q"), x, "Other")
f5a_d <- rbind(
  data.frame(label = "Genome", direction = "", sub = c("S", "Q", "Other"), n = c(bg, bg_other)),
  do.call(rbind, lapply(c("1", "1710", "3622"), function(st) {
    d <- degs[[st]]; d <- d[d$DE %in% c("up", "down"), ]
    tb <- table(ifelse(d$DE == "up", "Up", "Down"), sub_lab(d$subgenome))
    x <- as.data.frame(tb, stringsAsFactors = FALSE); names(x) <- c("direction", "sub", "n")
    data.frame(label = STRAINS$label[STRAINS$strain == st], x)
  })))
f5a_d$pct <- ave(f5a_d$n, paste(f5a_d$label, f5a_d$direction), FUN = function(v) 100 * v / sum(v))
f5a_d$x <- factor(ifelse(f5a_d$label == "Genome", "Genome", paste(sub("^A. latus ", "", f5a_d$label), f5a_d$direction)),
                  levels = c("Genome", as.vector(t(outer(c("1", "1710", "3622"), c("Up", "Down"), paste)))))
f5a_d$sub <- factor(f5a_d$sub, levels = c("Other", "Q", "S"))
f5a <- ggplot(f5a_d, aes(x, pct, fill = sub)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 100 * bg[1] / (sum(bg) + bg_other), linetype = 2, linewidth = 0.3) +
  geom_text(aes(label = ifelse(sub != "Other", sprintf("%.0f", pct), "")),
            position = position_stack(vjust = 0.5), size = BASE_SIZE / 3.2, colour = "white") +
  scale_fill_manual(values = c(S = SUB_COL[["S"]], Q = SUB_COL[["Q"]], Other = "grey80"), name = "Subgenome",
                    breaks = c("S", "Q", "Other")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression(paste(italic("A. latus"), " strain and direction")), y = "Genes (%)") +
  theme_paper() + theme(axis.text.x = element_text(angle = 40, hjust = 1))

pairs <- orth[orth$Al_rel %in% "1:2", c("A_fumigatus", "A_latus_S", "A_latus_Q")]
pairs <- pairs[!is.na(pairs$A_latus_S) & !is.na(pairs$A_latus_Q) &
               !grepl(",", pairs$A_latus_S) & !grepl(",", pairs$A_latus_Q) &
               !(pairs$A_latus_S %in% "None") & !(pairs$A_latus_Q %in% "None"), ]
note("One-to-one S-Q homeolog pairs (A. fumigatus 1:2): %d", nrow(pairs))
f5b_d <- do.call(rbind, lapply(c("1", "1710", "3622"), function(st) {
  f <- al_fc[[st]]; s <- f[pairs$A_latus_S]; q <- f[pairs$A_latus_Q]
  cls <- ifelse(!is.na(s) & !is.na(q), ifelse(sign(s) == sign(q), "Both, same direction", "Opposite directions"),
         ifelse(!is.na(s), "S homeolog only", ifelse(!is.na(q), "Q homeolog only", NA)))
  tb <- table(cls)
  x <- data.frame(label = STRAINS$label[STRAINS$strain == st], class = names(tb), n = as.integer(tb))
  note("%s: %d pairs with >= 1 DE homeolog; both same direction %d; S only %d; Q only %d; opposite %d (one homeolog only: %.1f %%)",
       x$label[1], sum(x$n), sum(x$n[x$class == "Both, same direction"]), sum(x$n[x$class == "S homeolog only"]),
       sum(x$n[x$class == "Q homeolog only"]), sum(x$n[x$class == "Opposite directions"]),
       100 * sum(x$n[grepl("only", x$class)]) / sum(x$n))
  x
}))
pl <- c("Both, same direction", "S homeolog only", "Q homeolog only", "Opposite directions")
f5b_d$class <- factor(f5b_d$class, levels = rev(pl))
f5b_d$label <- factor(f5b_d$label, levels = STRAINS$label[4:6])
f5b <- ggplot(f5b_d, aes(label, n, fill = class)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ifelse(n >= 60, scales::comma(n), "")),
            position = position_stack(vjust = 0.5), size = BASE_SIZE / 3.2, colour = "white") +
  scale_fill_manual(values = setNames(c("#117733", SUB_COL[["S"]], SUB_COL[["Q"]], "black"), pl),
                    name = NULL, guide = guide_legend(reverse = TRUE)) +
  scale_x_discrete(labels = it_strain) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::comma) +
  labs(x = NULL, y = "S-Q homeolog pairs with a DE copy") + theme_paper()
write_tab(f5a_d, "Fig5a_subgenome_DEGs"); write_tab(f5b_d, "Fig5b_homeolog_pairs")
fig5 <- f5a + f5b + plot_layout(widths = c(1.2, 1)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 1))
save_fig(fig5, "Fig5_hybrid", 70)

writeLines(c(paste("paper_figures.R -", format(Sys.time(), "%Y-%m-%d %H:%M")),
             paste("Table 2:", normalizePath(TABLE2)), paste("Table 3:", normalizePath(TABLE3)), "", NUMBERS),
           file.path(OUTDIR, "figure_numbers.txt"))
message("Done: ", normalizePath(OUTDIR))
