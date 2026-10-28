#!/usr/bin/env Rscript
# =============================================================================
# enrichment_figures.R  -  publication figures from the curated enrichment panel
#
# Input : Supplementary_Tables_3.xlsx (curated enrichment table), columns
#         species, strain, direction, collection, term, name, DEGs_in_term,
#         DEGs_annotated, term_size, universe, gene_ratio, fold_enrichment,
#         pvalue, padj, significant, genes
#         (numbers with decimal comma, e.g. "2,826" / "5,83E-40", are accepted)
#
# Run   : Rscript enrichment_figures.R
#         Rscript enrichment_figures.R <table.xlsx> <sheet name or number> <output folder>
#         or in RStudio: source("enrichment_figures.R")
#
# Figures (PDF vector + PNG 300 dpi + TIFF 600 dpi, LZW) in <output folder>:
#   Fig_enrichment_dotplot_all      every term x every strain/direction, all species
#                                   (rows = terms grouped by collection, columns =
#                                    species > strain > up/down; size = gene ratio,
#                                    colour = -log10 adjusted p)
#   Fig_enrichment_dotplot_<sp>     the same, one species per figure
#   Fig_enrichment_bars_<sp>        diverging bars per strain: up to the right,
#                                   down to the left, length = -log10 adjusted p
#   Fig_enrichment_counts           number of enriched terms per strain,
#                                   direction and collection
#   enrichment_panel_plotted.csv    the table exactly as plotted
#
# Needs: R >= 4.1, packages readxl, ggplot2, scales
# =============================================================================

# ---------------- settings (edit here) ----------------
INFILE  <- "Supplementary_Tables_3.xlsx"   # curated panel; or give the path as the first argument
SHEET   <- NA          # NA = first sheet that has the enrichment columns; or a name / number
OUTDIR  <- NA          # NA = folder "Enrichment_figures" next to the table

ONLY_SIGNIFICANT <- TRUE   # keep only rows with padj < PADJ_MAX (the panel is already curated)
PADJ_MAX         <- 0.05
WRAP             <- 48     # characters per line in term names
SPECIES_ORDER    <- c("A_fumigatus", "A_spinulosporus", "A_latus")
COLLECTION_LABEL <- c(GO_BP = "GO biological process", GO_MF = "GO molecular function",
                      GO_CC = "GO cellular component", KEGG_pathway = "KEGG pathway",
                      COG = "COG category")
DIR_COLOUR       <- c(Up = "#B2182B", Down = "#2166AC")
FILL_SCALE       <- c("#FEE8C8", "#FDBB84", "#E34A33", "#7F0000")   # -log10 padj, low -> high
BASE_SIZE        <- 9

# ---------------- command line ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) INFILE <- args[1]
if (length(args) >= 2) SHEET  <- if (grepl("^[0-9]+$", args[2])) as.integer(args[2]) else args[2]
if (length(args) >= 3) OUTDIR <- args[3]
if (is.na(OUTDIR)) OUTDIR <- file.path(dirname(INFILE), "Enrichment_figures")

# ---------------- packages ----------------
need <- c("readxl", "ggplot2", "scales")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss))
  stop("missing R packages: ", paste(miss, collapse = ", "),
       "\n  install with: install.packages(c(", paste0('"', miss, '"', collapse = ", "), "))",
       "\n  or (conda/micromamba env): micromamba install -c conda-forge ",
       paste0("r-", miss, collapse = " "), call. = FALSE)
suppressPackageStartupMessages({ library(ggplot2) })
if (!file.exists(INFILE)) stop("table not found: ", INFILE, call. = FALSE)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------- read ----------------
REQUIRED <- c("species", "strain", "direction", "collection", "name", "degs_in_term",
              "degs_annotated", "padj")
clean_names <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

read_table <- function(path, sheet) {
  if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE))
    return(read.delim(path, colClasses = "character", check.names = FALSE))
  if (grepl("\\.csv$", path, ignore.case = TRUE))
    return(read.csv(path, colClasses = "character", check.names = FALSE))
  sheets <- readxl::excel_sheets(path)
  if (is.na(sheet)) {
    ok <- vapply(sheets, function(s) {
      h <- clean_names(names(readxl::read_excel(path, sheet = s, n_max = 0)))
      all(REQUIRED %in% h)
    }, logical(1))
    if (!any(ok)) stop("no sheet with the enrichment columns (", paste(REQUIRED, collapse = ", "),
                       ") in ", path, "; sheets: ", paste(sheets, collapse = ", "), call. = FALSE)
    sheet <- sheets[which(ok)[1]]
  }
  message("Reading sheet '", if (is.numeric(sheet)) sheets[sheet] else sheet, "' of ", basename(path))
  as.data.frame(readxl::read_excel(path, sheet = sheet, col_types = "text"))
}

# numbers stored as text with decimal comma ("2,826", "5,83E-40") or as numbers
num <- function(x) {
  x <- gsub("[[:space:]]", "", as.character(x))
  both <- grepl("\\.", x) & grepl(",", x)             # "1.234,5" -> thousands point
  x[both] <- gsub("\\.", "", x[both])
  suppressWarnings(as.numeric(gsub(",", ".", x)))
}

d <- read_table(INFILE, SHEET)
names(d) <- clean_names(names(d))
miss_col <- setdiff(REQUIRED, names(d))
if (length(miss_col)) stop("columns missing: ", paste(miss_col, collapse = ", "), call. = FALSE)
d <- d[!is.na(d$name) & nzchar(trimws(d$name)), ]
for (v in intersect(c("degs_in_term", "degs_annotated", "term_size", "universe",
                      "fold_enrichment", "pvalue", "padj"), names(d))) d[[v]] <- num(d[[v]])
if (any(is.na(d$padj))) warning(sum(is.na(d$padj)), " rows with unreadable padj were removed")
d <- d[!is.na(d$padj) & !is.na(d$degs_in_term) & !is.na(d$degs_annotated), ]
if (ONLY_SIGNIFICANT) d <- d[d$padj < PADJ_MAX, ]
if (!nrow(d)) stop("no rows left to plot (check ONLY_SIGNIFICANT / PADJ_MAX)", call. = FALSE)

# ---------------- derived columns ----------------
d$species   <- trimws(d$species); d$strain <- trimws(d$strain); d$collection <- trimws(d$collection)
d$direction <- ifelse(tolower(trimws(d$direction)) %in% c("up", "up-regulated", "upregulated"), "Up", "Down")
d$gene_ratio_num <- d$degs_in_term / d$degs_annotated        # computed: Excel may turn "138/755" into a date
d$mlog10_padj    <- -log10(pmax(d$padj, 1e-300))
sp_levels <- c(intersect(SPECIES_ORDER, unique(d$species)), setdiff(unique(d$species), SPECIES_ORDER))
d$species_f <- factor(d$species, levels = sp_levels)
d$species_lab <- factor(sub("_", ". ", d$species), levels = sub("_", ". ", sp_levels))
coll_levels <- c(intersect(names(COLLECTION_LABEL), unique(d$collection)),
                 setdiff(unique(d$collection), names(COLLECTION_LABEL)))
lab_coll <- function(x) ifelse(x %in% names(COLLECTION_LABEL), COLLECTION_LABEL[x], x)
wrap_coll <- function(x) vapply(x, function(s) paste(strwrap(s, 14), collapse = "\n"), character(1))
d$collection_lab <- factor(wrap_coll(lab_coll(d$collection)), levels = unique(wrap_coll(lab_coll(coll_levels))))

# strains in natural order inside each species, then Up before Down
strain_key <- function(s) { n <- suppressWarnings(as.numeric(s)); ifelse(is.na(n), Inf, n) }
o <- order(d$species_f, strain_key(d$strain), d$strain, d$direction != "Up")
col_ids <- unique(paste(d$species[o], d$strain[o], d$direction[o], sep = "|"))
d$column <- factor(paste(d$species, d$strain, d$direction, sep = "|"), levels = col_ids)
col_labels <- setNames(sub("^[^|]*\\|([^|]*)\\|(.*)$", "\\1\n\\2", col_ids), col_ids)
d$strain_lab <- factor(paste("Strain", d$strain),
                       levels = unique(paste("Strain", d$strain[o])))

wrap_term <- function(x) vapply(x, function(s) paste(strwrap(s, WRAP), collapse = "\n"), character(1))
d$term_lab <- wrap_term(trimws(d$name))

# term order inside each collection: enriched in more strains first, then strongest
order_terms <- function(df) {
  key <- unique(df[, c("collection_lab", "term_lab")])
  key$n   <- mapply(function(c, t) sum(df$collection_lab == c & df$term_lab == t), key$collection_lab, key$term_lab)
  key$max <- mapply(function(c, t) max(df$mlog10_padj[df$collection_lab == c & df$term_lab == t]),
                    key$collection_lab, key$term_lab)
  key <- key[order(key$collection_lab, -key$n, -key$max), ]
  rev(unique(key$term_lab))                                   # first term on top
}
write.csv(d[order(d$column, d$collection_lab, d$padj),
            setdiff(names(d), c("column", "species_f"))],
          file.path(OUTDIR, "enrichment_panel_plotted.csv"), row.names = FALSE)

# ---------------- theme + saving ----------------
theme_pub <- theme_bw(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        strip.background = element_rect(fill = "grey95", colour = "grey60", linewidth = 0.3),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
        axis.text.y = element_text(lineheight = 0.85, colour = "black"),
        axis.text.x = element_text(colour = "black", lineheight = 0.9),
        legend.key.size = unit(0.4, "cm"), legend.title = element_text(size = BASE_SIZE - 1),
        plot.title = element_text(face = "bold", size = BASE_SIZE + 2))

save_fig <- function(p, name, w, h) {
  w <- min(max(w, 4), 18); h <- min(max(h, 3), 30)
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(file.path(OUTDIR, paste0(name, ".pdf")), p, width = w, height = h, units = "in", device = pdf_dev)
  ggsave(file.path(OUTDIR, paste0(name, ".png")), p, width = w, height = h, units = "in", dpi = 300)
  ggsave(file.path(OUTDIR, paste0(name, ".tiff")), p, width = w, height = h, units = "in", dpi = 600,
         compression = "lzw")
  message(sprintf("  %-40s %.1f x %.1f in", name, w, h))
}
term_width <- function(df) 0.055 * min(WRAP, max(nchar(trimws(df$name))))   # inches for the y labels

fill_scale <- scale_fill_gradientn(colours = FILL_SCALE, name = expression(-log[10]~italic(P)[adj]))
size_scale <- scale_size_area(max_size = 6, name = "Gene ratio\n(DEGs in term /\nannotated DEGs)",
                              labels = scales::percent_format(accuracy = 1))

dotplot <- function(df, title) {
  df$term_f <- factor(df$term_lab, levels = order_terms(df))
  # light grey band behind every down-regulated column, in every panel
  cols <- unique(df[order(df$column), c("species_lab", "column", "direction")])
  cols$i <- ave(seq_len(nrow(cols)), cols$species_lab, FUN = seq_along)
  bg <- merge(cols[cols$direction == "Down", ], data.frame(collection_lab = unique(df$collection_lab)))
  p <- ggplot(df, aes(column, term_f)) +
    geom_rect(data = bg, aes(xmin = i - 0.5, xmax = i + 0.5, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey93") +
    geom_point(aes(size = gene_ratio_num, fill = mlog10_padj), shape = 21, colour = "grey25", stroke = 0.3) +
    fill_scale + size_scale +
    scale_x_discrete(labels = col_labels) +
    labs(x = NULL, y = NULL, title = title) + theme_pub +
    theme(strip.text.x = element_text(face = "italic"))
  p + facet_grid(collection_lab ~ species_lab, scales = "free", space = "free")
}
dot_size <- function(df) {
  nt <- length(unique(paste(df$collection_lab, df$term_lab)))
  nl <- sum(vapply(unique(df$term_lab), function(s) length(strsplit(s, "\n")[[1]]), integer(1)))
  c(w = 2.6 + term_width(df) + 0.42 * length(unique(df$column)),
    h = 1.4 + 0.16 * nl + 0.06 * nt + 0.25 * length(unique(df$collection_lab)))
}

# ---------------- Figure 1: all species ----------------
message("Writing figures to ", OUTDIR)
s <- dot_size(d)
save_fig(dotplot(d, "Functional enrichment of DEGs (macrophage vs DMEM)"),
         "Fig_enrichment_dotplot_all", s["w"], s["h"])

# ---------------- Figure 2: one dot plot per species ----------------
for (sp in sp_levels) {
  ds <- d[d$species == sp, ]
  s <- dot_size(ds)
  save_fig(dotplot(ds, NULL) + ggtitle(sub("_", ". ", sp)) + theme(plot.title = element_text(face = "bold.italic")),
           paste0("Fig_enrichment_dotplot_", sp), s["w"], s["h"])
}

# ---------------- Figure 3: diverging bars per species ----------------
for (sp in sp_levels) {
  ds <- d[d$species == sp, ]
  ds$term_f <- factor(ds$term_lab, levels = order_terms(ds))
  ds$signed <- ifelse(ds$direction == "Up", 1, -1) * ds$mlog10_padj
  ds$direction <- factor(ds$direction, levels = c("Up", "Down"))
  lim <- max(abs(ds$signed)) * 1.18
  p <- ggplot(ds, aes(signed, term_f, fill = direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = c(-1, 1) * -log10(PADJ_MAX), linetype = "dashed", linewidth = 0.25, colour = "grey50") +
    geom_text(aes(label = degs_in_term, hjust = ifelse(signed > 0, -0.25, 1.25)), size = BASE_SIZE / 3.6) +
    scale_fill_manual(values = DIR_COLOUR, name = NULL, labels = c(Up = "Up-regulated", Down = "Down-regulated")) +
    scale_x_continuous(limits = c(-lim, lim), labels = function(x) abs(x),
                       expand = expansion(mult = 0.02)) +
    labs(x = expression(-log[10]~italic(P)[adj]~"(down-regulated  |  up-regulated)"), y = NULL,
         title = sub("_", ". ", sp), caption = "Numbers = DEGs in the term") +
    facet_grid(collection_lab ~ strain_lab, scales = "free_y", space = "free_y") +
    theme_pub + theme(plot.title = element_text(face = "bold.italic"), legend.position = "top")
  nl <- sum(vapply(unique(ds$term_lab), function(s) length(strsplit(s, "\n")[[1]]), integer(1)))
  save_fig(p, paste0("Fig_enrichment_bars_", sp),
           2.4 + term_width(ds) + 2.3 * length(unique(ds$strain)),
           1.6 + 0.16 * nl + 0.3 * length(unique(ds$collection_lab)))
}

# ---------------- Figure 4: number of enriched terms ----------------
cnt <- aggregate(list(terms = rep(1, nrow(d))),
                 by = list(species_lab = d$species_lab, strain = d$strain, direction = d$direction,
                           collection_lab = d$collection_lab, column = d$column), FUN = sum)
cnt$direction <- factor(cnt$direction, levels = c("Up", "Down"))
p <- ggplot(cnt, aes(column, terms, fill = collection_lab)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  scale_x_discrete(labels = col_labels) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E7298A", "#D95F02", "#666666")[seq_along(levels(d$collection_lab))],
                    name = NULL, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  facet_grid(. ~ species_lab, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Enriched terms in the panel") + theme_pub +
  theme(strip.text.x = element_text(face = "italic"), legend.position = "right")
save_fig(p, "Fig_enrichment_counts", 3.5 + 0.45 * length(unique(d$column)), 3.6)

message("Done: ", nrow(d), " enriched terms, ", length(unique(d$column)), " strain x direction columns, ",
        length(sp_levels), " species.")
