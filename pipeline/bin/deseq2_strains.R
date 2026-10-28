#!/usr/bin/env Rscript
# ============================================================
# deseq2_strains.R  counts.tsv  metadata.csv  species  outdir  [min_count=10]
#
# One DESeq2 model per species with all its samples:
#   design ~ group,  group = <strain>_<treatment>  (treatment: DMEN | Macrophage)
# Contrast, for EVERY strain: <strain>_Macrophage vs <strain>_DMEN (control).
# A single model shares the dispersion estimates across strains (more power,
# standard practice); each contrast only compares samples of the same strain.
#
#  - pre-filter: >= min_count reads in >= 3 samples (smallest group = 3)
#  - Wald test, BH-adjusted p (independent filtering on)
#  - log2FC shrinkage: apeglm (Zhu et al. 2018) if installed, else "normal";
#    the unshrunken MLE log2FC is also written
#  - vst(blind = FALSE) for PCA / heatmaps
# Outputs (outdir): DE_<strain>.tsv (all tested genes), normalized_counts.tsv,
#   vst.tsv, size_factors.tsv, samples_used.tsv, deseq2_info.txt
# ============================================================
suppressPackageStartupMessages(library(DESeq2))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 4) stop("usage: deseq2_strains.R counts.tsv metadata.csv species outdir [min_count]")
cnt_f <- a[1]; meta_f <- a[2]; sp <- a[3]; out <- a[4]
min_count <- if (length(a) >= 5) as.numeric(a[5]) else 10
dir.create(out, recursive = TRUE, showWarnings = FALSE)

counts <- as.matrix(read.delim(cnt_f, row.names = 1, check.names = FALSE))
meta <- read.csv(meta_f, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", fileEncoding = "UTF-8-BOM")
meta <- meta[meta$species_id == sp, ]
miss <- setdiff(meta$sample_old, colnames(counts))
if (length(miss)) stop("samples in metadata but not in the count matrix: ", paste(miss, collapse = ", "))
meta <- meta[order(meta$strain, meta$treatment != "DMEN", meta$replicate), ]
rownames(meta) <- meta$sample_old
counts <- counts[, rownames(meta)]
if (!all(meta$treatment %in% c("DMEN", "Macrophage"))) stop("treatment must be DMEN or Macrophage")
# syntactically valid factor levels (strain names can start with a digit)
meta$grp <- factor(make.names(paste0("s", meta$strain, "_", meta$treatment)))
strains <- unique(meta$strain)

keep <- rowSums(counts >= min_count) >= 3
cat(sprintf("%s: %d genes, %d kept (>= %d reads in >= 3 samples); %d samples, strains: %s\n",
            sp, nrow(counts), sum(keep), min_count, ncol(counts), paste(strains, collapse = ", ")))
dds <- DESeqDataSetFromMatrix(counts[keep, ], meta, ~ grp)
dds <- DESeq(dds, quiet = TRUE)
shrink <- if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else "normal"

norm <- counts(dds, normalized = TRUE)
write.table(data.frame(locus_tag = rownames(norm), round(norm, 3), check.names = FALSE),
            file.path(out, "normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
vsd <- vst(dds, blind = FALSE)
write.table(data.frame(locus_tag = rownames(vsd), round(assay(vsd), 4), check.names = FALSE),
            file.path(out, "vst.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = colnames(dds), size_factor = round(sizeFactors(dds), 4)),
            file.path(out, "size_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta, file.path(out, "samples_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

for (st in strains) {
  ctl <- make.names(paste0("s", st, "_DMEN")); trt <- make.names(paste0("s", st, "_Macrophage"))
  if (!all(c(ctl, trt) %in% levels(dds$grp))) { cat("  skip strain", st, "(needs DMEN and Macrophage)\n"); next }
  d2 <- dds; d2$grp <- relevel(d2$grp, ref = ctl); d2 <- nbinomWaldTest(d2, quiet = TRUE)
  coef <- paste0("grp_", trt, "_vs_", ctl)
  res <- results(d2, name = coef)
  shr <- lfcShrink(d2, coef = coef, type = shrink, quiet = TRUE)
  nc <- counts(d2, normalized = TRUE)
  s_ctl <- rownames(meta)[meta$strain == st & meta$treatment == "DMEN"]
  s_trt <- rownames(meta)[meta$strain == st & meta$treatment == "Macrophage"]
  df <- data.frame(locus_tag = rownames(res), baseMean = round(res$baseMean, 3),
                   log2FC = round(shr$log2FoldChange, 4), log2FC_MLE = round(res$log2FoldChange, 4),
                   lfcSE = round(res$lfcSE, 4), stat = round(res$stat, 4),
                   pvalue = signif(res$pvalue, 4), padj = signif(res$padj, 4),
                   mean_norm_DMEN = round(rowMeans(nc[, s_ctl, drop = FALSE]), 2),
                   mean_norm_Macrophage = round(rowMeans(nc[, s_trt, drop = FALSE]), 2))
  write.table(df, file.path(out, paste0("DE_", st, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  strain %-6s Macrophage vs DMEN: %d genes tested (padj not NA)\n", st, sum(!is.na(res$padj))))
}
writeLines(c(paste("DESeq2", as.character(packageVersion("DESeq2"))),
             paste("shrinkage", shrink, if (shrink == "apeglm") as.character(packageVersion("apeglm")) else ""),
             paste("R", R.version$major, R.version$minor, sep = "."),
             "design ~ group (strain_treatment); contrast per strain: Macrophage vs DMEN",
             sprintf("pre-filter: >= %d reads in >= 3 samples", min_count)),
           file.path(out, "deseq2_info.txt"))
