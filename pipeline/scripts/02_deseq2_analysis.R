#!/usr/bin/env Rscript
#
# DESeq2 analysis for RNA-seq time course experiment
# Reads experimental design from design.yaml
# Part A: Likelihood Ratio Test (overall time course effect)
# Part B: Pairwise Wald contrasts (auto-generated from design)
#

# ─── Snakemake integration ───
if (exists("snakemake")) {
  out_dir <- dirname(snakemake@output[["combined"]])
  alpha_val <- 0.05
} else {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args) >= 1) args[1] else "results"
  alpha_val <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
}

suppressPackageStartupMessages({
  library("DESeq2")
  library("ggplot2")
  library("pheatmap")
  library("RColorBrewer")
  library("dplyr")
  library("ggrepel")
  library("jsonlite")
})

`%||%` <- function(a, b) if (is.null(a) || is.na(a) || length(a) == 0) b else a

cat("=== DESeq2 Analysis ===\n")
cat(sprintf("Output dir: %s\n", out_dir))
# --- 1. Load data ---
counts <- as.matrix(read.table(file.path(out_dir, "counts_matrix.tsv"),
                                header = TRUE, row.names = 1, sep = "\t",
                                check.names = FALSE))
metadata <- read.table(file.path(out_dir, "metadata.tsv"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
annotations <- read.table(file.path(out_dir, "gene_annotations.tsv"),
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                          quote = "", fill = TRUE, comment.char = "")

cat(sprintf("Counts: %d genes x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("Metadata: %d samples\n", nrow(metadata)))

# Ensure sample order matches
stopifnot(all(colnames(counts) == metadata$sample_id))

# ─── Read experimental design ───
if (exists("snakemake")) {
  design <- jsonlite::fromJSON(snakemake@params$design)
} else {
  args <- commandArgs(trailingOnly = TRUE)
  design_file <- if (length(args) >= 2) args[2] else "../data/design.yaml"
  design <- yaml::read_yaml(design_file)
}

f <- design$factors
cl_list   <- f$cell_lines
tp_list   <- f$time_points
treat_id  <- f$treatment$id
comp      <- design$comparisons

cl_ids     <- sapply(cl_list, `[[`, "id")
tp_ids     <- sapply(tp_list, `[[`, "id")

# Helper: short label for a time point (keyed by original IDs)
tp_short <- setNames(sapply(tp_list, `[[`, "short"), tp_ids)
cl_short <- setNames(sapply(cl_list, `[[`, "short"), cl_ids)

ref_cl     <- f$reference_cell_line %||% cl_ids[1]
ref_tp     <- f$reference_time_point %||% tp_ids[1]

# Reorder so reference is first (DESeq2 uses first level as baseline)
cl_ids     <- c(ref_cl, setdiff(cl_ids, ref_cl))
tp_ids     <- c(ref_tp, setdiff(tp_ids, ref_tp))
nonref_cl  <- setdiff(cl_ids, ref_cl)[1]
nonref_tps <- setdiff(tp_ids, ref_tp)

metadata$cell_line <- factor(metadata$cell_line, levels = cl_ids)
metadata$time      <- factor(metadata$time, levels = tp_ids)
metadata$batch     <- factor(metadata$batch)
metadata$group     <- factor(paste(metadata$cell_line, metadata$time, sep = "_"))

rownames(metadata) <- metadata$sample_id

# --- 2. Build DESeq2 object ---
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = metadata,
                               design = ~ batch + cell_line + time + cell_line:time)

cat(sprintf("\nDesign: ~ batch + cell_line + time + cell_line:time\n"))
cat(sprintf("  cell_lines: %s (ref: %s)\n", paste(cl_ids, collapse=", "), ref_cl))
cat(sprintf("  time_points: %s (ref: %s)\n", paste(tp_ids, collapse=", "), ref_tp))

# --- 3. Part A: LRT ---
cat("\n=== Part A: Likelihood Ratio Test ===\n")

dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ batch + cell_line)
res_lrt <- results(dds_lrt, alpha = alpha_val)
res_lrt <- res_lrt[order(res_lrt$pvalue), ]

cat(sprintf("LRT: %d genes with padj < %.2g\n", sum(res_lrt$padj < alpha_val, na.rm = TRUE), alpha_val))
cat(sprintf("LRT: %d genes with padj < 0.01\n", sum(res_lrt$padj < 0.01, na.rm = TRUE)))

# --- 4. Part B: Pairwise Wald contrasts ---
dds_wald <- DESeq(dds, test = "Wald")

coef_names <- resultsNames(dds_wald)
cat("\nModel coefficients:\n")
print(coef_names)

# --- Dynamic coefficient extraction ---
# Time main effects
coef_time <- list()
for (tp in nonref_tps) {
  pat <- paste0("time_?", tp, "_vs_", ref_tp)
  coef_time[[tp]] <- grep(pat, coef_names, value = TRUE)[1]
  cat(sprintf("  time %s: %s\n", tp, coef_time[[tp]]))
}

# Cell line main effect
coef_cell_line <- grep(paste0("cell_line_?", nonref_cl), coef_names, value = TRUE)[1]
cat(sprintf("  cell_line %s: %s\n", nonref_cl, coef_cell_line))

# Interaction terms
coef_int <- list()
for (tp in nonref_tps) {
  pat <- paste0("cell.line.*", nonref_cl, ".*\\.time", ".*", tp)
  coef_int[[tp]] <- grep(pat, coef_names, value = TRUE)[1]
  cat(sprintf("  interact %s: %s\n", tp, coef_int[[tp]] %||% "(none)"))
}

# Helper: extract results with lfcShrink
extract_contrast <- function(dds, name = NULL, contrast = NULL, label = "") {
  if (!is.null(name)) {
    res <- lfcShrink(dds, coef = name, type = "ashr", quiet = TRUE)
  } else {
    res <- lfcShrink(dds, contrast = contrast, type = "ashr", quiet = TRUE)
  }
  cat(sprintf("  %-50s: %5d sig (padj < %.2g)\n",
              label, sum(res$padj < alpha_val, na.rm = TRUE), alpha_val))
  res
}

# --- Dynamic contrast generation ---
cat("\n=== Part B: Pairwise Contrasts ===\n")

contrast_list <- list()
contrast_labels <- c()
contrast_short  <- c()

# --- Within cell line: {cl}_{time}_vs_mock ---
if (isTRUE(comp$within_cell_line)) {
  for (cl in cl_ids) {
    cat(sprintf("\n-- Within %s --\n", cl))
    for (tp in nonref_tps) {
      cname <- paste0(cl, "_", tp, "_vs_mock")
      clabel <- paste0(cl, " ", tp, " vs mock")
      cshort <- paste0(cl_short[cl], tp_short[tp], "m")

      if (cl == ref_cl) {
        contrast_list[[cname]] <- extract_contrast(
          dds_wald, contrast = c("time", tp, ref_tp), label = clabel)
      } else {
        ct <- coef_time[[tp]]
        ci <- coef_int[[tp]]
        if (!is.na(ci) && ci != "") {
          contrast_list[[cname]] <- extract_contrast(
            dds_wald, contrast = list(c(ct, ci)), label = clabel)
        }
      }
      contrast_labels[cname] <- clabel
      contrast_short[cname]  <- cshort
    }
  }
}

# --- Progression: {cl}_{t2}_vs_{t1} (consecutive times) ---
if (isTRUE(comp$progression) && length(nonref_tps) >= 2) {
  cat("\n-- Time progression --\n")
  for (cl in cl_ids) {
    for (i in 2:length(nonref_tps)) {
      t1 <- nonref_tps[i-1]
      t2 <- nonref_tps[i]
      cname <- paste0(cl, "_", t2, "_vs_", t1)
      clabel <- paste0(cl, " ", t2, " vs ", t1)
      cshort <- paste0(cl_short[cl], tp_short[t2], "v", tp_short[t1])

      if (cl == ref_cl) {
        contrast_list[[cname]] <- extract_contrast(
          dds_wald, contrast = c("time", t2, t1), label = clabel)
      } else {
        ct1 <- coef_time[[t1]]; ci1 <- coef_int[[t1]]
        ct2 <- coef_time[[t2]]; ci2 <- coef_int[[t2]]
        if (!is.na(ci1) && ci1 != "" && !is.na(ci2) && ci2 != "") {
          contrast_list[[cname]] <- extract_contrast(
            dds_wald, contrast = list(c(ct2, ci2), c(ct1, ci1)), label = clabel)
        }
      }
      contrast_labels[cname] <- clabel
      contrast_short[cname]  <- cshort
    }
  }
}

# --- Between cell lines: {cl2}_vs_{cl1} at each time ---
if (isTRUE(comp$between_cell_lines)) {
  cat("\n-- Between cell lines --\n")
  for (tp in tp_ids) {
    if (tp == ref_tp) {
      cname <- paste0(nonref_cl, "_vs_", ref_cl, "_", tp)
    } else {
      cname <- paste0(nonref_cl, "_vs_", ref_cl, "_", tp)
    }
    clabel <- paste0(nonref_cl, " vs ", ref_cl, " @", tp)
    if (tp == "mock") tp_label <- "mock" else tp_label <- tp
    cshort <- paste0("B", tp_label)

    if (tp == ref_tp) {
      contrast_list[[cname]] <- extract_contrast(
        dds_wald, contrast = c("cell_line", nonref_cl, ref_cl), label = clabel)
    } else {
      ci <- coef_int[[tp]]
      if (!is.na(ci) && ci != "") {
        contrast_list[[cname]] <- extract_contrast(
          dds_wald, contrast = list(c(coef_cell_line, ci)), label = clabel)
      }
    }
    contrast_labels[cname] <- clabel
    contrast_short[cname]  <- cshort
  }
}

# --- Interactions: differential time response ---
if (isTRUE(comp$interactions)) {
  cat("\n-- Interaction (differential response) --\n")
  for (tp in nonref_tps) {
    cname <- paste0("interaction_", tp)
    clabel <- paste0("Interaction @", tp)
    cshort <- paste0("I", tp_short[tp])

    ci <- coef_int[[tp]]
    if (!is.na(ci) && ci != "") {
      contrast_list[[cname]] <- extract_contrast(
        dds_wald, name = ci, label = clabel)
    }
    contrast_labels[cname] <- clabel
    contrast_short[cname]  <- cshort
  }
}

# Convert labels/shorts to DESeq2-compatible format for results merging
contrast_pv_suffix <- sapply(names(contrast_short), function(cn) {
  gsub("_", ".", cn)
})

# --- 5. Merge all results ---
cat("\n=== Building combined results table ===\n")

# Start with LRT results
combined <- as.data.frame(res_lrt)
combined$gene_id <- rownames(combined)
colnames(combined)[colnames(combined) == "pvalue"] <- "lrt_pvalue"
colnames(combined)[colnames(combined) == "padj"] <- "lrt_padj"
colnames(combined)[colnames(combined) == "stat"] <- "lrt_stat"

keep_cols <- c("gene_id", "baseMean", "lrt_stat", "lrt_pvalue", "lrt_padj")
combined <- combined[, keep_cols, drop = FALSE]

# Add all contrast results
for (cname in names(contrast_list)) {
  cres <- as.data.frame(contrast_list[[cname]])
  cres$gene_id <- rownames(cres)
  colnames(cres)[colnames(cres) == "log2FoldChange"] <- paste0(cname, "_log2FC")
  colnames(cres)[colnames(cres) == "lfcSE"]        <- paste0(cname, "_lfcSE")
  colnames(cres)[colnames(cres) == "pvalue"]        <- paste0(cname, "_pvalue")
  colnames(cres)[colnames(cres) == "padj"]          <- paste0(cname, "_padj")

  merge_cols <- c("gene_id",
                  paste0(cname, "_log2FC"),
                  paste0(cname, "_lfcSE"),
                  paste0(cname, "_pvalue"),
                  paste0(cname, "_padj"))
  cres <- cres[, merge_cols, drop = FALSE]
  combined <- merge(combined, cres, by = "gene_id", all.x = TRUE)
}

# Add mock_is_DE flag
combined$mock_is_DE <- combined[[paste0(nonref_cl, "_vs_", ref_cl, "_", ref_tp, "_padj")]] < alpha_val
combined$mock_is_DE[is.na(combined$mock_is_DE)] <- FALSE

# Merge annotations (rename gene_id to match)
colnames(annotations)[1] <- "gene_id"
combined <- merge(combined, annotations, by = "gene_id", all.x = TRUE)

# Sort by LRT p-value
combined <- combined[order(combined$lrt_pvalue), ]

# Add overall significance flag
combined$lrt_signif <- combined$lrt_padj < alpha_val
combined$lrt_signif[is.na(combined$lrt_signif)] <- FALSE

# Write full combined table
write.table(combined,
            file = file.path(out_dir, "combined_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote combined_results.tsv (%d genes)\n", nrow(combined)))

# Write significant subsets
signif_lrt <- combined[combined$lrt_signif, ]
write.table(signif_lrt,
            file = file.path(out_dir, "signif_lrt.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote signif_lrt.tsv (%d genes with LRT padj < %.2g)\n",
            nrow(signif_lrt), alpha_val))

# Genes with LRT padj < threshold AND at least one contrast |log2FC| > 1
fc_cols <- grep("_log2FC$", names(combined), value = TRUE)
has_fc1 <- apply(combined[, fc_cols, drop = FALSE], 1, function(x) any(abs(x) > 1, na.rm = TRUE))
signif_fc1 <- combined[combined$lrt_signif & has_fc1, ]
write.table(signif_fc1,
            file = file.path(out_dir, "signif_lrt_foldchange1.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote signif_lrt_foldchange1.tsv (%d genes with padj < %.2g & |log2FC|>1)\n",
            nrow(signif_fc1), alpha_val))

# --- 6. Summary statistics ---
cat("\n=== Summary Statistics ===\n")
cat(sprintf("\nTotal genes tested: %d\n", nrow(combined)))
cat(sprintf("LRT significant (padj < %.2g): %d (%.1f%%)\n",
            alpha_val,
            sum(combined$lrt_signif),
            100 * sum(combined$lrt_signif) / nrow(combined)))
cat(sprintf("LRT + |log2FC| > 1: %d\n\n", nrow(signif_fc1)))

cat(sprintf("Significant genes per contrast (padj < %.2g):\n", alpha_val))
for (cname in names(contrast_list)) {
  padj_col <- paste0(cname, "_padj")
  n_sig <- sum(combined[[padj_col]] < alpha_val, na.rm = TRUE)
  n_up  <- sum(combined[[padj_col]] < alpha_val &
               combined[[paste0(cname, "_log2FC")]] > 0, na.rm = TRUE)
  n_dn  <- sum(combined[[padj_col]] < alpha_val &
               combined[[paste0(cname, "_log2FC")]] < 0, na.rm = TRUE)
  cat(sprintf("  %-30s %6d sig  (up: %5d, down: %5d)\n", cname, n_sig, n_up, n_dn))
}

cat(sprintf("\nBaseline differences (mock_is_DE): %d genes differ at mock (padj < %.2g)\n",
            sum(combined$mock_is_DE, na.rm = TRUE), alpha_val))

# --- 7. Visualizations ---
cat("\n=== Creating Plots ===\n")

# 7a. PCA plot
vsd <- vst(dds_wald, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = c("cell_line", "time"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_line, shape = time)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA: VST-transformed counts") +
  theme_minimal(base_size = 14)
ggsave(file.path(out_dir, "pca_plot.pdf"), pca_plot, width = 7, height = 5)
cat("Saved pca_plot.pdf\n")

# 7b. PCA colored by batch
pca_batch <- plotPCA(vsd, intgroup = "batch", returnData = TRUE)
pca_batch_plot <- ggplot(pca_batch, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA colored by batch") +
  theme_minimal(base_size = 14)
ggsave(file.path(out_dir, "pca_batch_plot.pdf"), pca_batch_plot, width = 7, height = 5)
cat("Saved pca_batch_plot.pdf\n")

# 7c. Heatmap of top 50 LRT-significant genes
top50 <- head(combined[order(combined$lrt_pvalue), ], 50)
if (nrow(top50) > 2) {
  top_ids <- top50$gene_id
  mat <- assay(vsd)[top_ids, , drop = FALSE]
  rownames(mat) <- top50$gene_symbol[match(top_ids, top50$gene_id)]

  # Z-score normalize rows
  mat <- t(scale(t(mat)))

  # Column annotation
  annotation_col <- metadata[, c("cell_line", "time", "batch"), drop = FALSE]

  # Colors (auto-generated from design)
  cl_colors <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
                          "#F781BF", "#999999", "#E7298A", "#66A61E", "#E6AB02")[seq_along(cl_ids)], cl_ids)
  tp_colors <- setNames(c("#4DAF4A", "#FF7F00", "#984EA3", "#377EB8", "#E41A1C",
                          "#A65628", "#F781BF", "#999999")[seq_along(tp_ids)], tp_ids)
  batch_levels <- levels(metadata$batch)
  bt_colors <- setNames(c("#A65628", "#F781BF", "#999999", "#E7298A",
                          "#66A61E", "#E6AB02")[seq_along(batch_levels)], batch_levels)
  ann_colors <- list(
    cell_line = cl_colors,
    time      = tp_colors,
    batch     = bt_colors
  )

  pdf(file.path(out_dir, "heatmap_top50.pdf"), width = 10, height = 12)
  pheatmap(mat,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 6,
           fontsize_col = 8,
           main = "Top 50 LRT-significant genes (Z-score)")
  dev.off()
  cat("Saved heatmap_top50.pdf\n")
}

# 7d. Volcano plots for key contrasts
# Pick volcano contrasts: first two within-cell-line for each cell, plus baseline and interaction
volcano_contrasts <- names(contrast_list)
# Prefer: within_cell_line first, then between, then interactions
# Limit to ~6 key contrasts to avoid too many plots
volcano_contrasts <- intersect(
  c(grep(paste0("^", ref_cl, "_", nonref_tps[1]), names(contrast_list), value = TRUE),
    grep(paste0("^", ref_cl, "_", tail(nonref_tps, 1)), names(contrast_list), value = TRUE),
    grep(paste0("^", nonref_cl, "_", nonref_tps[1]), names(contrast_list), value = TRUE),
    grep(paste0("^", nonref_cl, "_", tail(nonref_tps, 1)), names(contrast_list), value = TRUE),
    grep(paste0(nonref_cl, "_vs_", ref_cl, "_mock"), names(contrast_list), value = TRUE),
    grep("^interaction_", names(contrast_list), value = TRUE)[1]),
  names(contrast_list))
volcano_contrasts <- na.omit(volcano_contrasts)
if (length(volcano_contrasts) > 6) volcano_contrasts <- head(volcano_contrasts, 6)
for (vc in volcano_contrasts) {
  if (!vc %in% names(contrast_list)) next

  log2fc_col <- paste0(vc, "_log2FC")
  padj_col   <- paste0(vc, "_padj")

  vol_data <- data.frame(
    gene_id    = combined$gene_id,
    gene_symbol = combined$gene_symbol,
    log2FC     = combined[[log2fc_col]],
    padj       = combined[[padj_col]],
    stringsAsFactors = FALSE
  )
  vol_data <- vol_data[!is.na(vol_data$padj), ]
  vol_data$neg_log10_padj <- -log10(vol_data$padj)
  vol_data$signif <- vol_data$padj < alpha_val
  vol_data$label <- ifelse(vol_data$signif & abs(vol_data$log2FC) > 2,
                           vol_data$gene_symbol, "")

  p <- ggplot(vol_data, aes(x = log2FC, y = neg_log10_padj, color = signif)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 30,
                    box.padding = 0.3, point.padding = 0.2) +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
    geom_hline(yintercept = -log10(alpha_val), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
    xlab("log2 Fold Change") +
    ylab("-log10(adjusted p-value)") +
    ggtitle(paste0("Volcano: ", vc)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  ggsave(file.path(out_dir, paste0("volcano_", vc, ".pdf")), p, width = 7, height = 6)
}
cat("Saved volcano plots\n")

# 7e. Sample-to-sample distance heatmap
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(sample_dist_matrix) <- metadata$sample_id

pdf(file.path(out_dir, "sample_distance_heatmap.pdf"), width = 8, height = 7)
pheatmap(sample_dist_matrix,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         main = "Sample-to-sample distances (VST)")
dev.off()
cat("Saved sample_distance_heatmap.pdf\n")

# --- 8. Normalized counts output ---
# Write VST-normalized counts for external use
vst_counts <- assay(vsd)
vst_df <- as.data.frame(vst_counts)
vst_df$gene_id <- rownames(vst_df)
vst_df <- vst_df[, c("gene_id", metadata$sample_id)]
write.table(vst_df,
            file = file.path(out_dir, "vst_normalized_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved vst_normalized_counts.tsv\n")

cat("\n=== Analysis Complete ===\n")
cat(sprintf("All outputs in: %s/\n", out_dir))
