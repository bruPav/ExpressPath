#!/usr/bin/env Rscript
#
# DESeq2 analysis for RNA-seq time course experiment
# Reads experimental design from design.yaml
# Part A: Likelihood Ratio Test (overall time course effect)
# Part B: Pairwise Wald contrasts (auto-generated from design)
#

# ─── Snakemake integration ───
temporal_n_clusters <- 6
temporal_fc_thresh <- 1
if (exists("snakemake")) {
  out_dir <- dirname(snakemake@output[["combined"]])
  alpha_val <- 0.05
  if (!is.null(snakemake@params$temporal_n_clusters))
    temporal_n_clusters <- as.integer(snakemake@params$temporal_n_clusters)
  if (!is.null(snakemake@params$temporal_fc_thresh))
    temporal_fc_thresh <- as.numeric(snakemake@params$temporal_fc_thresh)
} else {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args) >= 1) args[1] else "results"
  alpha_val <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
  temporal_n_clusters <- if (length(args) >= 4) as.integer(args[4]) else 6
  temporal_fc_thresh   <- if (length(args) >= 5) as.numeric(args[5]) else 1
}

suppressPackageStartupMessages({
  library("DESeq2")
  library("ggplot2")
  library("pheatmap")
  library("RColorBrewer")
  library("dplyr")
  library("ggrepel")
  library("jsonlite")
  library("Mfuzz")
  library("ComplexUpset")
  library("VennDiagram")
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
tp_order   <- setNames(sapply(tp_list, function(x) x$temporal_order %||% which(tp_ids == x$id)), tp_ids)

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

# --- 9. Part C: Time-Series Clustering (Mfuzz) ---
cat("\n=== Part C: Time-Series Clustering (Mfuzz) ===\n")

tp_order_sorted <- sort(tp_order)
ordered_tps <- names(tp_order_sorted)
cat(sprintf("Temporal order: %s\n", paste(sprintf("%s(%d)", ordered_tps, tp_order_sorted), collapse = " -> ")))

lrt_signif_genes <- combined$gene_id[combined$lrt_signif]
cat(sprintf("LRT-significant genes: %d\n", length(lrt_signif_genes)))

cluster_assignments_list <- list()
cluster_profiles_list <- list()

if (length(lrt_signif_genes) >= 10) {
  n_clust <- if (is.character(temporal_n_clusters) && temporal_n_clusters == "auto") {
    min(c(6, floor(sqrt(length(lrt_signif_genes)))))
  } else {
    as.integer(temporal_n_clusters)
  }

  for (cl in cl_ids) {
    cl_samples <- rownames(metadata)[metadata$cell_line == cl]
    if (length(cl_samples) == 0) next

    tp_means <- list()
    for (tp in ordered_tps) {
      tp_samples <- cl_samples[metadata[cl_samples, "time"] == tp]
      if (length(tp_samples) == 0) next
      vst_sub <- assay(vsd)[lrt_signif_genes, tp_samples, drop = FALSE]
      tp_means[[tp]] <- rowMeans(vst_sub, na.rm = TRUE)
    }

    tp_matrix <- do.call(cbind, tp_means)
    colnames(tp_matrix) <- names(tp_means)

    if (nrow(tp_matrix) < 10) {
      cat(sprintf("  %s: too few genes (%d), skipping\n", cl, nrow(tp_matrix)))
      next
    }

    cat(sprintf("  Clustering %s (%d genes x %d timepoints)...\n", cl, nrow(tp_matrix), ncol(tp_matrix)))

    tmp_expr <- tryCatch({
      new("ExpressionSet", exprs = tp_matrix)
    }, error = function(e) NULL)
    if (is.null(tmp_expr)) next

    tmp_s <- standardise(tmp_expr)
    m1 <- mestimate(tmp_s)
    cl_result <- mfuzz(tmp_s, c = n_clust, m = m1)

    memb <- cl_result$membership
    colnames(memb) <- paste0("C", 1:ncol(memb))
    cluster_assign <- data.frame(
      gene_id = rownames(memb),
      cell_line = cl,
      stringsAsFactors = FALSE
    )
    cluster_assign <- cbind(cluster_assign, as.data.frame(memb))

    best_cluster <- max.col(memb)
    cluster_assign$cluster <- paste0("C", best_cluster)
    cluster_assign$membership_score <- apply(memb, 1, max)
    cluster_assignments_list[[cl]] <- cluster_assign

    # Mean profile per cluster
    for (cn in 1:ncol(memb)) {
      cluster_genes <- rownames(memb)[best_cluster == cn]
      if (length(cluster_genes) == 0) next
      mean_prof <- colMeans(tp_matrix[cluster_genes, , drop = FALSE], na.rm = TRUE)
      row <- data.frame(
        cell_line = cl,
        cluster = paste0("C", cn),
        n_genes = length(cluster_genes),
        t(mean_prof),
        stringsAsFactors = FALSE
      )
      cluster_profiles_list[[length(cluster_profiles_list) + 1]] <- row
    }

    pdf(file.path(out_dir, paste0("cluster_profiles_", cl, ".pdf")), width = 10, height = 8)
    mfuzz.plot2(tmp_s, cl = cl_result, mfrow = c(ceiling(n_clust / 3), min(3, n_clust)),
                time.labels = colnames(tp_matrix),
                xlab = "Time point", ylab = "Expression")
    dev.off()
    cat(sprintf("    Saved cluster_profiles_%s.pdf (%d clusters)\n", cl, n_clust))
  }

  cluster_assign_all <- do.call(rbind, cluster_assignments_list)
  if (!is.null(cluster_assign_all) && nrow(cluster_assign_all) > 0) {
    write.table(cluster_assign_all,
                file = file.path(out_dir, "cluster_assignments.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("Wrote cluster_assignments.tsv (%d genes)\n", nrow(cluster_assign_all)))
  }

  cluster_prof_all <- do.call(rbind, cluster_profiles_list)
  if (!is.null(cluster_prof_all) && nrow(cluster_prof_all) > 0) {
    write.table(cluster_prof_all,
                file = file.path(out_dir, "cluster_mean_profiles.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("Wrote cluster_mean_profiles.tsv (%d profiles)\n", nrow(cluster_prof_all)))
  }
} else {
  cat("Fewer than 10 LRT-significant genes; skipping clustering.\n")
  write.table(data.frame(gene_id = character(), cell_line = character(),
                         cluster = character(), membership_score = numeric(),
                         stringsAsFactors = FALSE),
              file = file.path(out_dir, "cluster_assignments.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(cell_line = character(), cluster = character(),
                         n_genes = integer(), stringsAsFactors = FALSE),
              file = file.path(out_dir, "cluster_mean_profiles.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# --- 10. Part D: Velocity of Response ---
cat("\n=== Part D: Velocity of Response ===\n")

velocity_rows <- list()
for (cl in cl_ids) {
  for (tp in nonref_tps) {
    cname <- paste0(cl, "_", tp, "_vs_mock")
    padj_col <- paste0(cname, "_padj")
    lfc_col  <- paste0(cname, "_log2FC")

    if (!padj_col %in% colnames(combined)) next

    is_sig <- combined[[padj_col]] < alpha_val & !is.na(combined[[padj_col]])
    n_sig <- sum(is_sig, na.rm = TRUE)
    n_up  <- sum(is_sig & combined[[lfc_col]] > 0, na.rm = TRUE)
    n_down <- sum(is_sig & combined[[lfc_col]] < 0, na.rm = TRUE)
    mean_abs_lfc <- if (n_sig > 0) mean(abs(combined[[lfc_col]][is_sig]), na.rm = TRUE) else 0

    velocity_rows[[length(velocity_rows) + 1]] <- data.frame(
      cell_line = cl,
      timepoint = tp,
      n_up = n_up,
      n_down = n_down,
      n_total = n_sig,
      mean_abs_log2FC = round(mean_abs_lfc, 4),
      stringsAsFactors = FALSE
    )
  }
}

velocity_summary <- do.call(rbind, velocity_rows)
write.table(velocity_summary,
            file = file.path(out_dir, "velocity_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote velocity_summary.tsv\n"))

if (nrow(velocity_summary) > 0) {
  p_vel <- ggplot(velocity_summary, aes(x = timepoint, y = n_total, fill = cell_line)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = n_total), position = position_dodge(0.7), vjust = -0.3, size = 3) +
    labs(x = "Time Point", y = "Number of DEGs",
         title = "Response Velocity: DEGs per Time Point",
         fill = "Cell Line") +
    theme_minimal(base_size = 14)
  ggsave(file.path(out_dir, "velocity_barplot.pdf"), p_vel, width = 8, height = 5)
  cat("Saved velocity_barplot.pdf\n")

  fc_data <- list()
  for (cl in cl_ids) {
    for (tp in nonref_tps) {
      cname <- paste0(cl, "_", tp, "_vs_mock")
      lfc_col  <- paste0(cname, "_log2FC")
      padj_col <- paste0(cname, "_padj")
      if (!padj_col %in% colnames(combined)) next
      is_sig <- combined[[padj_col]] < alpha_val & !is.na(combined[[padj_col]])
      if (sum(is_sig) > 0) {
        fc_data[[length(fc_data) + 1]] <- data.frame(
          cell_line = cl,
          timepoint = tp,
          log2FC = combined[[lfc_col]][is_sig],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(fc_data) > 0) {
    fc_df <- do.call(rbind, fc_data)
    fc_df$tp_cl <- paste(fc_df$timepoint, fc_df$cell_line, sep = "_")
    p_box <- ggplot(fc_df, aes(x = timepoint, y = log2FC, fill = cell_line)) +
      geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      labs(x = "Time Point", y = "log2 Fold Change",
           title = "Response Magnitude: log2FC Distribution per Time Point") +
      theme_minimal(base_size = 14)
    ggsave(file.path(out_dir, "velocity_fc_boxplot.pdf"), p_box, width = 8, height = 5)
    cat("Saved velocity_fc_boxplot.pdf\n")
  }
}

# --- 11. Part E: Persistence Check (Venn/UpSet + Heatmap) ---
cat("\n=== Part E: Persistence Check (Venn/UpSet + Heatmap) ===\n")

if (length(nonref_tps) >= 2) {
  venn_genelists <- list()
  persistence_rows <- list()
  gene_act_rows <- list()

  for (cl in cl_ids) {
    deg_sets <- list()
    lfc_vals <- list()
    padj_vals <- list()
    for (tp in nonref_tps) {
      cname <- paste0(cl, "_", tp, "_vs_mock")
      padj_col <- paste0(cname, "_padj")
      lfc_col  <- paste0(cname, "_log2FC")
      if (!padj_col %in% colnames(combined) || !lfc_col %in% colnames(combined)) next

      is_sig <- combined[[padj_col]] < alpha_val & !is.na(combined[[padj_col]])
      if (temporal_fc_thresh > 0) {
        is_sig <- is_sig & abs(combined[[lfc_col]]) >= temporal_fc_thresh
      }
      deg_genes <- combined$gene_id[is_sig]
      deg_sets[[tp]] <- deg_genes

      # Store log2FC and padj for these genes (named vector)
      tmp_lfc <- combined[[lfc_col]]
      tmp_padj <- combined[[padj_col]]
      names(tmp_lfc) <- names(tmp_padj) <- combined$gene_id
      lfc_vals[[tp]] <- tmp_lfc
      padj_vals[[tp]] <- tmp_padj

      for (g in deg_genes) {
        venn_genelists[[length(venn_genelists) + 1]] <- data.frame(
          cell_line = cl,
          timepoint = tp,
          gene_id = g,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(deg_sets) < 2) next

    all_deg_genes <- unique(unlist(deg_sets))
    for (gene in all_deg_genes) {
      tps_present <- names(deg_sets)[sapply(deg_sets, function(s) gene %in% s)]
      first_tp <- tps_present[1]
      last_tp  <- tps_present[length(tps_present)]
      n_tps    <- length(tps_present)

      if (n_tps == 1) {
        if (first_tp == nonref_tps[1]) {
          category <- "Transient"
        } else if (first_tp == tail(nonref_tps, 1)) {
          category <- "Secondary_Deferred"
        } else {
          category <- "Transient_Mid"
        }
      } else {
        if (first_tp == nonref_tps[1] && last_tp == tail(nonref_tps, 1)) {
          category <- "Sustained"
        } else if (first_tp == nonref_tps[1]) {
          category <- "Partially_Sustained"
        } else {
          category <- "Complex"
        }
      }

      persistence_rows[[length(persistence_rows) + 1]] <- data.frame(
        gene_id = gene,
        cell_line = cl,
        category = category,
        first_timepoint = first_tp,
        last_timepoint = last_tp,
        n_timepoints = n_tps,
        stringsAsFactors = FALSE
      )

      # Build gene activity row
      gsym <- combined$gene_symbol[match(gene, combined$gene_id)]
      row <- data.frame(gene_id = gene, gene_symbol = gsym, cell_line = cl,
                        category = category, stringsAsFactors = FALSE)
      for (tp in nonref_tps) {
        lfc_v <- if (gene %in% names(lfc_vals[[tp]])) lfc_vals[[tp]][gene] else NA
        padj_v <- if (gene %in% names(padj_vals[[tp]])) padj_vals[[tp]][gene] else NA
        is_s <- gene %in% deg_sets[[tp]]
        row[[paste0("sig_", tp)]] <- is_s
        row[[paste0("log2FC_", tp)]] <- round(lfc_v, 4)
        row[[paste0("padj_", tp)]] <- if (is.na(padj_v)) NA else round(padj_v, 6)
      }
      gene_act_rows[[length(gene_act_rows) + 1]] <- row
    }

    # --- Venn or UpSet plot ---
    n_sets <- length(deg_sets)
    if (n_sets <= 3) {
      # Venn diagram
      venn_sets <- list()
      for (tp in names(deg_sets)) venn_sets[[tp]] <- deg_sets[[tp]]

      pdf(file.path(out_dir, paste0("venn_plot_", cl, ".pdf")), width = 7, height = 7)
      if (n_sets == 2) {
        grid.newpage()
        draw.pairwise.venn(
          area1 = length(venn_sets[[1]]), area2 = length(venn_sets[[2]]),
          cross.area = length(intersect(venn_sets[[1]], venn_sets[[2]])),
          category = names(venn_sets),
          fill = c("#E41A1C", "#377EB8"), alpha = 0.5,
          cex = 1.5, cat.cex = 1.3, cat.pos = c(-30, 30),
          margin = 0.05
        )
      } else {
        a12 <- length(intersect(venn_sets[[1]], venn_sets[[2]]))
        a13 <- length(intersect(venn_sets[[1]], venn_sets[[3]]))
        a23 <- length(intersect(venn_sets[[2]], venn_sets[[3]]))
        a123 <- length(Reduce(intersect, venn_sets))
        grid.newpage()
        draw.triple.venn(
          area1 = length(venn_sets[[1]]), area2 = length(venn_sets[[2]]),
          area3 = length(venn_sets[[3]]),
          n12 = a12, n13 = a13, n23 = a23, n123 = a123,
          category = names(venn_sets),
          fill = c("#E41A1C", "#377EB8", "#4DAF4A"), alpha = 0.5,
          cex = 1.5, cat.cex = 1.3, margin = 0.05
        )
      }
      dev.off()
      png(file.path(out_dir, paste0("venn_plot_", cl, ".png")), width = 7, height = 7,
          units = "in", res = 150)
      if (n_sets == 2) {
        grid.newpage()
        draw.pairwise.venn(
          area1 = length(venn_sets[[1]]), area2 = length(venn_sets[[2]]),
          cross.area = length(intersect(venn_sets[[1]], venn_sets[[2]])),
          category = names(venn_sets),
          fill = c("#E41A1C", "#377EB8"), alpha = 0.5,
          cex = 1.5, cat.cex = 1.3, cat.pos = c(-30, 30),
          margin = 0.05
        )
      } else {
        grid.newpage()
        draw.triple.venn(
          area1 = length(venn_sets[[1]]), area2 = length(venn_sets[[2]]),
          area3 = length(venn_sets[[3]]),
          n12 = a12, n13 = a13, n23 = a23, n123 = a123,
          category = names(venn_sets),
          fill = c("#E41A1C", "#377EB8", "#4DAF4A"), alpha = 0.5,
          cex = 1.5, cat.cex = 1.3, margin = 0.05
        )
      }
      dev.off()
      cat(sprintf("Saved venn_plot_%s.pdf / .png\n", cl))
    } else {
      # UpSet plot for 4+ timepoints
      upset_genes <- unique(unlist(deg_sets))
      upset_matrix <- as.data.frame(sapply(names(deg_sets), function(x) {
        as.integer(upset_genes %in% deg_sets[[x]])
      }))
      colnames(upset_matrix) <- paste0(cl, "_", names(deg_sets))
      rownames(upset_matrix) <- upset_genes

      pdf(file.path(out_dir, paste0("upset_plot_", cl, ".pdf")), width = 10, height = 6)
      print(upset(upset_matrix, intersect = colnames(upset_matrix),
                  name = paste0("DEG Overlaps: ", cl),
                  width_ratio = 0.3))
      dev.off()
      png(file.path(out_dir, paste0("upset_plot_", cl, ".png")), width = 10, height = 6,
          units = "in", res = 150)
      print(upset(upset_matrix, intersect = colnames(upset_matrix),
                  name = paste0("DEG Overlaps: ", cl),
                  width_ratio = 0.3))
      dev.off()
      cat(sprintf("Saved upset_plot_%s.pdf / .png\n", cl))
    }

    # --- Gene activity heatmap ---
    if (length(all_deg_genes) >= 3) {
      cl_act_rows <- gene_act_rows[sapply(gene_act_rows, function(x) x$cell_line == cl)]
      cl_act_df <- do.call(rbind, cl_act_rows)
      cl_act_df <- cl_act_df[order(cl_act_df$category, cl_act_df$gene_id), ]

      # Build log2FC matrix
      lfc_cols <- grep("^log2FC_", names(cl_act_df), value = TRUE)
      lfc_mat <- as.matrix(cl_act_df[, lfc_cols, drop = FALSE])
      lfc_mat[is.na(lfc_mat)] <- 0
      rownames(lfc_mat) <- ifelse(is.na(cl_act_df$gene_symbol) | cl_act_df$gene_symbol == "--" | cl_act_df$gene_symbol == "",
                                  cl_act_df$gene_id, cl_act_df$gene_symbol)
      colnames(lfc_mat) <- sub("^log2FC_", "", lfc_cols)

      # Category annotation
      cat_colors <- c("Transient" = "#4DAF4A", "Sustained" = "#FF7F00",
                      "Secondary_Deferred" = "#377EB8", "Partially_Sustained" = "#984EA3",
                      "Transient_Mid" = "#F781BF", "Complex" = "#999999")
      ann_row <- data.frame(Category = cl_act_df$category, row.names = rownames(lfc_mat))
      ann_colors <- list(Category = cat_colors[intersect(names(cat_colors), unique(cl_act_df$category))])

      abs_vals <- abs(lfc_mat[is.finite(lfc_mat) & lfc_mat != 0])
      lim <- if (length(abs_vals) > 0) max(3, quantile(abs_vals, 0.90, na.rm = TRUE)) else 3
      lfc_mat_clamped <- lfc_mat
      lfc_mat_clamped[lfc_mat_clamped >  lim] <-  lim
      lfc_mat_clamped[lfc_mat_clamped < -lim] <- -lim

      pdf(file.path(out_dir, paste0("gene_activity_heatmap_", cl, ".pdf")),
          width = max(6, 2 + length(lfc_cols) * 1.2),
          height = max(6, nrow(lfc_mat) * 0.25))
      pheatmap(lfc_mat_clamped, annotation_row = ann_row, annotation_colors = ann_colors,
               cluster_rows = FALSE, cluster_cols = FALSE,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-lim, lim, length.out = 101),
               main = paste0("Gene Activity: ", cl),
               fontsize_row = 7, fontsize_col = 10,
               display_numbers = nrow(lfc_mat) <= 30,
               number_format = "%.2f", number_color = "black",
               border_color = NA, legend = TRUE)
      dev.off()
      png(file.path(out_dir, paste0("gene_activity_heatmap_", cl, ".png")),
          width = max(6, 2 + length(lfc_cols) * 1.2),
          height = max(6, nrow(lfc_mat) * 0.25),
          units = "in", res = 150)
      pheatmap(lfc_mat_clamped, annotation_row = ann_row, annotation_colors = ann_colors,
               cluster_rows = FALSE, cluster_cols = FALSE,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-lim, lim, length.out = 101),
               main = paste0("Gene Activity: ", cl),
               fontsize_row = 7, fontsize_col = 10,
               display_numbers = nrow(lfc_mat) <= 30,
               number_format = "%.2f", number_color = "black",
               border_color = NA, legend = TRUE)
      dev.off()
      cat(sprintf("Saved gene_activity_heatmap_%s.pdf / .png\n", cl))
    }
  }

  venn_df <- do.call(rbind, venn_genelists)
  write.table(venn_df,
              file = file.path(out_dir, "venn_genelists.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("Wrote venn_genelists.tsv (%d entries)\n", nrow(venn_df)))

  persist_df <- do.call(rbind, persistence_rows)
  write.table(persist_df,
              file = file.path(out_dir, "persistence_classes.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("Wrote persistence_classes.tsv (%d entries)\n", nrow(persist_df)))

  gene_act_df <- do.call(rbind, gene_act_rows)
  write.table(gene_act_df,
              file = file.path(out_dir, "gene_activity.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("Wrote gene_activity.tsv (%d entries)\n", nrow(gene_act_df)))

  if (!is.null(persist_df) && nrow(persist_df) > 0) {
    cat("\nPersistence classification summary:\n")
    for (cat_name in c("Transient", "Secondary_Deferred", "Sustained",
                       "Partially_Sustained", "Transient_Mid", "Complex")) {
      n <- sum(persist_df$category == cat_name)
      if (n > 0) cat(sprintf("  %-25s %6d\n", cat_name, n))
    }
  }
} else {
  cat("Only one non-reference timepoint; skipping Venn/persistence analysis.\n")
  write.table(data.frame(cell_line = character(), timepoint = character(),
                         gene_id = character(), stringsAsFactors = FALSE),
              file = file.path(out_dir, "venn_genelists.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(gene_id = character(), cell_line = character(),
                         category = character(), first_timepoint = character(),
                         last_timepoint = character(), n_timepoints = integer(),
                         stringsAsFactors = FALSE),
              file = file.path(out_dir, "persistence_classes.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(gene_id = character(), gene_symbol = character(),
                         cell_line = character(), category = character(),
                         stringsAsFactors = FALSE),
              file = file.path(out_dir, "gene_activity.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# --- 12. Part F: Cross-Cell-Line Comparison ---
cat("\n=== Part F: Cross-Cell-Line Comparison ===\n")

if (length(cl_ids) >= 2 && length(nonref_tps) >= 1) {

  # Step F1: Build DEG sets for all cell lines × timepoints from combined results
  all_deg_sets <- list()
  all_lfc      <- list()
  all_padj     <- list()

  for (cl in cl_ids) {
    all_deg_sets[[cl]] <- list()
    all_lfc[[cl]]      <- list()
    all_padj[[cl]]     <- list()
    for (tp in nonref_tps) {
      cname    <- paste0(cl, "_", tp, "_vs_mock")
      padj_col <- paste0(cname, "_padj")
      lfc_col  <- paste0(cname, "_log2FC")

      if (!padj_col %in% colnames(combined) || !lfc_col %in% colnames(combined)) next

      is_sig <- combined[[padj_col]] < alpha_val & !is.na(combined[[padj_col]])
      if (temporal_fc_thresh > 0) {
        is_sig <- is_sig & abs(combined[[lfc_col]]) >= temporal_fc_thresh
      }

      all_deg_sets[[cl]][[tp]] <- combined$gene_id[is_sig]

      lfc_v  <- combined[[lfc_col]];  names(lfc_v)  <- combined$gene_id
      padj_v <- combined[[padj_col]]; names(padj_v) <- combined$gene_id
      all_lfc[[cl]][[tp]]  <- lfc_v
      all_padj[[cl]][[tp]] <- padj_v
    }
  }

  # Step F2: Process each timepoint
  cross_shared_rows   <- list()
  cross_specific_rows <- list()

  for (tp in nonref_tps) {
    cat(sprintf("\n  Timepoint: %s\n", tp))

    cl_pairs <- combn(cl_ids, 2, simplify = FALSE)

    for (pair in cl_pairs) {
      cl_a <- pair[1]; cl_b <- pair[2]

      deg_a <- all_deg_sets[[cl_a]][[tp]]
      deg_b <- all_deg_sets[[cl_b]][[tp]]

      if (is.null(deg_a) || is.null(deg_b)) next
      if (length(deg_a) == 0 && length(deg_b) == 0) next

      shared_genes <- intersect(deg_a, deg_b)
      a_only <- setdiff(deg_a, deg_b)
      b_only <- setdiff(deg_b, deg_a)

      lfc_a_all  <- all_lfc[[cl_a]][[tp]]
      lfc_b_all  <- all_lfc[[cl_b]][[tp]]
      padj_a_all <- all_padj[[cl_a]][[tp]]
      padj_b_all <- all_padj[[cl_b]][[tp]]

      # --- Shared genes: concordance + magnitude classification ---
      for (g in shared_genes) {
        lfc_a  <- lfc_a_all[g]
        lfc_b  <- lfc_b_all[g]
        padj_a <- padj_a_all[g]
        padj_b <- padj_b_all[g]
        sym <- combined$gene_symbol[match(g, combined$gene_id)]

        if (!is.na(lfc_a) && !is.na(lfc_b)) {
          if (lfc_a > 0 && lfc_b > 0) {
            conc <- "Concordant_Up"
          } else if (lfc_a < 0 && lfc_b < 0) {
            conc <- "Concordant_Down"
          } else if (lfc_a > 0 && lfc_b < 0) {
            conc <- paste0(cl_a, "_Up_", cl_b, "_Down")
          } else if (lfc_a < 0 && lfc_b > 0) {
            conc <- paste0(cl_a, "_Down_", cl_b, "_Up")
          } else {
            conc <- "Zero"
          }
        } else {
          conc <- "NA"
        }

        abs_a <- abs(lfc_a); abs_b <- abs(lfc_b)
        if (!is.na(abs_a) && !is.na(abs_b) && is.finite(abs_a) && is.finite(abs_b) &&
            min(abs_a, abs_b) > 0) {
          mag_ratio <- max(abs_a, abs_b) / min(abs_a, abs_b)
        } else {
          mag_ratio <- NA
        }
        mag_div <- !is.na(mag_ratio) && mag_ratio > 2

        cross_shared_rows[[length(cross_shared_rows) + 1]] <- data.frame(
          gene_id = g, gene_symbol = sym, timepoint = tp,
          cell_line_1 = cl_a, cell_line_2 = cl_b,
          log2FC_1 = round(lfc_a, 4), padj_1 = round(padj_a, 6),
          log2FC_2 = round(lfc_b, 4), padj_2 = round(padj_b, 6),
          concordance = conc,
          magnitude_ratio = if (is.na(mag_ratio)) NA else round(mag_ratio, 3),
          magnitude_divergent = mag_div,
          stringsAsFactors = FALSE
        )
      }

      # --- Cell-line-specific genes ---
      for (g in a_only) {
        lfc  <- lfc_a_all[g]
        padj <- padj_a_all[g]
        sym  <- combined$gene_symbol[match(g, combined$gene_id)]
        dir  <- if (!is.na(lfc) && lfc > 0) "Up" else if (!is.na(lfc) && lfc < 0) "Down" else "Zero"
        cross_specific_rows[[length(cross_specific_rows) + 1]] <- data.frame(
          gene_id = g, gene_symbol = sym, timepoint = tp,
          cell_line = cl_a, log2FC = round(lfc, 4), padj = round(padj, 6),
          direction = dir, stringsAsFactors = FALSE
        )
      }
      for (g in b_only) {
        lfc  <- lfc_b_all[g]
        padj <- padj_b_all[g]
        sym  <- combined$gene_symbol[match(g, combined$gene_id)]
        dir  <- if (!is.na(lfc) && lfc > 0) "Up" else if (!is.na(lfc) && lfc < 0) "Down" else "Zero"
        cross_specific_rows[[length(cross_specific_rows) + 1]] <- data.frame(
          gene_id = g, gene_symbol = sym, timepoint = tp,
          cell_line = cl_b, log2FC = round(lfc, 4), padj = round(padj, 6),
          direction = dir, stringsAsFactors = FALSE
        )
      }

      # --- Console summary for this pair ---
      n_up   <- sum(sapply(shared_genes, function(g) {
        la <- lfc_a_all[g]; lb <- lfc_b_all[g]
        !is.na(la) && !is.na(lb) && la > 0 && lb > 0 }))
      n_down <- sum(sapply(shared_genes, function(g) {
        la <- lfc_a_all[g]; lb <- lfc_b_all[g]
        !is.na(la) && !is.na(lb) && la < 0 && lb < 0 }))
      n_disc <- length(shared_genes) - n_up - n_down
      n_mag  <- sum(sapply(shared_genes, function(g) {
        la <- lfc_a_all[g]; lb <- lfc_b_all[g]
        aa <- abs(la); ab <- abs(lb)
        if (is.na(aa) || is.na(ab) || !is.finite(aa) || !is.finite(ab) || min(aa, ab) <= 0)
          return(FALSE)
        max(aa, ab) / min(aa, ab) > 2
      }))

      cat(sprintf("    %s vs %s: %d shared (%d down, %d up, %d discordant",
                  cl_a, cl_b, length(shared_genes), n_down, n_up, n_disc))
      if (n_mag > 0) cat(sprintf("; %d magnitude-divergent", n_mag))
      cat(sprintf(")\n                %d %s-specific, %d %s-specific\n",
                  length(a_only), cl_a, length(b_only), cl_b))
    }
  }

  # Step F3: Write cross-cell-line output tables
  if (length(cross_shared_rows) > 0) {
    cross_shared_df <- do.call(rbind, cross_shared_rows)
    write.table(cross_shared_df,
                file = file.path(out_dir, "cross_cellline_shared.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("\nWrote cross_cellline_shared.tsv (%d rows)\n", nrow(cross_shared_df)))
  } else {
    write.table(data.frame(gene_id = character(), gene_symbol = character(),
                           timepoint = character(), cell_line_1 = character(),
                           cell_line_2 = character(), log2FC_1 = numeric(),
                           padj_1 = numeric(), log2FC_2 = numeric(),
                           padj_2 = numeric(), concordance = character(),
                           magnitude_ratio = numeric(), magnitude_divergent = logical(),
                           stringsAsFactors = FALSE),
                file = file.path(out_dir, "cross_cellline_shared.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Wrote cross_cellline_shared.tsv (empty)\n")
  }

  if (length(cross_specific_rows) > 0) {
    cross_specific_df <- do.call(rbind, cross_specific_rows)
    write.table(cross_specific_df,
                file = file.path(out_dir, "cross_cellline_specific.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("Wrote cross_cellline_specific.tsv (%d rows)\n", nrow(cross_specific_df)))
  } else {
    write.table(data.frame(gene_id = character(), gene_symbol = character(),
                           timepoint = character(), cell_line = character(),
                           log2FC = numeric(), padj = numeric(),
                           direction = character(), stringsAsFactors = FALSE),
                file = file.path(out_dir, "cross_cellline_specific.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Wrote cross_cellline_specific.tsv (empty)\n")
  }

  # Step F4: Cross-cell-line Venn + UpSet plots per timepoint
  for (tp in nonref_tps) {
    venn_sets <- lapply(cl_ids, function(cl) all_deg_sets[[cl]][[tp]])
    names(venn_sets) <- cl_ids
    venn_sets <- venn_sets[!sapply(venn_sets, is.null)]
    venn_sets <- venn_sets[sapply(venn_sets, length) > 0]

    if (length(venn_sets) < 2) next

    if (length(venn_sets) == 2) {
      cl_a <- names(venn_sets)[1]; cl_b <- names(venn_sets)[2]
      shared <- intersect(venn_sets[[1]], venn_sets[[2]])
      lfc_a <- all_lfc[[cl_a]][[tp]]; lfc_b <- all_lfc[[cl_b]][[tp]]
      n_up <- sum(sapply(shared, function(g) {
        la <- lfc_a[g]; lb <- lfc_b[g]
        !is.na(la) && !is.na(lb) && la > 0 && lb > 0 }))
      n_down <- sum(sapply(shared, function(g) {
        la <- lfc_a[g]; lb <- lfc_b[g]
        !is.na(la) && !is.na(lb) && la < 0 && lb < 0 }))

      pdf(file.path(out_dir, paste0("cross_venn_", tp, ".pdf")), width = 8, height = 7)
      grid.newpage()
      draw.pairwise.venn(
        area1 = length(venn_sets[[1]]), area2 = length(venn_sets[[2]]),
        cross.area = length(shared),
        category = c(cl_a, cl_b),
        fill = c("#E41A1C", "#377EB8"), alpha = 0.5,
        cex = 1.5, cat.cex = 1.3, cat.pos = c(-30, 30),
        margin = 0.05
      )
      dev.off()

      png(file.path(out_dir, paste0("cross_venn_", tp, ".png")),
          width = 8, height = 7, units = "in", res = 150)
      grid.newpage()
      draw.pairwise.venn(
        area1 = length(venn_sets[[1]]), area2 = length(venn_sets[[2]]),
        cross.area = length(shared),
        category = c(cl_a, cl_b),
        fill = c("#E41A1C", "#377EB8"), alpha = 0.5,
        cex = 1.5, cat.cex = 1.3, cat.pos = c(-30, 30),
        margin = 0.05
      )
      dev.off()

      cat(sprintf("Saved cross_venn_%s.pdf / .png\n", tp))
    } else {
      # UpSet for >2 cell lines
      all_deg_genes <- unique(unlist(venn_sets))
      upset_matrix <- as.data.frame(sapply(names(venn_sets), function(x) {
        as.integer(all_deg_genes %in% venn_sets[[x]])
      }))
      rownames(upset_matrix) <- all_deg_genes

      pdf(file.path(out_dir, paste0("cross_upset_", tp, ".pdf")), width = 10, height = 6)
      print(upset(upset_matrix, intersect = colnames(upset_matrix),
                  name = paste0("Cross-Cell-Line DEGs: ", tp),
                  width_ratio = 0.3))
      dev.off()
      png(file.path(out_dir, paste0("cross_upset_", tp, ".png")),
          width = 10, height = 6, units = "in", res = 150)
      print(upset(upset_matrix, intersect = colnames(upset_matrix),
                  name = paste0("Cross-Cell-Line DEGs: ", tp),
                  width_ratio = 0.3))
      dev.off()
      cat(sprintf("Saved cross_upset_%s.pdf / .png\n", tp))
    }
  }

  # Step F4b: Cross-cell-line log2FC scatter plots (one per timepoint)
  for (tp in nonref_tps) {
    cl_pairs <- combn(cl_ids, 2, simplify = FALSE)
    for (pair in cl_pairs) {
      cl_a <- pair[1]; cl_b <- pair[2]
      deg_a <- all_deg_sets[[cl_a]][[tp]]
      deg_b <- all_deg_sets[[cl_b]][[tp]]
      if (is.null(deg_a) || is.null(deg_b)) next
      if (length(deg_a) == 0 && length(deg_b) == 0) next

      lfc_a_all  <- all_lfc[[cl_a]][[tp]]
      lfc_b_all  <- all_lfc[[cl_b]][[tp]]
      padj_a_all <- all_padj[[cl_a]][[tp]]
      padj_b_all <- all_padj[[cl_b]][[tp]]

      # Build scatter data: all DEGs in either cell line
      all_degs <- union(deg_a, deg_b)
      scatter_df <- data.frame(
        gene_id = all_degs,
        log2FC_a = lfc_a_all[all_degs],
        log2FC_b = lfc_b_all[all_degs],
        stringsAsFactors = FALSE
      )
      scatter_df$gene_symbol <- combined$gene_symbol[match(scatter_df$gene_id, combined$gene_id)]
      scatter_df$shared  <- scatter_df$gene_id %in% intersect(deg_a, deg_b)
      scatter_df$a_only  <- scatter_df$gene_id %in% setdiff(deg_a, deg_b)
      scatter_df$b_only  <- scatter_df$gene_id %in% setdiff(deg_b, deg_a)

      # For shared genes, determine concordance
      scatter_df$concordance <- "Specific"
      scatter_df$magnitude_divergent <- FALSE
      scatter_df$label <- ""

      shared_ids <- intersect(deg_a, deg_b)
      for (g in shared_ids) {
        la <- lfc_a_all[g]; lb <- lfc_b_all[g]
        if (!is.na(la) && !is.na(lb)) {
          if (la > 0 && lb > 0) scatter_df$concordance[scatter_df$gene_id == g] <- "Up"
          else if (la < 0 && lb < 0) scatter_df$concordance[scatter_df$gene_id == g] <- "Down"
          else scatter_df$concordance[scatter_df$gene_id == g] <- "Discordant"
        }
        # Magnitude
        aa <- abs(la); ab <- abs(lb)
        if (!is.na(aa) && !is.na(ab) && is.finite(aa) && is.finite(ab) && min(aa, ab) > 0) {
          mr <- max(aa, ab) / min(aa, ab)
          scatter_df$magnitude_divergent[scatter_df$gene_id == g] <- (mr > 2)
        }
      }

      # Label magnitude-divergent shared genes
      scatter_df$label[scatter_df$magnitude_divergent] <- scatter_df$gene_symbol[scatter_df$magnitude_divergent]

      # Counts for legend
      n_shared_up   <- sum(scatter_df$concordance == "Up")
      n_shared_down <- sum(scatter_df$concordance == "Down")
      n_discordant  <- sum(scatter_df$concordance == "Discordant")
      n_specific_a  <- sum(scatter_df$a_only)
      n_specific_b  <- sum(scatter_df$b_only)

      # Axis limits: symmetric around 0
      max_abs <- max(abs(c(scatter_df$log2FC_a, scatter_df$log2FC_b)), na.rm = TRUE)
      lim <- max_abs * 1.15
      if (!is.finite(lim) || lim == 0) lim <- 5

      scatter_df$pt_color <- NA
      scatter_df$pt_color[scatter_df$concordance == "Up"]   <- "#E41A1C"
      scatter_df$pt_color[scatter_df$concordance == "Down"] <- "#377EB8"
      scatter_df$pt_color[scatter_df$concordance == "Discordant"] <- "#FF7F00"
      scatter_df$pt_color[scatter_df$concordance == "Specific"]   <- "grey70"
      scatter_df$pt_color[is.na(scatter_df$pt_color)] <- "grey70"

      scatter_df$pt_size <- ifelse(scatter_df$concordance == "Specific", 0.8, 1.5)
      scatter_df$pt_alpha <- ifelse(scatter_df$concordance == "Specific", 0.4, 0.7)

      scatter_df <- scatter_df[order(scatter_df$concordance == "Specific"), ]

      p <- ggplot(scatter_df, aes(x = log2FC_a, y = log2FC_b)) +
        geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
        geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
        geom_abline(slope = 1, intercept = 0, color = "grey40", linetype = "dashed", linewidth = 0.5) +
        geom_point(aes(color = pt_color), size = scatter_df$pt_size, alpha = scatter_df$pt_alpha) +
        scale_color_identity(
          guide = guide_legend(title = NULL),
          labels = c(
            "#E41A1C" = paste0("Both up (n=", n_shared_up, ")"),
            "#377EB8" = paste0("Both down (n=", n_shared_down, ")"),
            "#FF7F00" = paste0("Discordant (n=", n_discordant, ")"),
            "grey70" = paste0(cl_a, "-specific (", n_specific_a, ") / ", cl_b, "-specific (", n_specific_b, ")")
          ),
          breaks = c("#E41A1C", "#377EB8", "#FF7F00", "grey70")
        ) +
        geom_text_repel(aes(label = label), size = 3, max.overlaps = 25,
                        box.padding = 0.5, point.padding = 0.3, na.rm = TRUE) +
        coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
        labs(
          x = paste0(cl_a, " log2FC"), y = paste0(cl_b, " log2FC"),
          title = paste0("Cross-Cell-Line DE: ", tp),
          subtitle = paste0("Shared: ", n_shared_up + n_shared_down + n_discordant,
                            " genes (", n_shared_up, " up, ", n_shared_down, " down",
                            if (n_discordant > 0) paste0(", ", n_discordant, " discordant") else "",
                            ")")
        ) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "bottom",
              panel.grid.minor = element_blank())

      pdf(file.path(out_dir, paste0("cross_scatter_", tp, ".pdf")), width = 8, height = 7.5)
      print(p)
      dev.off()
      png(file.path(out_dir, paste0("cross_scatter_", tp, ".png")),
          width = 8, height = 7.5, units = "in", res = 150)
      print(p)
      dev.off()
      cat(sprintf("Saved cross_scatter_%s.pdf / .png\n", tp))
    }
  }

  # Step F5: Cross-cell-line persistence table (from existing persist_df)
  if (exists("persist_df") && !is.null(persist_df) && nrow(persist_df) > 0) {
    cross_persist <- persist_df[, c("gene_id", "cell_line", "category")]
    cross_persist$gene_symbol <- combined$gene_symbol[match(cross_persist$gene_id, combined$gene_id)]
    cross_persist <- cross_persist[, c("gene_id", "gene_symbol", "cell_line", "category")]

    write.table(cross_persist,
                file = file.path(out_dir, "cross_persistence.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("Wrote cross_persistence.tsv (%d rows)\n", nrow(cross_persist)))

    # Summary: concordance of persistence categories between cell lines
    cl_u <- unique(cross_persist$cell_line)
    if (length(cl_u) >= 2) {
      persist_wide <- reshape(cross_persist, idvar = c("gene_id", "gene_symbol"),
                              timevar = "cell_line", direction = "wide")
      cat_cols <- grep("^category\\.", colnames(persist_wide), value = TRUE)
      if (length(cat_cols) >= 2) {
        valid <- complete.cases(persist_wide[, cat_cols, drop = FALSE])
        n_same <- sum(persist_wide[valid, cat_cols[1]] == persist_wide[valid, cat_cols[2]])
        n_diff <- sum(valid) - n_same
        cat(sprintf("Cross persistence agreement: %d same, %d different (out of %d genes with categories in both)\n",
                    n_same, n_diff, sum(valid)))

        # --- Cross-persistence contingency heatmap ---
        if (sum(valid) >= 3) {
          tab <- table(
            factor(persist_wide[valid, cat_cols[1]],
                   levels = c("Transient", "Secondary_Deferred", "Sustained",
                              "Partially_Sustained", "Transient_Mid", "Complex")),
            factor(persist_wide[valid, cat_cols[2]],
                   levels = c("Transient", "Secondary_Deferred", "Sustained",
                              "Partially_Sustained", "Transient_Mid", "Complex"))
          )
          # Remove all-zero rows and columns
          tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]

          if (nrow(tab) >= 1 && ncol(tab) >= 1) {
            cl_a_name <- cl_u[1]; cl_b_name <- cl_u[2]
            rownames(tab) <- paste0(cl_a_name, ": ", rownames(tab))
            colnames(tab) <- paste0(cl_b_name, ": ", colnames(tab))

            # pheatmap with counts displayed
            pdf(file.path(out_dir, "cross_persistence_heatmap.pdf"),
                width = max(5, 1 + ncol(tab) * 1.8),
                height = max(3.5, 1 + nrow(tab) * 0.8))
            pheatmap(tab,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     display_numbers = TRUE, number_format = "%d",
                     number_color = "black",
                     color = colorRampPalette(c("white", "#4DAF4A", "#377EB8"))(100),
                     main = "Cross-Cell-Line Persistence Categories",
                     fontsize = 12, fontsize_number = 11,
                     legend = TRUE, legend_breaks = seq(0, max(tab), length.out = 5),
                     legend_labels = as.character(round(seq(0, max(tab), length.out = 5))),
                     border_color = "grey80")
            dev.off()
            png(file.path(out_dir, "cross_persistence_heatmap.png"),
                width = max(5, 1 + ncol(tab) * 1.8),
                height = max(3.5, 1 + nrow(tab) * 0.8),
                units = "in", res = 150)
            pheatmap(tab,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     display_numbers = TRUE, number_format = "%d",
                     number_color = "black",
                     color = colorRampPalette(c("white", "#4DAF4A", "#377EB8"))(100),
                     main = "Cross-Cell-Line Persistence Categories",
                     fontsize = 12, fontsize_number = 11,
                     legend = TRUE, legend_breaks = seq(0, max(tab), length.out = 5),
                     legend_labels = as.character(round(seq(0, max(tab), length.out = 5))),
                     border_color = "grey80")
            dev.off()
            cat("Saved cross_persistence_heatmap.pdf / .png\n")
          }
        }
      }
    }
  }

} else {
  cat("Fewer than 2 cell lines or no non-reference timepoints; skipping cross-cell-line analysis.\n")
  write.table(data.frame(gene_id = character(), gene_symbol = character(),
                         timepoint = character(), cell_line_1 = character(),
                         cell_line_2 = character(), log2FC_1 = numeric(),
                         padj_1 = numeric(), log2FC_2 = numeric(),
                         padj_2 = numeric(), concordance = character(),
                         magnitude_ratio = numeric(), magnitude_divergent = logical(),
                         stringsAsFactors = FALSE),
              file = file.path(out_dir, "cross_cellline_shared.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(gene_id = character(), gene_symbol = character(),
                         timepoint = character(), cell_line = character(),
                         log2FC = numeric(), padj = numeric(),
                         direction = character(), stringsAsFactors = FALSE),
              file = file.path(out_dir, "cross_cellline_specific.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(gene_id = character(), gene_symbol = character(),
                         cell_line = character(), category = character(),
                         stringsAsFactors = FALSE),
              file = file.path(out_dir, "cross_persistence.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("\n=== Analysis Complete ===\n")
cat(sprintf("All outputs in: %s/\n", out_dir))
