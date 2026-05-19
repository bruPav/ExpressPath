#!/usr/bin/env Rscript
#
# Pathway analysis: GSEA + Pathview + GSVA
# For all 11 contrasts from DESeq2 time-course analysis
#

# ─── Snakemake integration ───
if (exists("snakemake")) {
  results_dir <- dirname(dirname(snakemake@output[["gsea_kegg"]]))
  out_dir     <- dirname(snakemake@output[["gsea_kegg"]])
  gsea_min    <- 15
  gsea_max    <- 500
  design <- jsonlite::fromJSON(snakemake@params$design)
} else {
  args <- commandArgs(trailingOnly = TRUE)
  results_dir <- if (length(args) >= 1) args[1] else "results"
  out_dir     <- if (length(args) >= 2) args[2] else file.path(results_dir, "pathway")
  gsea_min    <- if (length(args) >= 3) as.integer(args[3]) else 15L
  gsea_max    <- if (length(args) >= 4) as.integer(args[4]) else 500L
  design_file <- if (length(args) >= 5) args[5] else "../data/design.yaml"
  if (file.exists(design_file)) {
    suppressPackageStartupMessages(library("yaml"))
    design <- yaml::read_yaml(design_file)
  } else {
    design <- list()
  }
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a) || length(a) == 0) b else a

suppressPackageStartupMessages({
  library("clusterProfiler")
  library("org.Hs.eg.db")
  library("pathview")
  library("GSVA")
  library("ggplot2")
  library("pheatmap")
  library("dplyr")
  library("msigdbr")
})

cat("=== Pathview Analysis: GSEA + Pathview + GSVA ===\n")
cat(sprintf("Results dir: %s\n", results_dir))
cat(sprintf("Pathway dir: %s\n", out_dir))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pathview_dir <- normalizePath(file.path(out_dir, "pathview_maps"), mustWork = FALSE)
dir.create(pathview_dir, showWarnings = FALSE, recursive = TRUE)
# --- Load DESeq2 results ---
cat("\nLoading DESeq2 combined results...\n")
combined <- read.delim(file.path(results_dir, "tables", "combined_results.tsv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("  %d genes loaded\n", nrow(combined)))

# --- ID Mapping: ENSG -> ENTREZ ---
cat("Running ENSG -> ENTREZ mapping...\n")
ensg_ids <- combined$gene_id
id_map <- bitr(ensg_ids, fromType = "ENSEMBL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
cat(sprintf("  %d of %d genes mapped to ENTREZ\n", nrow(id_map), length(ensg_ids)))

# Add ENTREZ to combined table
combined$entrez <- id_map$ENTREZID[match(combined$gene_id, id_map$ENSEMBL)]

# --- Define contrasts (read from combined_results column names) ---
combined_cols <- colnames(combined)
contrasts <- unique(gsub("_log2FC$", "",
  grep("_log2FC$", combined_cols, value = TRUE)))
cat(sprintf("  %d contrasts found\n", length(contrasts)))

# ==============================
# Part 1: GSEA (KEGG + GO) for all contrasts
# ==============================
cat("\n=== Part 1: GSEA ===\n")

set.seed(42)

gsea_kegg_all <- list()
gsea_go_all <- list()

for (cname in contrasts) {
  log2fc_col <- paste0(cname, "_log2FC")
  pval_col   <- paste0(cname, "_pvalue")

  if (!log2fc_col %in% names(combined)) {
    cat(sprintf("  SKIP %s: columns not found\n", cname))
    next
  }

  cat(sprintf("  %-30s ", cname))

  df <- combined[, c("gene_id", "entrez", log2fc_col, pval_col)]
  df <- df[!is.na(df$entrez) & !is.na(df[[log2fc_col]]) & !is.na(df[[pval_col]]), ]

  if (nrow(df) < 10) { cat("too few genes\n"); next }

  # Rank: -log10(pvalue) * sign(log2FC)
  df$rank_stat <- -log10(df[[pval_col]]) * sign(df[[log2fc_col]])
  df <- df[order(df$rank_stat, decreasing = TRUE), ]

  gene_list <- setNames(df$rank_stat, df$entrez)
  gene_list <- gene_list[!duplicated(names(gene_list))]

  # --- KEGG GSEA ---
  kegg_res <- tryCatch({
    gseKEGG(geneList = gene_list, organism = "hsa", pvalueCutoff = 0.05,
            minGSSize = gsea_min, maxGSSize = gsea_max, eps = 0, seed = 42)
  }, error = function(e) { cat("KEGG-err "); NULL })

  if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
    kr <- kegg_res@result
    kr$contrast <- cname
    gsea_kegg_all[[cname]] <- kr
    cat(sprintf("KEGG:%d ", sum(kr$p.adjust < 0.05)))
  } else {
    cat("KEGG:0 ")
  }

  # --- GO BP GSEA ---
  go_res <- tryCatch({
    gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, ont = "BP",
          pvalueCutoff = 0.05, minGSSize = gsea_min, maxGSSize = gsea_max,
          eps = 0, seed = 42)
  }, error = function(e) { cat("GO-err "); NULL })

  if (!is.null(go_res) && nrow(go_res@result) > 0) {
    gr <- go_res@result
    gr$contrast <- cname
    gsea_go_all[[cname]] <- gr
    cat(sprintf("GO:%d", sum(gr$p.adjust < 0.05)))
  } else {
    cat("GO:0")
  }
  cat("\n")
}

# --- Combine and save ---
cat("\nSaving GSEA results...\n")

if (length(gsea_kegg_all) > 0) {
  gsea_kegg_combined <- do.call(rbind, gsea_kegg_all)
  write.table(gsea_kegg_combined,
              file = file.path(out_dir, "gsea_kegg_all.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  gsea_kegg_all.tsv: %d rows\n", nrow(gsea_kegg_combined)))

  gsea_kegg_sig <- gsea_kegg_combined[gsea_kegg_combined$p.adjust < 0.05, ]
  write.table(gsea_kegg_sig,
              file = file.path(out_dir, "gsea_kegg_signif.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  gsea_kegg_signif.tsv: %d significant enrichments\n",
              nrow(gsea_kegg_sig)))
}

if (length(gsea_go_all) > 0) {
  gsea_go_combined <- do.call(rbind, gsea_go_all)
  write.table(gsea_go_combined,
              file = file.path(out_dir, "gsea_go_all.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  gsea_go_all.tsv: %d rows\n", nrow(gsea_go_combined)))

  gsea_go_sig <- gsea_go_combined[gsea_go_combined$p.adjust < 0.05, ]
  write.table(gsea_go_sig,
              file = file.path(out_dir, "gsea_go_signif.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  gsea_go_signif.tsv: %d significant enrichments\n",
              nrow(gsea_go_sig)))
}

# --- GSEA Dot Plot ---
if (length(gsea_kegg_all) > 0 && nrow(gsea_kegg_sig) > 2) {
  cat("\nGenerating GSEA dot plot...\n")
  plot_data <- gsea_kegg_sig %>%
    dplyr::group_by(contrast) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup()

  dp <- ggplot(plot_data, aes(x = NES, y = reorder(Description, NES),
                              size = -log10(p.adjust), color = p.adjust)) +
    geom_point() +
    facet_wrap(~ contrast, scales = "free_y", ncol = 3) +
    scale_color_gradient(low = "red", high = "blue") +
    theme_minimal(base_size = 9) +
    labs(x = "NES", y = "", title = "GSEA KEGG: Top 10 per contrast (padj < 0.05)") +
    theme(strip.text = element_text(size = 7),
          axis.text.y = element_text(size = 6))

  ggsave(file.path(out_dir, "gsea_dotplot_kegg.pdf"), dp,
         width = 16, height = max(8, nrow(gsea_kegg_sig) * 0.15),
         limitsize = FALSE)
  cat("  Saved gsea_dotplot_kegg.pdf\n")
}

# ==============================
# Part 2: Pathview Maps
# ==============================
cat("\n=== Part 2: Pathview ===\n")

# Select key contrasts for pathview (auto-pick: between-cell-line + first/last time contrasts)
pathview_contrasts <- intersect(c(
  grep("_vs_", contrasts, value = TRUE),            # between-cell-line contrasts
  grep("_vs_mock$", contrasts, value = TRUE),       # time vs mock
  grep("interaction_", contrasts, value = TRUE)     # interactions
), contrasts)
# Limit to max 7 contrasts (keep the most interesting ones)
if (length(pathview_contrasts) > 7) pathview_contrasts <- head(pathview_contrasts, 7)

# Get top KEGG pathways across all significant enrichments
if (exists("gsea_kegg_sig") && nrow(gsea_kegg_sig) > 0) {

  # Pick pathways with strong enrichments (lowest padj)
  top_kegg <- gsea_kegg_sig %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(min_padj = min(p.adjust), desc = dplyr::first(Description)) %>%
    dplyr::arrange(min_padj) %>%
    dplyr::slice_head(n = 20)

  top_kegg_ids <- top_kegg$ID

  # Create a root directory for pathview output images
  # (pathview writes final images to working directory, not kegg.dir)
  pv_root <- normalizePath(file.path(out_dir, "pathview_output"), mustWork = FALSE)
  dir.create(pv_root, showWarnings = FALSE, recursive = TRUE)
  main_wd <- getwd()

  for (cid in top_kegg_ids) {
    cid_desc <- top_kegg$desc[top_kegg$ID == cid]
    cid_desc <- gsub("[^A-Za-z0-9_-]", "_", cid_desc)
    cat(sprintf("  Pathview: %s (%s)\n", cid, cid_desc))

    # Sub-directory for this pathway
    cid_dir <- file.path(pv_root, paste0(cid, "_", cid_desc))
    dir.create(cid_dir, showWarnings = FALSE, recursive = TRUE)

    for (pc in pathview_contrasts) {
      log2fc_col <- paste0(pc, "_log2FC")
      if (!log2fc_col %in% names(combined)) next

      # Build ENTREZ-named log2FC vector
      fc_df <- combined[!is.na(combined$entrez) & !is.na(combined[[log2fc_col]]),
                        c("entrez", log2fc_col)]
      gene_fc <- setNames(fc_df[[log2fc_col]], as.character(fc_df$entrez))
      gene_fc <- gene_fc[!duplicated(names(gene_fc))]

      # Skip contrasts with very few meaningful fold changes
      n_meaningful <- sum(abs(gene_fc) > 0.3)
      if (n_meaningful < 5) next

      # Switch working directory so pathview writes images here
      setwd(cid_dir)

      out_suffix <- gsub("_", ".", pc)

      tryCatch({
        pathview(gene.data = gene_fc, pathway.id = cid, species = "hsa",
                 gene.idtype = "ENTREZ", kegg.dir = pathview_dir,
                 out.suffix = out_suffix, limit = list(gene = 2, cpd = 1),
                 low = list(gene = "blue", cpd = "blue"),
                 mid = list(gene = "grey", cpd = "grey"),
                 high = list(gene = "red", cpd = "yellow"),
                 kegg.native = TRUE, same.layer = FALSE)
      }, error = function(e) {
        cat(sprintf("    WARNING: pathview failed for %s/%s: %s\n",
                    cid, pc, e$message))
      })

      setwd(main_wd)
    }
  }

  # Clean up stray files in working directory
  stray_pngs <- list.files(main_wd, pattern = "^hsa.*\\.png$", full.names = TRUE)
  if (length(stray_pngs) > 0) file.remove(stray_pngs)
  stray_xmls <- list.files(main_wd, pattern = "^hsa.*\\.xml$", full.names = TRUE)
  if (length(stray_xmls) > 0) file.remove(stray_xmls)

  cat(sprintf("  Pathview maps saved to %s\n", pv_root))
} else {
  cat("  No significant KEGG enrichments for pathview\n")
}

# ==============================
# Part 3: GSVA Pathway Activity
# ==============================
cat("\n=== Part 3: GSVA ===\n")

# Load metadata
metadata <- read.delim(file.path(results_dir, "metadata.tsv"),
                       stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample_id
metadata$cell_line <- factor(metadata$cell_line)
metadata$time <- factor(metadata$time)
metadata$group <- factor(paste(metadata$cell_line, metadata$time, sep = "_"))

# Load VST counts
cat("Loading VST counts...\n")
vst_counts <- as.matrix(read.table(file.path(results_dir, "vst_normalized_counts.tsv"),
                                    header = TRUE, row.names = 1, sep = "\t",
                                    check.names = FALSE))

# Get Hallmark gene sets (ENTREZ)
cat("Loading Hallmark gene sets...\n")
msigdb_h <- msigdbr(species = "Homo sapiens", category = "H")
msigdb_h$entrez <- as.character(msigdb_h$entrez_gene)
h_gene_sets <- split(msigdb_h$entrez, msigdb_h$gs_name)
cat(sprintf("  %d Hallmark gene sets loaded\n", length(h_gene_sets)))

# Filter and map VST matrix to ENTREZ IDs
ensg_to_entrez <- setNames(combined$entrez, combined$gene_id)
ensg_to_entrez <- ensg_to_entrez[!is.na(ensg_to_entrez)]

common_genes <- intersect(rownames(vst_counts), names(ensg_to_entrez))
vst_entrez <- vst_counts[common_genes, , drop = FALSE]
rownames(vst_entrez) <- ensg_to_entrez[common_genes]

# Remove duplicate ENTREZ (keep first)
dup_entrez <- duplicated(rownames(vst_entrez))
vst_entrez <- vst_entrez[!dup_entrez, , drop = FALSE]

cat(sprintf("  Expression matrix: %d genes x %d samples\n",
            nrow(vst_entrez), ncol(vst_entrez)))

# Run GSVA (Hallmark)
cat("Running GSVA (Hallmark)...\n")
gsva_param <- gsvaParam(exprData = vst_entrez, geneSets = h_gene_sets,
                         minSize = gsea_min, maxSize = gsea_max, kcdf = "Gaussian")
gsva_res <- gsva(gsva_param)
cat(sprintf("  GSVA complete: %d gene sets x %d samples\n",
            nrow(gsva_res), ncol(gsva_res)))

# Save GSVA scores
gsva_df <- as.data.frame(gsva_res)
gsva_df$pathway <- rownames(gsva_df)
gsva_df <- gsva_df[, c("pathway", colnames(gsva_res))]
write.table(gsva_df, file = file.path(out_dir, "gsva_scores.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved gsva_scores.tsv\n")

# --- Differential Pathway Activity ---
cat("Testing differential pathway activity...\n")

# Auto-generate tests from metadata factors
cl_levs <- levels(metadata$cell_line)
tp_levs <- levels(metadata$time)
ref_cl <- design$factors$reference_cell_line %||% cl_levs[1]
nonref_cl <- setdiff(cl_levs, ref_cl)[1]
ref_tp <- design$factors$reference_time_point %||% tp_levs[1]
nonref_tps <- setdiff(tp_levs, ref_tp)

diff_tests <- list()
# Within-cell-line: each cell × each nonref time vs ref time
for (cl in cl_levs) {
  for (tp in nonref_tps) {
    tname <- paste0(cl, "_", tp, "_vs_", ref_tp)
    diff_tests[[tname]] <- list(
      samples_a = metadata$sample_id[metadata$cell_line == cl & metadata$time == tp],
      samples_b = metadata$sample_id[metadata$cell_line == cl & metadata$time == ref_tp]
    )
  }
}
# Between cell lines: at ref time and last nonref time
if (length(cl_levs) >= 2) {
  for (tp in c(ref_tp, tail(nonref_tps, 1))) {
    tname <- paste0(nonref_cl, "_vs_", ref_cl, "_", tp)
    diff_tests[[tname]] <- list(
      samples_a = metadata$sample_id[metadata$cell_line == nonref_cl & metadata$time == tp],
      samples_b = metadata$sample_id[metadata$cell_line == ref_cl & metadata$time == tp]
    )
  }
}

all_diff <- list()

for (test_name in names(diff_tests)) {
  ga <- diff_tests[[test_name]]$samples_a
  gb <- diff_tests[[test_name]]$samples_b

  if (length(ga) < 2 || length(gb) < 2) next

  for (i in seq_len(nrow(gsva_res))) {
    vals_a <- gsva_res[i, ga, drop = TRUE]
    vals_b <- gsva_res[i, gb, drop = TRUE]

    if (sd(c(vals_a, vals_b)) < 1e-6) next

    wt <- wilcox.test(vals_a, vals_b)
    tmean_diff <- mean(vals_a) - mean(vals_b)

    all_diff[[length(all_diff) + 1]] <- data.frame(
      test = test_name,
      pathway = rownames(gsva_res)[i],
      mean_diff = tmean_diff,
      pvalue = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
}

if (length(all_diff) > 0) {
  diff_table <- do.call(rbind, all_diff)
  diff_table$padj <- p.adjust(diff_table$pvalue, method = "BH")

  write.table(diff_table, file = file.path(out_dir, "gsva_diff_results.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  n_sig <- sum(diff_table$padj < 0.05)
  cat(sprintf("  gsva_diff_results.tsv: %d rows, %d significant (padj < 0.05)\n",
              nrow(diff_table), n_sig))

  # --- GSVA Heatmap ---
  sig_pathways <- diff_table$pathway[diff_table$padj < 0.05]
  if (length(sig_pathways) > 2) {
    sig_pathways <- unique(sig_pathways)
    sig_pathways <- head(sig_pathways, 50)

    mat <- gsva_res[sig_pathways, , drop = FALSE]
    rownames(mat) <- gsub("HALLMARK_", "", rownames(mat))

    annotation_col <- metadata[, c("cell_line", "time"), drop = FALSE]
    ann_colors <- list(
      cell_line = setNames(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628")[seq_along(cl_levs)], cl_levs),
      time = setNames(c("#4DAF4A","#FF7F00","#984EA3","#377EB8","#E41A1C","#A65628","#F781BF","#999999")[seq_along(tp_levs)], tp_levs)
    )

    pdf(file.path(out_dir, "gsva_heatmap.pdf"), width = 12, height = max(8, nrow(mat) * 0.3))
    pheatmap(mat, annotation_col = annotation_col, annotation_colors = ann_colors,
             scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
             show_rownames = TRUE, show_colnames = TRUE,
             fontsize_row = 7, fontsize_col = 8,
             main = "Significantly changed Hallmark pathways (GSVA)")
    dev.off()
    cat("  Saved gsva_heatmap.pdf\n")
  }
}

# --- GSVA Dot/Bar Plot ---
if (exists("diff_table") && nrow(diff_table) > 0) {
  top_diff <- diff_table %>%
    dplyr::group_by(test) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup()

  top_diff$neg_log10_padj <- -log10(top_diff$padj)
  top_diff$neg_log10_padj[is.infinite(top_diff$neg_log10_padj)] <- max(
    top_diff$neg_log10_padj[is.finite(top_diff$neg_log10_padj)], na.rm = TRUE)

  dbp <- ggplot(top_diff, aes(x = mean_diff, y = reorder(pathway, mean_diff),
                               fill = mean_diff, size = neg_log10_padj)) +
    geom_point(shape = 21) +
    facet_wrap(~ test, scales = "free_y", ncol = 2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme_minimal(base_size = 9) +
    labs(x = "Mean GSVA score difference", y = "",
         title = "Top 10 differentially active Hallmark pathways per test",
         fill = "Diff", size = "-log10(padj)") +
    theme(strip.text = element_text(size = 7),
          axis.text.y = element_text(size = 6))

  ggsave(file.path(out_dir, "gsva_diff_dotplot.pdf"), dbp,
         width = 14, height = max(8, nrow(top_diff) * 0.15),
         limitsize = FALSE)
  cat("  Saved gsva_diff_dotplot.pdf\n")
}

cat("\n=== Pathway Analysis Complete ===\n")
cat(sprintf("All outputs in: %s/\n", out_dir))
