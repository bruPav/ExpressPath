#!/usr/bin/env Rscript
#
# TF Target Enrichment via enrichR
# Reads cross-cell-line DEG outputs and persistence classes,
# runs TF target enrichment against multiple Enrichr libraries.
#
# Usage: Rscript 03b_tf_enrichment.R <results_dir>
#

# Self-healing: install CRAN-only packages if missing
for (pkg in c("enrichR", "visNetwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

suppressPackageStartupMessages({
  library("enrichR")
  library("ggplot2")
  library("dplyr")
})

# --- Argument handling ---
if (exists("snakemake")) {
  out_dir <- dirname(dirname(snakemake@output[["results"]]))
  dbs_str <- snakemake@params$databases
  enrichr_dbs <- strsplit(dbs_str, ",")[[1]]
} else {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args) >= 1) args[1] else "results"
  dbs_str <- if (length(args) >= 2) args[2] else
    "ChEA_2016,ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X,TRANSFAC_and_JASPAR_PWMs,ARCHS4_TFs_Coexp"
  enrichr_dbs <- strsplit(dbs_str, ",")[[1]]
}

cat("=== TF Target Enrichment (enrichR) ===\n")
cat(sprintf("Output dir: %s\n", out_dir))
cat(sprintf("Databases: %s\n", paste(enrichr_dbs, collapse = ", ")))

# Create output directories
dir.create(file.path(out_dir, "tf"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tf", "cross_cellline"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tf", "cross_temporal"), showWarnings = FALSE, recursive = TRUE)
min_genes <- 5       # minimum valid symbols per gene set
top_n_tfs <- 5       # top TFs per gene set for dotplot

# --- 1. Load data ---
cross_shared <- read.table(file.path(out_dir, "cross_temporal", "cross_cellline_shared.tsv"),
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
persist_df <- read.table(file.path(out_dir, "cross_temporal", "persistence_classes.tsv"),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
annot <- read.table(file.path(out_dir, "tables", "gene_annotations.tsv"),
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    quote = "", fill = TRUE, comment.char = "")

# Cross-temporal persistence (Step G)
cross_tpersist_file <- file.path(out_dir, "cross_cellline", "cross_temporal_persistence.tsv")
if (file.exists(cross_tpersist_file)) {
  cross_tpersist <- read.table(cross_tpersist_file,
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat(sprintf("Loaded cross_tpersist: %d rows\n", nrow(cross_tpersist)))
} else {
  cross_tpersist <- data.frame(gene_id = character(), pair = character(),
                                category = character(), n_timepoints = integer(),
                                stringsAsFactors = FALSE)
}

colnames(annot)[1] <- "gene_id"
cat(sprintf("Loaded: cross_shared=%d, persist=%d, cross_tpersist=%d, annotations=%d\n",
            nrow(cross_shared), nrow(persist_df), nrow(cross_tpersist), nrow(annot)))

# --- 2. Build gene_id → gene_symbol mapping ---
id2sym <- setNames(annot$gene_symbol, annot$gene_id)

is_valid_symbol <- function(s) {
  !is.na(s) & s != "" & s != "--" & !grepl("^NewGene", s)
}

# --- 3. Build gene sets ---
gene_sets <- list()

# 3a. Shared concordant up / down (from cross_cellline_shared)
if (nrow(cross_shared) > 0) {
  shared_up_genes <- unique(cross_shared$gene_id[cross_shared$concordance == "Concordant_Up"])
  shared_down_genes <- unique(cross_shared$gene_id[cross_shared$concordance == "Concordant_Down"])

  gene_sets[["Shared_Up"]]   <- shared_up_genes
  gene_sets[["Shared_Down"]] <- shared_down_genes
}

# 3b. Per cell_line × persistence category
if (nrow(persist_df) > 0) {
  for (cl in unique(persist_df$cell_line)) {
    for (cat_val in unique(persist_df$category)) {
      set_name <- paste0(cl, "_", cat_val)
      gene_sets[[set_name]] <- persist_df$gene_id[
        persist_df$cell_line == cl & persist_df$category == cat_val
      ]
    }
  }
}

# 3c. Cross-Sustained: genes where ALL cell lines have "Sustained"
if (nrow(persist_df) > 0) {
  cl_list <- unique(persist_df$cell_line)
  if (length(cl_list) >= 2) {
    sustained_per_cl <- lapply(cl_list, function(cl) {
      persist_df$gene_id[persist_df$cell_line == cl & persist_df$category == "Sustained"]
    })
    cross_sustained <- Reduce(intersect, sustained_per_cl)
    gene_sets[["Cross_Sustained"]] <- cross_sustained
  }
}

# 3d. Cross-temporal persistence categories (Step G)
if (nrow(cross_tpersist) > 0) {
  for (cat_val in unique(cross_tpersist$category)) {
    set_name <- paste0("Divergence_", cat_val)
    gene_sets[[set_name]] <- cross_tpersist$gene_id[cross_tpersist$category == cat_val]
  }

  # Cross-Emergent: all emergent patterns combined
  emergent_genes <- cross_tpersist$gene_id[
    grepl("^Emergent", cross_tpersist$category)
  ]
  if (length(emergent_genes) > 0) {
    gene_sets[["Divergence_Emergent_All"]] <- emergent_genes
  }

  # Divergence up/down: emergent genes with known direction
  cross_ga_file <- file.path(out_dir, "cross_cellline", "cross_temporal_gene_activity.tsv")
  if (file.exists(cross_ga_file)) {
    cross_ga <- read.table(cross_ga_file,
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    for (dset in c("Divergence_Emergent_All", "Divergence_Constitutive")) {
      orig <- if (dset == "Divergence_Emergent_All") "Emergent" else "Constitutive"
      ga_sub <- cross_ga[grepl(paste0("^", orig), cross_ga$category), ]
      if (nrow(ga_sub) > 0) {
        # Find the treatment timepoint with the largest abs log2FC to determine direction
        lfc_cols <- grep("^log2FC_", names(ga_sub), value = TRUE)
        if (length(lfc_cols) > 0) {
          ga_sub$max_lfc <- apply(ga_sub[, lfc_cols, drop = FALSE], 1,
                                   function(x) {
                                     abs_x <- abs(as.numeric(x))
                                     if (all(is.na(abs_x))) return(NA)
                                     as.numeric(x)[which.max(abs_x)]
                                   })
          up_genes   <- ga_sub$gene_id[ga_sub$max_lfc > 0 & !is.na(ga_sub$max_lfc)]
          down_genes <- ga_sub$gene_id[ga_sub$max_lfc < 0 & !is.na(ga_sub$max_lfc)]
          if (length(up_genes) > 0) {
            gene_sets[[paste0(dset, "_Up")]] <- up_genes
          }
          if (length(down_genes) > 0) {
            gene_sets[[paste0(dset, "_Down")]] <- down_genes
          }
        }
      }
    }
  }

  # 3e. Per-timepoint divergence gene sets (Divmock, Div1h, Div3h)
  if (file.exists(cross_ga_file) && nrow(cross_ga) > 0) {
    cross_ga <- read.table(cross_ga_file,
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    sig_cols <- grep("^sig_", names(cross_ga), value = TRUE)
    for (sc in sig_cols) {
      tp_name <- sub("^sig_", "", sc)
      lfc_col <- paste0("log2FC_", tp_name)
      is_sig <- isTRUE(cross_ga[[sc]]) | cross_ga[[sc]] == TRUE
      if (!all(is.na(is_sig))) is_sig[is.na(is_sig)] <- FALSE
      if (sum(is_sig) == 0) next

      # Short tag for clean gene set names
      tag <- gsub("^mock$", "mock", tp_name)
      tag <- gsub("h$", "h", tag)
      gene_sets[[paste0("Div", tag)]] <- cross_ga$gene_id[is_sig]

      # Directional splits
      if (lfc_col %in% names(cross_ga)) {
        up_genes <- cross_ga$gene_id[is_sig & cross_ga[[lfc_col]] > 0 & !is.na(cross_ga[[lfc_col]])]
        down_genes <- cross_ga$gene_id[is_sig & cross_ga[[lfc_col]] < 0 & !is.na(cross_ga[[lfc_col]])]
        if (length(up_genes) > 0) {
          gene_sets[[paste0("Div", tag, "_Up")]] <- up_genes
        }
        if (length(down_genes) > 0) {
          gene_sets[[paste0("Div", tag, "_Down")]] <- down_genes
        }
      }
    }
  }
}

# --- 4. Map to gene symbols and filter ---
cat("\nGene sets:\n")
gene_set_lists <- list()
for (nm in names(gene_sets)) {
  ids <- gene_sets[[nm]]
  syms <- id2sym[ids]
  names(syms) <- NULL
  valid <- syms[is_valid_symbol(syms)]

  n_total  <- length(ids)
  n_valid  <- length(valid)
  n_unique <- length(unique(valid))

  cat(sprintf("  %-25s %3d IDs  -> %3d valid symbols (%d unique)",
              nm, n_total, n_valid, n_unique))

  if (n_unique >= min_genes) {
    gene_set_lists[[nm]] <- unique(valid)
    cat(sprintf("  -> RUN\n"))
  } else {
    cat(sprintf("  -> SKIP (< %d valid symbols)\n", min_genes))
  }
}

if (length(gene_set_lists) == 0) {
  cat("WARNING: No gene sets with enough valid symbols. Writing empty output.\n")
  empty <- data.frame(gene_set = character(), library = character(), TF = character(),
                      n_overlap = integer(), pvalue = numeric(), adj_pvalue = numeric(),
                      odds_ratio = numeric(), combined_score = numeric(), genes = character(),
                      stringsAsFactors = FALSE)
  write.table(empty, file = file.path(out_dir, "tf", "tf_enrichment_results.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  quit(save = "no", status = 0)
}

# --- 5. Set Enrichr organism ---
listEnrichrSites()
Sys.sleep(0.5)  # rate-limit courtesy

# --- 6. Run enrichment ---
all_results <- list()

for (set_name in names(gene_set_lists)) {
  gene_list <- gene_set_lists[[set_name]]
  cat(sprintf("\nEnriching: %s (%d genes)...\n", set_name, length(gene_list)))

  enriched <- tryCatch({
    enrichr(gene_list, enrichr_dbs)
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    return(NULL)
  })

  if (is.null(enriched)) next

  for (db_name in names(enriched)) {
    res <- enriched[[db_name]]
    if (is.null(res) || nrow(res) == 0) next

    res$gene_set  <- set_name
    res$library   <- db_name
    res$n_overlap <- as.integer(sapply(strsplit(res$Overlap, "/"), `[`, 1))

    # Standardize column names
    colnames(res)[colnames(res) == "Term"] <- "TF"
    colnames(res)[colnames(res) == "P.value"] <- "pvalue"
    colnames(res)[colnames(res) == "Adjusted.P.value"] <- "adj_pvalue"
    colnames(res)[colnames(res) == "Odds.Ratio"] <- "odds_ratio"
    colnames(res)[colnames(res) == "Combined.Score"] <- "combined_score"
    colnames(res)[colnames(res) == "Genes"] <- "genes"

    keep_cols <- c("gene_set", "library", "TF", "n_overlap", "pvalue",
                   "adj_pvalue", "odds_ratio", "combined_score", "genes")
    res <- res[, intersect(keep_cols, colnames(res)), drop = FALSE]
    all_results[[length(all_results) + 1]] <- res
  }
  Sys.sleep(1)  # rate-limit between queries
}

if (length(all_results) == 0) {
  cat("WARNING: No enrichment results returned. Writing empty output.\n")
  empty <- data.frame(gene_set = character(), library = character(), TF = character(),
                      n_overlap = integer(), pvalue = numeric(), adj_pvalue = numeric(),
                      odds_ratio = numeric(), combined_score = numeric(), genes = character(),
                      stringsAsFactors = FALSE)
  write.table(empty, file = file.path(out_dir, "tf", "tf_enrichment_results.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  quit(save = "no", status = 0)
}

# --- 7. Combine and apply global FDR ---
results_df <- do.call(rbind, all_results)
cat(sprintf("\nTotal enrichment results: %d rows\n", nrow(results_df)))

results_df <- results_df[!is.na(results_df$pvalue), ]
results_df$adj_pvalue_global <- p.adjust(results_df$pvalue, method = "BH")
results_df <- results_df[order(results_df$adj_pvalue_global), ]

# Write TSV
write.table(results_df, file = file.path(out_dir, "tf", "tf_enrichment_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote tf_enrichment_results.tsv (%d rows)\n", nrow(results_df)))

# --- 7b. TF Activity Inference: cross-reference TFs with their own expression ---
cat("\n=== TF Activity Inference ===\n")

combined_expr <- read.table(file.path(out_dir, "tables", "combined_results.tsv"),
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                            quote = "", fill = TRUE, comment.char = "")

# Dynamically discover all {cell_line}_{timepoint}_vs_mock comparisons
lfc_cols  <- grep("_vs_mock_log2FC$", colnames(combined_expr), value = TRUE)
padj_cols <- grep("_vs_mock_padj$",  colnames(combined_expr), value = TRUE)
cat(sprintf("Found %d comparisons in combined_results\n", length(lfc_cols)))

# Clean TF name → gene symbol (first word before any space/suffix)
all_tf_names <- unique(results_df$TF)
tf_gene_symbols <- sub(" .*", "", all_tf_names)

# Build TF name → cleaned symbol lookup
tf2sym <- setNames(tf_gene_symbols, all_tf_names)

# For each unique TF, look up its expression in combined_results
tf_expr_info <- data.frame(
  TF = all_tf_names,
  gene_symbol = tf_gene_symbols,
  TF_is_DEG  = FALSE,
  TF_best_log2FC  = NA_real_,
  TF_best_padj    = NA_real_,
  TF_best_comparison = "",
  TF_n_sig  = 0L,
  stringsAsFactors = FALSE
)

n_tf_found <- 0
for (i in seq_len(nrow(tf_expr_info))) {
  sym <- tf_expr_info$gene_symbol[i]
  row <- combined_expr[combined_expr$gene_symbol == sym, ]
  if (nrow(row) == 0) next

  n_tf_found <- n_tf_found + 1

  # Extract all padj values for _vs_mock comparisons
  padj_vals <- as.numeric(unlist(row[1, padj_cols]))
  lfc_vals  <- as.numeric(unlist(row[1, lfc_cols]))

  if (all(is.na(padj_vals))) next

  best_idx <- which.min(padj_vals)
  if (length(best_idx) == 0 || is.na(padj_vals[best_idx])) next

  tf_expr_info$TF_is_DEG[i]       <- any(padj_vals < 0.05, na.rm = TRUE)
  tf_expr_info$TF_best_log2FC[i]  <- lfc_vals[best_idx]
  tf_expr_info$TF_best_padj[i]    <- padj_vals[best_idx]
  tf_expr_info$TF_best_comparison[i] <- sub("_padj$", "", padj_cols[best_idx])
  tf_expr_info$TF_n_sig[i]        <- sum(padj_vals < 0.05, na.rm = TRUE)
}

cat(sprintf("TF expression found for: %d / %d unique TFs (%d are DEGs in >=1 comparison)\n",
            n_tf_found, length(all_tf_names), sum(tf_expr_info$TF_is_DEG)))

# Join expression info back to enrichment results
results_df <- merge(results_df,
                    tf_expr_info[, c("TF", "TF_is_DEG", "TF_best_log2FC",
                                      "TF_best_padj", "TF_best_comparison", "TF_n_sig")],
                    by = "TF", all.x = TRUE)

# Re-write TSV with TF expression columns
write.table(results_df, file = file.path(out_dir, "tf", "tf_enrichment_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Re-wrote tf_enrichment_results.tsv with TF activity columns (%d rows)\n", nrow(results_df)))

# --- 8. Enrichment Heatmap: all significant TFs × all gene sets ---
suppressPackageStartupMessages({
  library("pheatmap")
  library("RColorBrewer")
})

# Filter to TFs with library-level adj_pvalue < 0.05 in at least one gene set
sig_tfs <- results_df[results_df$adj_pvalue < 0.05, ]
cat(sprintf("\nTFs with library-level adj.p < 0.05: %d rows from %d unique TFs\n",
            nrow(sig_tfs), length(unique(sig_tfs$TF))))

if (nrow(sig_tfs) > 0) {

  # Collapse across libraries: per TF × gene_set, keep best library by adj_pvalue
  sig_best <- sig_tfs %>%
    group_by(gene_set, TF) %>%
    slice_min(adj_pvalue, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Clean TF names and preserve full dataset before top-50 capping
  sig_best$TF_clean <- gsub(" [Hh]omo [Ss]apiens", "", sig_best$TF)
  sig_best$TF_clean <- gsub("_HUMAN$", "", sig_best$TF_clean)
  sig_best_full <- sig_best

  # Cap at top 50 TFs by minimum adj_pvalue
  tf_order_df <- sig_best %>%
    group_by(TF) %>%
    summarise(min_padj = min(adj_pvalue), .groups = "drop") %>%
    arrange(min_padj)
  top_tfs <- head(tf_order_df$TF, 50)
  if (nrow(tf_order_df) > 50) {
    cat(sprintf("Capping at top 50 TFs (out of %d)\n", nrow(tf_order_df)))
  }
  sig_best <- sig_best[sig_best$TF %in% top_tfs, ]
  sig_best$TF <- factor(sig_best$TF, levels = rev(top_tfs))

  # Build matrix: rows = TFs, cols = gene_sets, value = -log10(adj_pvalue)
  gs_levels <- sort(unique(sig_best$gene_set))
  tf_levels <- levels(sig_best$TF)

  heat_mat <- matrix(0, nrow = length(tf_levels), ncol = length(gs_levels))
  rownames(heat_mat) <- tf_levels
  colnames(heat_mat) <- gs_levels

  # Overlap count matrix for cell labels
  overlap_mat <- matrix("", nrow = length(tf_levels), ncol = length(gs_levels))
  rownames(overlap_mat) <- tf_levels
  colnames(overlap_mat) <- gs_levels

  for (i in seq_len(nrow(sig_best))) {
    tf_idx <- match(as.character(sig_best$TF[i]), tf_levels)
    gs_idx <- match(sig_best$gene_set[i], gs_levels)
    heat_mat[tf_idx, gs_idx] <- -log10(sig_best$adj_pvalue[i])
    overlap_mat[tf_idx, gs_idx] <- as.character(sig_best$n_overlap[i])
  }

  # Use clean TF names for the plot
  tf_clean_names <- sig_best$TF_clean[match(tf_levels, sig_best$TF)]
  tf_clean_names[is.na(tf_clean_names)] <- tf_levels[is.na(tf_clean_names)]
  rownames(heat_mat) <- tf_clean_names
  rownames(overlap_mat) <- tf_clean_names

  # Row annotation: which library each TF came from
  tf_lib <- sig_best %>%
    group_by(TF_clean) %>%
    slice_min(adj_pvalue, n = 1, with_ties = FALSE) %>%
    ungroup()
  tf_lib_vec <- setNames(tf_lib$library, tf_lib$TF_clean)
  tf_lib_vec <- tf_lib_vec[tf_clean_names]
  # Shorten library names
  short_lib <- c(
    "ChEA_2016" = "ChEA",
    "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X" = "ENCODE+ChEA",
    "TRANSFAC_and_JASPAR_PWMs" = "TRANSFAC/JASPAR",
    "ARCHS4_TFs_Coexp" = "ARCHS4"
  )
  # Row annotation: TF_is_DEG (is the TF itself differentially expressed?)
  tf_is_deg_vec <- sapply(tf_clean_names, function(nm) {
    orig_tf <- tf_levels[match(nm, tf_clean_names)]
    if (is.na(orig_tf)) return(FALSE)
    idx <- match(orig_tf, tf_expr_info$TF)
    if (is.na(idx)) return(FALSE)
    tf_expr_info$TF_is_DEG[idx]
  })

  ann_row <- data.frame(
    Library = short_lib[tf_lib_vec],
    TF_DEG   = ifelse(tf_is_deg_vec, "DEG", "Not_DEG"),
    row.names = tf_clean_names
  )
  lib_colors <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"),
                         c("ChEA", "ENCODE+ChEA", "TRANSFAC/JASPAR", "ARCHS4"))
  deg_colors <- setNames(c("#FF7F00", "grey80"), c("DEG", "Not_DEG"))

  # Column annotation: gene set type
  gs_type <- sapply(colnames(heat_mat), function(gs) {
    if (grepl("^Shared_", gs)) "Shared"
    else if (grepl("^Cross_", gs)) "Cross"
    else if (grepl("^Div", gs)) "Divergence"
    else "Per_CellLine"
  })
  ann_col <- data.frame(Type = gs_type, row.names = colnames(heat_mat))
  type_colors <- setNames(c("#FF7F00", "#377EB8", "#4DAF4A", "#E41A1C"),
                          c("Shared", "Per_CellLine", "Cross", "Divergence"))

  ann_colors_full <- list(
    Library = lib_colors[intersect(names(lib_colors), unique(ann_row$Library))],
    TF_DEG  = deg_colors[intersect(names(deg_colors), unique(ann_row$TF_DEG))],
    Type    = type_colors[intersect(names(type_colors), unique(ann_col$Type))]
  )

  # Color scale: white → green → orange → red
  max_val <- max(heat_mat, na.rm = TRUE)
  if (!is.finite(max_val) || max_val == 0) max_val <- 3
  breaks <- seq(0, ceiling(max_val), length.out = 101)
  heat_colors <- colorRampPalette(c("white", "#4DAF4A", "#FF7F00", "#E41A1C"))(100)

  # Plot dimensions
  plot_w <- max(6, 1 + ncol(heat_mat) * 1.3)
  plot_h <- max(5, 1 + nrow(heat_mat) * 0.35)

  pdf(file.path(out_dir, "tf", "tf_enrichment_heatmap.pdf"), width = plot_w, height = plot_h)
  pheatmap(heat_mat,
           annotation_row    = ann_row,
           annotation_col    = ann_col,
           annotation_colors = ann_colors_full,
           cluster_rows      = (nrow(heat_mat) > 2),
           cluster_cols      = (ncol(heat_mat) > 2),
           display_numbers   = overlap_mat,
           number_format     = "%s",
           number_color      = "black",
           fontsize_number   = 8,
           color             = heat_colors,
           breaks            = breaks,
           main              = "TF Target Enrichment\nColor = -log10(adj.p), Numbers = overlap count",
           fontsize          = 10,
           fontsize_row      = 9,
           fontsize_col      = 10,
           border_color      = "grey80",
           legend_breaks     = seq(0, ceiling(max_val), by = 1),
           name              = "-log10(adj.p)")
  dev.off()

  png(file.path(out_dir, "tf", "tf_enrichment_heatmap.png"),
      width = plot_w, height = plot_h, units = "in", res = 150)
  pheatmap(heat_mat,
           annotation_row    = ann_row,
           annotation_col    = ann_col,
           annotation_colors = ann_colors_full,
           cluster_rows      = (nrow(heat_mat) > 2),
           cluster_cols      = (ncol(heat_mat) > 2),
           display_numbers   = overlap_mat,
           number_format     = "%s",
           number_color      = "black",
           fontsize_number   = 8,
           color             = heat_colors,
           breaks            = breaks,
           main              = "TF Target Enrichment\nColor = -log10(adj.p), Numbers = overlap count",
           fontsize          = 10,
           fontsize_row      = 9,
           fontsize_col      = 10,
           border_color      = "grey80",
           legend_breaks     = seq(0, ceiling(max_val), by = 1),
           name              = "-log10(adj.p)")
  dev.off()

  png(file.path(out_dir, "tf", "tf_enrichment_heatmap.png"),
      width = plot_w, height = plot_h, units = "in", res = 150)
  pheatmap(heat_mat,
           annotation_row    = ann_row,
           annotation_col    = ann_col,
           annotation_colors = ann_colors_full,
           cluster_rows      = (nrow(heat_mat) > 2),
           cluster_cols      = (ncol(heat_mat) > 2),
           display_numbers   = overlap_mat,
           number_format     = "%s",
           number_color      = "black",
           fontsize_number   = 8,
           color             = heat_colors,
           breaks            = breaks,
           main              = "TF Target Enrichment\nColor = -log10(adj.p), Numbers = overlap count",
           fontsize          = 10,
           fontsize_row      = 9,
           fontsize_col      = 10,
           border_color      = "grey80",
           legend_breaks     = seq(0, ceiling(max_val), by = 1),
           name              = "-log10(adj.p)")
  dev.off()
  cat(sprintf("Saved tf_enrichment_heatmap.pdf / .png\n"))

  # --- 8c. Focused heatmaps per analysis dimension ---
  build_focused_heatmap <- function(subset_type, subdir, title, sig_data = sig_best_full) {
    all_gs <- unique(sig_data$gene_set)
    gs_allowed <- all_gs[vapply(all_gs, function(gs) {
      if (grepl("^Shared_", gs)) "Shared" %in% subset_type
      else if (grepl("^Cross_", gs)) "Cross" %in% subset_type
      else if (grepl("^Div", gs)) "Divergence" %in% subset_type
      else "Per_CellLine" %in% subset_type
    }, logical(1))]
    if (length(gs_allowed) == 0) return(FALSE)
    sb <- sig_data[sig_data$gene_set %in% gs_allowed, ]
    if (nrow(sb) == 0) return(FALSE)

    sb$gene_set <- factor(sb$gene_set, levels = sort(unique(as.character(sb$gene_set))))
    tf_focus <- sb %>%
      group_by(TF) %>%
      summarise(min_padj = min(adj_pvalue), .groups = "drop") %>%
      arrange(min_padj)
    top_focus <- head(tf_focus$TF, min(50, nrow(tf_focus)))
    sb <- sb[sb$TF %in% top_focus, ]
    sb$TF <- factor(sb$TF, levels = rev(top_focus))

    gs_focus <- sort(unique(as.character(sb$gene_set)))
    tf_focus_levels <- levels(sb$TF)

    fmat <- matrix(0, nrow = length(tf_focus_levels), ncol = length(gs_focus))
    rownames(fmat) <- tf_focus_levels
    colnames(fmat) <- gs_focus
    omat <- matrix("", nrow = length(tf_focus_levels), ncol = length(gs_focus))
    rownames(omat) <- tf_focus_levels
    colnames(omat) <- gs_focus

    for (i in seq_len(nrow(sb))) {
      ti <- match(as.character(sb$TF[i]), tf_focus_levels)
      gi <- match(sb$gene_set[i], gs_focus)
      if (!is.na(ti) && !is.na(gi)) {
        fmat[ti, gi] <- -log10(sb$adj_pvalue[i])
        omat[ti, gi] <- as.character(sb$n_overlap[i])
      }
    }

    tf_clean_f <- sb$TF_clean[match(tf_focus_levels, sb$TF)]
    tf_clean_f[is.na(tf_clean_f)] <- tf_focus_levels[is.na(tf_clean_f)]
    rownames(fmat) <- tf_clean_f
    rownames(omat) <- tf_clean_f

    fmax <- max(fmat, na.rm = TRUE)
    if (!is.finite(fmax) || fmax == 0) fmax <- 3

    fw <- min(14, 3 + length(gs_focus) * 1.2)
    fh <- max(5, 2 + length(tf_focus_levels) * 0.3)

    fbreaks <- seq(0, fmax + 0.5, length.out = 100)
    fcolors <- colorRampPalette(c("white", "#4DAF4A", "#FF7F00", "#E41A1C"))(100)

    pdf(file.path(out_dir, "tf", subdir, "tf_enrichment_heatmap.pdf"),
        width = fw, height = fh)
    pheatmap(fmat,
             cluster_rows = (nrow(fmat) > 2), cluster_cols = (ncol(fmat) > 2),
             display_numbers = omat, number_format = "%s",
             number_color = "black", fontsize_number = 8,
             color = fcolors, breaks = fbreaks,
             main = title, fontsize = 10, fontsize_row = 9,
             fontsize_col = 10, border_color = "grey80",
             legend_breaks = seq(0, ceiling(fmax), by = 1),
             name = "-log10(adj.p)")
    dev.off()

    png(file.path(out_dir, "tf", subdir, "tf_enrichment_heatmap.png"),
        width = fw, height = fh, units = "in", res = 150)
    pheatmap(fmat,
             cluster_rows = (nrow(fmat) > 2), cluster_cols = (ncol(fmat) > 2),
             display_numbers = omat, number_format = "%s",
             number_color = "black", fontsize_number = 8,
             color = fcolors, breaks = fbreaks,
             main = title, fontsize = 10, fontsize_row = 9,
             fontsize_col = 10, border_color = "grey80",
             legend_breaks = seq(0, ceiling(fmax), by = 1),
             name = "-log10(adj.p)")
    dev.off()
    cat(sprintf("Saved %s/tf_enrichment_heatmap.pdf / .png\n", subdir))
    return(TRUE)
  }

  build_focused_heatmap(c("Shared", "Per_CellLine", "Cross"),
    "cross_temporal", "TF Enrichment: Cross-Temporal Gene Sets")
  build_focused_heatmap(c("Divergence"),
    "cross_cellline", "TF Enrichment: Cross-Cell-Line Gene Sets")
} else {
  cat("No TFs with library-level adj.p < 0.05; skipping heatmap.\n")
}

# --- 8b. TF-Gene Regulatory Network ---
cat("\n=== TF-Gene Regulatory Networks ===\n")

sig_net <- results_df[results_df$adj_pvalue < 0.05, ]
n_sig_tfs <- length(unique(sig_net$TF))

if (nrow(sig_net) > 0 && n_sig_tfs >= 1) {
  suppressPackageStartupMessages({
    library("igraph")
    library("ggraph")
    library("tidygraph")
    library("visNetwork")
  })

  # Extract cell lines from gene_set prefixes (dynamic, any number)
  all_gs <- unique(sig_net$gene_set)
  cell_lines <- setdiff(unique(sub("_.*", "", all_gs)), c("Shared", "Cross"))
  cat(sprintf("Cell lines in enrichment: %s\n", paste(cell_lines, collapse = ", ")))

  # Define networks: one per cell line + Shared/Cross
  network_prefixes <- c(setNames(paste0("^", cell_lines, "_"), cell_lines),
                        Shared = "^Shared_|^Cross_")

  # Helper: build persistence shapes lookup
  gene_sym_to_persist <- list()
  if (exists("persist_df") && nrow(persist_df) > 0) {
    sym_to_id <- setNames(annot$gene_id, annot$gene_symbol)
    for (gene_sym in names(sym_to_id)) {
      gid <- sym_to_id[gene_sym]
      if (is.na(gid)) next
      cats <- unique(persist_df$category[persist_df$gene_id == gid])
      if (length(cats) == 1) gene_sym_to_persist[[gene_sym]] <- cats[1]
      else if (length(cats) > 1) gene_sym_to_persist[[gene_sym]] <- paste(cats, collapse = "/")
    }
  }

  direction_colors <- c("Up" = "#E41A1C", "Down" = "#377EB8",
                        "Mixed" = "grey70", "Unknown" = "grey70")
  # TF fill color
  tf_fill <- "#FF7F00"
  short_lib_net <- c("ChEA_2016" = "ChEA",
                     "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X" = "ENCODE+ChEA",
                     "TRANSFAC_and_JASPAR_PWMs" = "JASPAR",
                     "ARCHS4_TFs_Coexp" = "ARCHS4")

  for (net_name in names(network_prefixes)) {
    cat(sprintf("\n--- %s Regulatory Network ---\n", net_name))

    # Filter enrichments for this network
    net_edges_raw <- sig_net[grepl(network_prefixes[net_name], sig_net$gene_set), ]
    if (nrow(net_edges_raw) == 0) {
      cat("  No significant enrichments for this network; skipping.\n")
      next
    }

    # Build edge list
    edge_list <- list()
    for (i in seq_len(nrow(net_edges_raw))) {
      tf_name <- net_edges_raw$TF[i]
      gene_str <- net_edges_raw$genes[i]
      if (is.na(gene_str) || gene_str == "") next
      gene_syms <- trimws(strsplit(gene_str, ";")[[1]])
      gene_syms <- gene_syms[gene_syms != ""]
      for (gs in gene_syms) {
        edge_list[[length(edge_list) + 1]] <- data.frame(
          from = tf_name, to = gs,
          gene_set = net_edges_raw$gene_set[i],
          combined_score = net_edges_raw$combined_score[i],
          adj_pvalue = net_edges_raw$adj_pvalue[i],
          stringsAsFactors = FALSE
        )
      }
    }
    edges_df <- do.call(rbind, edge_list)
    cat(sprintf("  Edges: %d TF→gene\n", nrow(edges_df)))

    # Node lists
    tf_nodes   <- unique(net_edges_raw$TF)
    gene_nodes <- setdiff(unique(edges_df$to), tf_nodes)

    # TF attributes
    tf_best <- net_edges_raw %>%
      group_by(TF) %>%
      summarise(best_library = first(library), best_adjp = min(adj_pvalue),
                best_score = max(combined_score, na.rm = TRUE), .groups = "drop")
    tf_df <- data.frame(name = tf_nodes, type = "TF", stringsAsFactors = FALSE)
    tf_df$library <- tf_best$best_library[match(tf_df$name, tf_best$TF)]
    tf_df$adj_pvalue <- tf_best$best_adjp[match(tf_df$name, tf_best$TF)]
    tf_df$combined_score <- tf_best$best_score[match(tf_df$name, tf_best$TF)]
    tf_df$node_size <- -log10(tf_df$adj_pvalue)
    tf_df$node_size[!is.finite(tf_df$node_size)] <- 3
    tf_df$is_DEG <- tf_expr_info$TF_is_DEG[match(tf_df$name, tf_expr_info$TF)]
    tf_df$is_DEG[is.na(tf_df$is_DEG)] <- FALSE

    # Gene attributes
    gene_df <- data.frame(name = gene_nodes, type = "Gene", stringsAsFactors = FALSE)
    gene_df$direction <- "Mixed"
    gene_df$is_DEG    <- FALSE
    gene_df$log2FC    <- NA_real_
    gene_df$persistence <- ""

    for (i in seq_len(nrow(gene_df))) {
      sym <- gene_df$name[i]
      row <- combined_expr[combined_expr$gene_symbol == sym, ]
      if (nrow(row) == 0) next
      padj_vals <- as.numeric(unlist(row[1, padj_cols]))
      lfc_vals  <- as.numeric(unlist(row[1, lfc_cols]))
      if (all(is.na(padj_vals))) next
      best_idx <- which.min(padj_vals)
      if (length(best_idx) == 0 || is.na(padj_vals[best_idx])) next
      gene_df$log2FC[i] <- lfc_vals[best_idx]
      gene_df$is_DEG[i] <- padj_vals[best_idx] < 0.05
      if (lfc_vals[best_idx] > 0) gene_df$direction[i] <- "Up"
      else if (lfc_vals[best_idx] < 0) gene_df$direction[i] <- "Down"
      gene_df$persistence[i] <- if (!is.null(gene_sym_to_persist[[sym]]))
        gene_sym_to_persist[[sym]] else ""
    }
    gene_df$persistence[gene_df$persistence == ""] <- "Unknown"
    gene_df$persistence[grepl("/", gene_df$persistence)] <- "Divergent"

    # Degree
    gene_deg <- table(edges_df$to); tf_deg <- table(edges_df$from)
    gene_df$degree <- as.integer(gene_deg[gene_df$name]); gene_df$degree[is.na(gene_df$degree)] <- 0
    tf_df$degree   <- as.integer(tf_deg[tf_df$name]);   tf_df$degree[is.na(tf_df$degree)] <- 0

    # Combine
    tf_df$direction <- NA_character_; tf_df$log2FC <- NA_real_; tf_df$persistence <- ""
    gene_df$library <- NA_character_; gene_df$node_size <- 3
    gene_df$adj_pvalue <- NA_real_; gene_df$combined_score <- NA_real_

    all_nodes <- rbind(
      tf_df[, c("name","type","library","direction","is_DEG","log2FC","persistence","node_size","adj_pvalue","combined_score","degree")],
      gene_df[, c("name","type","library","direction","is_DEG","log2FC","persistence","node_size","adj_pvalue","combined_score","degree")]
    )

    # Build graph
    g <- graph_from_data_frame(edges_df[, c("from","to","gene_set","combined_score")],
                               vertices = all_nodes, directed = TRUE)
    if (vcount(g) < 2) { cat("  Too few nodes; skipping.\n"); next }

    comm <- cluster_louvain(as_undirected(g))
    V(g)$community <- comm$membership
    n_comm <- length(unique(comm$membership))
    cat(sprintf("  Nodes: %d (TFs=%d, genes=%d) in %d communities\n",
                vcount(g), sum(V(g)$type == "TF"), sum(V(g)$type == "Gene"), n_comm))

    for (cid in sort(unique(comm$membership))) {
      members <- V(g)$name[V(g)$community == cid]
      tfs_in_c   <- intersect(members, tf_nodes)
      genes_in_c <- intersect(members, gene_nodes)
      cat(sprintf("    C%d: TFs=[%s] + %d genes\n",
                  cid, paste(tfs_in_c, collapse=", "), length(genes_in_c)))
    }

    # --- Static plot ---
    set.seed(42)
    g_tidy <- as_tbl_graph(g)

    p_net <- ggraph(g_tidy, layout = "fr") +
      geom_edge_link(aes(width = combined_score, alpha = combined_score), color = "grey60") +
      scale_edge_width(range = c(0.15, 1.2), guide = "none") +
      scale_edge_alpha(range = c(0.15, 0.6), guide = "none") +
      geom_node_point(data = function(x) x[x$type == "Gene", ],
        aes(fill = direction, shape = persistence, size = pmin(pmax(abs(log2FC), 0.5, na.rm = TRUE), 2, na.rm = TRUE),
            stroke = ifelse(is_DEG, 1.2, 0.3)), color = "grey40") +
      geom_node_point(data = function(x) x[x$type == "TF", ],
        aes(fill = "TF", size = node_size, stroke = ifelse(is_DEG, 1.5, 0.3)),
        shape = 21, color = "black") +
      geom_node_text(aes(label = name, filter = (type == "TF" | degree >= 2 | (!is.na(log2FC) & abs(log2FC) >= 1))),
        size = 2.2, repel = TRUE, max.overlaps = 100,
        box.padding = 0.15, point.padding = 0.15, segment.color = "grey70") +
      scale_fill_manual(
        values = c(direction_colors, "TF" = tf_fill),
        breaks = c("Up", "Down", "Mixed", "TF"),
        drop = FALSE,
        na.value = "grey50",
        guide = guide_legend(title = "Direction", override.aes = list(size = 4))) +
      scale_shape_manual(
        values = c("Transient" = 21, "Sustained" = 24, "Secondary_Deferred" = 22,
                   "Partially_Sustained" = 25, "Transient_Mid" = 23, "Complex" = 18,
                   "Divergent" = 8, "Unknown" = 1),
        guide = guide_legend(title = "Persistence", override.aes = list(fill = "grey50"))) +
      scale_size_continuous(range = c(1, 8), guide = "none") +
      labs(title = paste0("TF → Target Gene Network: ", net_name),
           subtitle = paste0(length(tf_nodes), " enriched TFs, ", length(gene_nodes),
                             " target genes, ", nrow(edges_df), " edges, ", n_comm, " modules"),
            caption = "Fill=Direction · Shape=Persistence (see legend) · Size=|log2FC|") +
      theme_graph(base_family = "sans") + theme(legend.position = "bottom")

    tag <- tolower(net_name)

    # Route: cell-line + shared → cross_cellline/; divergence → cross_temporal/
    net_subdir <- if (grepl("^div", tag)) "cross_cellline" else "cross_temporal"

    ggsave(file.path(out_dir, "tf", net_subdir, paste0("tf_regulatory_network_", tag, ".pdf")),
           p_net, width = 16, height = 14)
    ggsave(file.path(out_dir, "tf", net_subdir, paste0("tf_regulatory_network_", tag, ".png")),
           p_net, width = 16, height = 14, dpi = 150)
    cat(sprintf("  Saved tf_regulatory_network_%s.pdf / .png\n", tag))

    # --- Interactive plot ---
    vis_nodes <- data.frame(
      id = V(g)$name, label = V(g)$name,
      group = ifelse(V(g)$type == "TF", "TF", V(g)$direction),
      title = ifelse(V(g)$type == "TF",
        paste0(V(g)$name, "\nLibrary: ", short_lib_net[V(g)$library],
               "\np-value: ", formatC(V(g)$adj_pvalue, format="e", digits=2),
               "\nTargets: ", V(g)$degree),
        paste0(V(g)$name, "\nDirection: ", V(g)$direction,
               ifelse(!is.na(V(g)$log2FC), paste0("\nlog2FC: ", round(V(g)$log2FC, 3)), ""),
               "\nPersistence: ", V(g)$persistence,
               ifelse(V(g)$is_DEG, "\n[IS DEG]", ""), "\nDegree: ", V(g)$degree)),
      shape = ifelse(V(g)$type == "TF", "dot",
              ifelse(V(g)$persistence == "Transient", "triangle",
              ifelse(V(g)$persistence == "Secondary_Deferred", "square",
              ifelse(V(g)$persistence == "Sustained", "diamond",
              ifelse(V(g)$persistence == "Divergent", "star", "dot"))))),
      size = ifelse(V(g)$type == "TF", V(g)$node_size * 3,
                    pmin(pmax(abs(V(g)$log2FC), 1, na.rm = TRUE), 2) * 8),
      stringsAsFactors = FALSE
    )

    vis_edges <- data.frame(
      from = edges_df$from, to = edges_df$to,
      title = paste0("Score: ", round(edges_df$combined_score, 1),
                     "\nGene set: ", edges_df$gene_set),
      value = pmin(edges_df$combined_score / 10, 5),
      stringsAsFactors = FALSE
    )

    vis_net <- visNetwork(vis_nodes, vis_edges, width = "100%", height = "800px") %>%
      visGroups(groupname = "TF", color = list(background = "#FF7F00", border = "black"),
                shape = "dot") %>%
      visGroups(groupname = "Up", color = list(background = "#E41A1C", border = "grey40")) %>%
      visGroups(groupname = "Down", color = list(background = "#377EB8", border = "grey40")) %>%
      visGroups(groupname = "Mixed", color = list(background = "grey70", border = "grey40")) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                 nodesIdSelection = TRUE) %>%
      visEdges(arrows = "to", smooth = FALSE,
               color = list(color = "#C0C0C0", highlight = "#666666"),
               scaling = list(min = 0.3, max = 2)) %>%
      visPhysics(solver = "barnesHut",
                 barnesHut = list(gravitationalConstant = -1200, centralGravity = 0.1,
                                  springLength = 250, springConstant = 0.02),
                 stabilization = list(iterations = 300, fit = TRUE)) %>%
      visInteraction(navigationButtons = TRUE, dragNodes = TRUE) %>%
      visLayout(randomSeed = 42) %>%
      visLegend(width = 0.22, position = "right",
                main = "Color = direction · Shape = persistence",
                useGroups = FALSE,
                addNodes = list(
                  list(label = "TF (enriched)", shape = "dot", color = "#FF7F00", size = 18),
                  list(label = "Gene up (log2FC > 0)", shape = "dot", color = "#E41A1C", size = 14),
                  list(label = "Gene down (log2FC < 0)", shape = "dot", color = "#377EB8", size = 14),
                  list(label = "Mixed / unknown", shape = "dot", color = "grey70", size = 14),
                  list(label = "Transient", shape = "triangle", color = "grey70", size = 14),
                  list(label = "Secondary Deferred", shape = "square", color = "grey70", size = 14),
                  list(label = "Sustained", shape = "diamond", color = "grey70", size = 14),
                  list(label = "Divergent (differs per cell line)", shape = "star", color = "grey70", size = 14)
                ))

    visSave(vis_net, file.path(out_dir, "tf", net_subdir, paste0("tf_regulatory_network_", tag, ".html")),
            selfcontained = TRUE)
    cat(sprintf("  Saved tf_regulatory_network_%s.html\n", tag))
  }

} else {
  cat("Not enough significant TFs for network visualization.\n")
}

# --- 9. Summary ---
cat("\n=== TF Enrichment Summary ===\n")

# Top hits per gene set using library-level adjusted p-value
sig_library <- results_df[results_df$adj_pvalue < 0.05, ]
if (nrow(sig_library) > 0) {
  cat(sprintf("TFs with library-level adj.p < 0.05: %d\n", nrow(sig_library)))
}

# Show top 3 TFs per gene set (best by library-level adj_pvalue)
best_per_set <- results_df %>%
  group_by(gene_set) %>%
  slice_min(adj_pvalue, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(gene_set, adj_pvalue)

for (gs in unique(best_per_set$gene_set)) {
  gs_rows <- best_per_set[best_per_set$gene_set == gs, ]
  cat(sprintf("\n  %s (%d genes):\n", gs, length(gene_set_lists[[gs]])))
  for (i in seq_len(nrow(gs_rows))) {
    r <- gs_rows[i, ]
    score_str <- if (!is.na(r$combined_score))
      sprintf(", score=%.1f", r$combined_score) else ""
    cat(sprintf("    %-30s  (%-45s) n=%d, p=%.1e%s\n",
                r$TF, r$library, r$n_overlap, r$adj_pvalue, score_str))
  }
}

# --- 9b. Active TF candidates: enriched AND differentially expressed ---
cat("\n=== Active TF Candidates (enriched + differentially expressed) ===\n")

active_tfs <- results_df[results_df$TF_is_DEG & results_df$adj_pvalue < 0.05, ]
if (nrow(active_tfs) > 0) {
  active_summary <- active_tfs %>%
    group_by(TF) %>%
    summarise(
      gene_sets    = paste(unique(gene_set), collapse = ", "),
      best_pvalue  = min(adj_pvalue),
      tf_log2FC    = first(na.omit(TF_best_log2FC)),
      tf_comparison = first(na.omit(TF_best_comparison)),
      tf_direction = if (!is.na(first(na.omit(TF_best_log2FC))) && first(na.omit(TF_best_log2FC)) > 0) "Up" else "Down",
      .groups = "drop"
    ) %>%
    arrange(best_pvalue)

  for (i in seq_len(nrow(active_summary))) {
    r <- active_summary[i, ]
    cat(sprintf("  %-12s enriched in %-30s | TF %s in %s (log2FC=%.2f, p=%.1e)\n",
                r$TF, r$gene_sets, r$tf_direction, r$tf_comparison,
                r$tf_log2FC, r$best_pvalue))
  }
} else {
  cat("  No TFs are both enriched AND differentially expressed.\n")
}

cat("\n=== TF Enrichment Complete ===\n")
