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

# Create output directory
dir.create(file.path(out_dir, "tf"), showWarnings = FALSE, recursive = TRUE)
min_genes <- 5       # minimum valid symbols per gene set
top_n_tfs <- 5       # top TFs per gene set for dotplot

# --- 1. Load data ---
cross_shared <- read.table(file.path(out_dir, "cross", "cross_cellline_shared.tsv"),
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
persist_df <- read.table(file.path(out_dir, "temporal", "persistence_classes.tsv"),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
annot <- read.table(file.path(out_dir, "tables", "gene_annotations.tsv"),
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    quote = "", fill = TRUE, comment.char = "")
colnames(annot)[1] <- "gene_id"
cat(sprintf("Loaded: cross_shared=%d, persist=%d, annotations=%d\n",
            nrow(cross_shared), nrow(persist_df), nrow(annot)))

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

  # Clean TF names
  sig_best$TF_clean <- gsub(" [Hh]omo [Ss]apiens", "", sig_best$TF)
  sig_best$TF_clean <- gsub("_HUMAN$", "", sig_best$TF_clean)

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
    else "Per_CellLine"
  })
  ann_col <- data.frame(Type = gs_type, row.names = colnames(heat_mat))
  type_colors <- setNames(c("#FF7F00", "#377EB8", "#4DAF4A"),
                          c("Shared", "Per_CellLine", "Cross"))

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
} else {
  cat("No TFs with library-level adj.p < 0.05; skipping heatmap.\n")
}

# --- 8b. TF-Gene Regulatory Network ---
cat("\n=== TF-Gene Regulatory Network ===\n")

sig_net <- results_df[results_df$adj_pvalue < 0.05, ]
n_sig_tfs <- length(unique(sig_net$TF))

if (nrow(sig_net) > 0 && n_sig_tfs >= 1) {
  suppressPackageStartupMessages({
    library("igraph")
    library("ggraph")
    library("tidygraph")
    library("visNetwork")
  })

  # Build edge list: TF → gene_symbol from enrichment genes column
  edge_list <- list()
  for (i in seq_len(nrow(sig_net))) {
    tf_name <- sig_net$TF[i]
    gene_str <- sig_net$genes[i]
    if (is.na(gene_str) || gene_str == "") next
    gene_syms <- strsplit(gene_str, ";")[[1]]
    gene_syms <- trimws(gene_syms)
    gene_syms <- gene_syms[gene_syms != ""]
    for (gs in gene_syms) {
      edge_list[[length(edge_list) + 1]] <- data.frame(
        from = tf_name,
        to = gs,
        gene_set = sig_net$gene_set[i],
        combined_score = sig_net$combined_score[i],
        adj_pvalue = sig_net$adj_pvalue[i],
        stringsAsFactors = FALSE
      )
    }
  }
  edges_df <- do.call(rbind, edge_list)
  cat(sprintf("Edges: %d TF→gene relationships\n", nrow(edges_df)))

  # Build node list: all unique TFs and target genes
  tf_nodes <- unique(sig_net$TF)
  gene_nodes <- setdiff(unique(edges_df$to), tf_nodes)

  # TF node attributes
  tf_df <- data.frame(
    name = tf_nodes,
    type = "TF",
    stringsAsFactors = FALSE
  )
  # Best library and significance per TF
  tf_best <- sig_net %>%
    group_by(TF) %>%
    summarise(
      best_library = first(library),
      best_adjp    = min(adj_pvalue),
      best_score   = max(combined_score, na.rm = TRUE),
      .groups = "drop"
    )
  tf_df$library <- tf_best$best_library[match(tf_df$name, tf_best$TF)]
  tf_df$adj_pvalue <- tf_best$best_adjp[match(tf_df$name, tf_best$TF)]
  tf_df$combined_score <- tf_best$best_score[match(tf_df$name, tf_best$TF)]
  tf_df$node_size <- -log10(tf_df$adj_pvalue)
  tf_df$node_size[!is.finite(tf_df$node_size)] <- 3

  # Look up TF expression
  tf_df$is_DEG <- tf_expr_info$TF_is_DEG[match(tf_df$name, tf_expr_info$TF)]
  tf_df$is_DEG[is.na(tf_df$is_DEG)] <- FALSE

  # Gene node attributes: look up expression in combined_results
  gene_df <- data.frame(
    name = gene_nodes,
    type = "Gene",
    stringsAsFactors = FALSE
  )

  # For each gene, find its best direction and persistence
  gene_df$direction <- "Mixed"
  gene_df$is_DEG    <- FALSE
  gene_df$log2FC     <- NA_real_
  gene_df$persistence <- ""

  for (i in seq_len(nrow(gene_df))) {
    sym <- gene_df$name[i]
    row <- combined_expr[combined_expr$gene_symbol == sym, ]
    if (nrow(row) == 0) next

    # Best comparison by padj
    padj_vals <- as.numeric(unlist(row[1, padj_cols]))
    lfc_vals  <- as.numeric(unlist(row[1, lfc_cols]))
    if (all(is.na(padj_vals))) next

    best_idx <- which.min(padj_vals)
    if (length(best_idx) == 0 || is.na(padj_vals[best_idx])) next

    gene_df$log2FC[i]  <- lfc_vals[best_idx]
    gene_df$is_DEG[i]  <- padj_vals[best_idx] < 0.05

    if (lfc_vals[best_idx] > 0) {
      gene_df$direction[i] <- "Up"
    } else if (lfc_vals[best_idx] < 0) {
      gene_df$direction[i] <- "Down"
    }
  }

  # Persistence from persist_df
  if (exists("persist_df") && nrow(persist_df) > 0) {
    gene_ids <- annot$gene_id[match(gene_df$name, annot$gene_symbol)]
    names(gene_ids) <- gene_df$name
    for (i in seq_len(nrow(gene_df))) {
      gid <- gene_ids[i]
      if (!is.na(gid)) {
        cats <- unique(persist_df$category[persist_df$gene_id == gid])
        if (length(cats) == 1) gene_df$persistence[i] <- cats[1]
        else if (length(cats) > 1) gene_df$persistence[i] <- paste(cats, collapse = "/")
      }
    }
  }
  gene_df$persistence[gene_df$persistence == ""] <- "Unknown"

  # Node degree (hub score)
  gene_degree <- table(edges_df$to)
  gene_df$degree <- as.integer(gene_degree[gene_df$name])
  gene_df$degree[is.na(gene_df$degree)] <- 0

  tf_degree <- table(edges_df$from)
  tf_df$degree <- as.integer(tf_degree[tf_df$name])
  tf_df$degree[is.na(tf_df$degree)] <- 0

  # Combine nodes
  tf_df$direction   <- NA_character_
  tf_df$log2FC      <- NA_real_
  tf_df$persistence <- ""
  gene_df$library <- NA_character_
  gene_df$node_size <- 3
  gene_df$adj_pvalue <- NA_real_
  gene_df$combined_score <- NA_real_

  all_nodes <- rbind(
    tf_df[, c("name", "type", "library", "direction", "is_DEG", "log2FC",
              "persistence", "node_size", "adj_pvalue", "combined_score", "degree")],
    gene_df[, c("name", "type", "library", "direction", "is_DEG", "log2FC",
                "persistence", "node_size", "adj_pvalue", "combined_score", "degree")]
  )

  # Build igraph
  g <- graph_from_data_frame(edges_df[, c("from", "to", "gene_set", "combined_score")],
                             vertices = all_nodes, directed = TRUE)

  # Community detection
  comm <- cluster_louvain(as.undirected(g))
  V(g)$community <- comm$membership
  n_comm <- length(unique(comm$membership))
  cat(sprintf("Nodes: %d (TFs=%d, genes=%d) in %d communities\n",
              vcount(g), sum(V(g)$type == "TF"), sum(V(g)$type == "Gene"), n_comm))

  # Community summary
  for (cid in sort(unique(comm$membership))) {
    members <- V(g)$name[V(g)$community == cid]
    tfs_in_c   <- intersect(members, tf_nodes)
    genes_in_c <- intersect(members, gene_nodes)
    cat(sprintf("  Community %d: TFs=[%s] + %d genes\n",
                cid, paste(tfs_in_c, collapse = ", "), length(genes_in_c)))
  }

  # --- Static network plot (ggraph) ---
  g_tidy <- as_tbl_graph(g)

  direction_colors <- c("Up" = "#E41A1C", "Down" = "#377EB8",
                        "Mixed" = "grey70", "Unknown" = "grey70")
  lib_colors_net <- c("ChEA_2016" = "#E41A1C",
                      "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X" = "#377EB8",
                      "TRANSFAC_and_JASPAR_PWMs" = "#4DAF4A",
                      "ARCHS4_TFs_Coexp" = "#FF7F00")

  short_lib_net <- c("ChEA_2016" = "ChEA",
                     "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X" = "ENCODE+ChEA",
                     "TRANSFAC_and_JASPAR_PWMs" = "JASPAR",
                     "ARCHS4_TFs_Coexp" = "ARCHS4")

  V(g_tidy)$library_short <- short_lib_net[V(g_tidy)$library]

  set.seed(42)
  p_net <- ggraph(g_tidy, layout = "fr") +
    # Edges
    geom_edge_link(
      aes(width = combined_score, alpha = combined_score),
      color = "grey60"
    ) +
    scale_edge_width(range = c(0.15, 1.2), guide = "none") +
    scale_edge_alpha(range = c(0.15, 0.6), guide = "none") +
    # Gene nodes
    geom_node_point(
      data = function(x) x[x$type == "Gene", ],
      aes(fill = direction, shape = persistence, size = degree + 1,
          stroke = ifelse(is_DEG, 1.2, 0.3)),
      color = "grey40"
    ) +
    # TF nodes
    geom_node_point(
      data = function(x) x[x$type == "TF", ],
      aes(fill = library_short, size = node_size, stroke = ifelse(is_DEG, 1.5, 0.3)),
      shape = 21, color = "black"
    ) +
    # Labels
    geom_node_text(
      aes(label = name, filter = (type == "TF" | degree >= 3)),
      size = 2.5, repel = TRUE, max.overlaps = 50,
      box.padding = 0.3, point.padding = 0.2, segment.color = "grey70"
    ) +
    scale_fill_manual(
      values = c(direction_colors, lib_colors_net),
      breaks = c("Up", "Down", "ChEA", "ENCODE+ChEA", "JASPAR", "ARCHS4"),
      na.value = "grey50",
      guide = guide_legend(title = "Direction / Library", override.aes = list(size = 4))
    ) +
    scale_shape_manual(
      values = c("Transient" = 21, "Sustained" = 24, "Secondary_Deferred" = 22,
                 "Unknown" = 1),
      guide = guide_legend(title = "Persistence", override.aes = list(fill = "grey50"))
    ) +
    scale_size_continuous(range = c(1, 8), guide = "none") +
    labs(
      title = "TF → Target Gene Regulatory Network",
      subtitle = paste0(n_sig_tfs, " enriched TFs, ", length(gene_nodes),
                       " target genes, ", nrow(edges_df), " edges, ", n_comm, " modules"),
      caption = "Gene: red=up, blue=down · Shape: persistence · TF: library color · Community: see console"
    ) +
    theme_graph(base_family = "sans") +
    theme(legend.position = "bottom")

  ggsave(file.path(out_dir, "tf", "tf_regulatory_network.pdf"), p_net,
         width = 14, height = 12)
  ggsave(file.path(out_dir, "tf", "tf_regulatory_network.png"), p_net,
         width = 14, height = 12, dpi = 150)
  cat("Saved tf_regulatory_network.pdf / .png\n")

  # --- Interactive network (visNetwork) ---
  vis_nodes <- data.frame(
    id = V(g)$name,
    label = V(g)$name,
    group = ifelse(V(g)$type == "TF", "TF", V(g)$direction),
    title = ifelse(V(g)$type == "TF",
      paste0(V(g)$name, "\n",
             "Library: ", short_lib_net[V(g)$library],
             "\np-value: ", formatC(V(g)$adj_pvalue, format = "e", digits = 2),
             "\nTargets: ", V(g)$degree),
      paste0(V(g)$name, "\n",
             "Direction: ", V(g)$direction,
             ifelse(!is.na(V(g)$log2FC), paste0("\nlog2FC: ", round(V(g)$log2FC, 3)), ""),
             "\nPersistence: ", V(g)$persistence,
             ifelse(V(g)$is_DEG, "\n[IS DEG]", ""),
             "\nDegree: ", V(g)$degree)
    ),
    shape = ifelse(V(g)$type == "TF", "dot", "triangle"),
    size = ifelse(V(g)$type == "TF", V(g)$node_size * 3, (V(g)$degree + 1) * 3),
    stringsAsFactors = FALSE
  )

  vis_edges <- data.frame(
    from = edges_df$from,
    to   = edges_df$to,
    title = paste0("Score: ", round(edges_df$combined_score, 1),
                   "\nGene set: ", edges_df$gene_set),
    value = pmin(edges_df$combined_score / 10, 5),
    stringsAsFactors = FALSE
  )

  vis_net <- visNetwork(vis_nodes, vis_edges, width = "100%", height = "800px") %>%
    visGroups(groupname = "TF", color = list(background = "#FF7F00", border = "black"),
              shape = "dot") %>%
    visGroups(groupname = "Up", color = list(background = "#E41A1C", border = "grey40"),
              shape = "triangle") %>%
    visGroups(groupname = "Down", color = list(background = "#377EB8", border = "grey40"),
              shape = "triangle") %>%
    visGroups(groupname = "Mixed", color = list(background = "grey70", border = "grey40"),
              shape = "triangle") %>%
    visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
               nodesIdSelection = TRUE) %>%
    visEdges(arrows = "to", smooth = FALSE,
             scaling = list(min = 1, max = 5)) %>%
    visPhysics(solver = "barnesHut",
               barnesHut = list(gravitationalConstant = -3000, centralGravity = 0.3,
                                springLength = 150, springConstant = 0.04),
               stabilization = list(iterations = 300, fit = TRUE)) %>%
    visInteraction(navigationButtons = TRUE, dragNodes = TRUE) %>%
    visLayout(randomSeed = 42) %>%
    visLegend(width = 0.2, position = "right",
              main = "Node types",
              useGroups = FALSE,
              addNodes = list(
                list(label = "TF (enriched)", shape = "dot", color = "#FF7F00", size = 20),
                list(label = "Gene up", shape = "triangle", color = "#E41A1C", size = 15),
                list(label = "Gene down", shape = "triangle", color = "#377EB8", size = 15),
                list(label = "Mixed/unknown", shape = "triangle", color = "grey70", size = 15)
              ))

  visSave(vis_net, file.path(out_dir, "tf", "tf_regulatory_network.html"),
          selfcontained = TRUE)
  cat("Saved tf_regulatory_network.html\n")

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
