#!/usr/bin/env Rscript
#
# TF Target Enrichment via enrichR
# Reads cross-cell-line DEG outputs and persistence classes,
# runs TF target enrichment against multiple Enrichr libraries.
#
# Usage: Rscript 03b_tf_enrichment.R <results_dir>
#

suppressPackageStartupMessages({
  library("enrichR")
  library("ggplot2")
  library("dplyr")
})

# --- Argument handling ---
if (exists("snakemake")) {
  out_dir <- dirname(snakemake@output[["cross_shared"]])
} else {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args) >= 1) args[1] else "results"
}

cat("=== TF Target Enrichment (enrichR) ===\n")
cat(sprintf("Output dir: %s\n", out_dir))

# --- Config ---
enrichr_dbs <- c("ChEA_2016",
                 "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
                 "TRANSFAC_and_JASPAR_PWMs",
                 "ARCHS4_TFs_Coexp")
min_genes <- 5       # minimum valid symbols per gene set
top_n_tfs <- 5       # top TFs per gene set for dotplot

# --- 1. Load data ---
cross_shared <- read.table(file.path(out_dir, "cross_cellline_shared.tsv"),
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
persist_df <- read.table(file.path(out_dir, "persistence_classes.tsv"),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
annot <- read.table(file.path(out_dir, "gene_annotations.tsv"),
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
  write.table(empty, file = file.path(out_dir, "tf_enrichment_results.tsv"),
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
  write.table(empty, file = file.path(out_dir, "tf_enrichment_results.tsv"),
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
write.table(results_df, file = file.path(out_dir, "tf_enrichment_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote tf_enrichment_results.tsv (%d rows)\n", nrow(results_df)))

# --- 8. Dotplot: top N TFs per gene set ---
# For each gene set, pick top TFs by global adjusted p-value
# If a TF appears in multiple libraries, keep the most significant one
top_per_set <- results_df %>%
  group_by(gene_set, TF) %>%
  summarise(
    adj_pvalue_global = min(adj_pvalue_global),
    n_overlap    = max(n_overlap),
    pvalue       = min(pvalue),
    .groups = "drop"
  ) %>%
  group_by(gene_set) %>%
  slice_min(adj_pvalue_global, n = top_n_tfs) %>%
  ungroup()

# Compute overlap ratio (overlap / gene_set_size)
gs_sizes <- sapply(gene_set_lists, length)
top_per_set$set_size <- gs_sizes[top_per_set$gene_set]
top_per_set$overlap_ratio <- top_per_set$n_overlap / top_per_set$set_size

# Clean TF names
top_per_set$TF <- gsub(" [Hh]omo [Ss]apiens", "", top_per_set$TF)
top_per_set$TF <- gsub("_HUMAN$", "", top_per_set$TF)

# Order gene sets by total number of hits
gs_order <- names(sort(sapply(gene_set_lists, length), decreasing = TRUE))
top_per_set$gene_set <- factor(top_per_set$gene_set, levels = intersect(gs_order, unique(top_per_set$gene_set)))

# Plot
p <- ggplot(top_per_set, aes(x = gene_set, y = reorder(TF, -log10(adj_pvalue_global)))) +
  geom_point(aes(size = overlap_ratio, color = -log10(adj_pvalue_global))) +
  scale_color_gradientn(
    colors = c("grey80", "#4DAF4A", "#FF7F00", "#E41A1C"),
    name = "-log10(adj.p)"
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = "Overlap ratio",
    labels = scales::percent_format()
  ) +
  labs(
    x = "Gene Set",
    y = "Transcription Factor",
    title = "TF Target Enrichment (enrichR)",
    subtitle = paste0("Top ", top_n_tfs, " TFs per gene set by adjusted p-value")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "tf_enrichment_dotplot.pdf"), p, width = 12, height = 10)
ggsave(file.path(out_dir, "tf_enrichment_dotplot.png"), p, width = 12, height = 10, dpi = 150)
cat(sprintf("Saved tf_enrichment_dotplot.pdf / .png\n"))

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

cat("\n=== TF Enrichment Complete ===\n")
