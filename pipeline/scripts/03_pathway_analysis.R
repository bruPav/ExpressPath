#!/usr/bin/env Rscript
#
# Pathway analysis: GSEA (KEGG, GO BP, Hallmark, Reactome, custom) + Pathview + GSVA
# Data-driven: all gene-set collections are built from MSigDB or a GMT file;
# no contrasts, timepoints, or cell lines are hard-coded.
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
  library("RColorBrewer")
})

cat("=== Pathway Analysis: GSEA + Pathview + GSVA ===\n")
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

combined$entrez <- id_map$ENTREZID[match(combined$gene_id, id_map$ENSEMBL)]

# --- Define contrasts (read from combined_results column names) ---
combined_cols <- colnames(combined)
contrasts <- unique(gsub("_log2FC$", "",
  grep("_log2FC$", combined_cols, value = TRUE)))
cat(sprintf("  %d contrasts found\n", length(contrasts)))

# ==============================
# Part 1: Build ranked gene lists
# ==============================
cat("\n=== Part 1: Build ranked lists ===\n")

set.seed(42)
EPS <- 1e-10

ranked_lists <- list()

for (cname in contrasts) {
  log2fc_col <- paste0(cname, "_log2FC")
  pval_col   <- paste0(cname, "_pvalue")

  if (!log2fc_col %in% names(combined) || !pval_col %in% names(combined)) next

  df <- combined[, c("gene_id", "entrez", log2fc_col, pval_col)]
  df <- df[!is.na(df$entrez) & !is.na(df[[log2fc_col]]) & !is.na(df[[pval_col]]), ]
  if (nrow(df) < 10) next

  # Cap p-values away from zero to avoid Inf ranks
  pvals <- pmax(df[[pval_col]], EPS)
  df$rank_stat <- -log10(pvals) * sign(df[[log2fc_col]])
  df <- df[order(df$rank_stat, decreasing = TRUE), ]

  gene_list <- setNames(df$rank_stat, df$entrez)
  gene_list <- gene_list[!duplicated(names(gene_list))]
  ranked_lists[[cname]] <- gene_list
}

# ==============================
# Part 2: Helper to run and save GSEA for any collection
# ==============================

empty_gsea <- function() {
  data.frame(ID = character(), Description = character(), setSize = integer(),
             enrichmentScore = numeric(), NES = numeric(), pvalue = numeric(),
             p.adjust = numeric(), qvalue = numeric(), rank = numeric(),
             leading_edge = character(), core_enrichment = character(),
             contrast = character(), stringsAsFactors = FALSE)
}

run_gsea_collection <- function(gene_lists, term2gene, collection_name) {
  cat(sprintf("\nRunning GSEA: %s\n", collection_name))
  all_res <- list()

  for (cname in names(gene_lists)) {
    gl <- gene_lists[[cname]]
    res <- tryCatch({
      GSEA(geneList = gl, TERM2GENE = term2gene, pvalueCutoff = 0.05,
           minGSSize = gsea_min, maxGSSize = gsea_max, eps = EPS, seed = 42)
    }, error = function(e) { NULL })

    if (!is.null(res) && nrow(res@result) > 0) {
      rr <- res@result
      rr$contrast <- cname
      all_res[[cname]] <- rr
    }
  }

  if (length(all_res) > 0) {
    all_df <- do.call(rbind, all_res)
    sig_df <- all_df[all_df$p.adjust < 0.05, ]
    write.table(all_df, file.path(out_dir, sprintf("gsea_%s_all.tsv", collection_name)),
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(sig_df, file.path(out_dir, sprintf("gsea_%s_signif.tsv", collection_name)),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  %s: %d total, %d significant\n",
                collection_name, nrow(all_df), nrow(sig_df)))
    return(list(all = all_df, signif = sig_df))
  } else {
    empty <- empty_gsea()
    write.table(empty, file.path(out_dir, sprintf("gsea_%s_all.tsv", collection_name)),
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(empty, file.path(out_dir, sprintf("gsea_%s_signif.tsv", collection_name)),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  %s: 0 enrichments\n", collection_name))
    return(list(all = empty, signif = empty))
  }
}

plot_gsea_dot <- function(sig_df, collection_name, n_top = 10) {
  plot_pdf <- file.path(out_dir, sprintf("gsea_dotplot_%s.pdf", collection_name))
  plot_png <- file.path(out_dir, sprintf("gsea_dotplot_%s.png", collection_name))
  if (nrow(sig_df) > 2) {
    plot_data <- sig_df %>%
      dplyr::group_by(contrast) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = n_top) %>%
      dplyr::ungroup()

    # Truncate long descriptions
    plot_data$Description <- substr(plot_data$Description, 1, 80)

    dp <- ggplot(plot_data, aes(x = NES, y = reorder(Description, NES),
                                size = -log10(p.adjust), color = p.adjust)) +
      geom_point() +
      facet_wrap(~ contrast, scales = "free_y", ncol = 3) +
      scale_color_gradient(low = "red", high = "blue") +
      theme_minimal(base_size = 9) +
      labs(x = "NES", y = "",
           title = sprintf("GSEA %s: Top %d per contrast (padj < 0.05)",
                           collection_name, n_top)) +
      theme(strip.text = element_text(size = 7),
            axis.text.y = element_text(size = 6))

    ggsave(plot_pdf, dp,
           width = 16, height = max(8, nrow(plot_data) * 0.12),
           limitsize = FALSE)
    ggsave(plot_png, dp,
           width = 16, height = max(8, nrow(plot_data) * 0.12),
           dpi = 150, limitsize = FALSE)
    cat(sprintf("  Saved %s and %s\n", basename(plot_pdf), basename(plot_png)))
  } else {
    p_blank <- ggplot() + theme_void() + labs(
      title = sprintf("GSEA %s: No significant enrichments", collection_name))
    pdf(plot_pdf, width = 6, height = 4)
    print(p_blank)
    dev.off()
    ggsave(plot_png, p_blank, width = 6, height = 4, dpi = 150)
    cat(sprintf("  Saved %s and %s (blank)\n", basename(plot_pdf), basename(plot_png)))
  }
}

# --- 2a. KEGG ---
kegg_res <- list()
for (cname in names(ranked_lists)) {
  gl <- ranked_lists[[cname]]
  res <- tryCatch({
    gseKEGG(geneList = gl, organism = "hsa", pvalueCutoff = 0.05,
            minGSSize = gsea_min, maxGSSize = gsea_max, eps = EPS, seed = 42)
  }, error = function(e) { NULL })
  if (!is.null(res) && nrow(res@result) > 0) {
    rr <- res@result
    rr$contrast <- cname
    kegg_res[[cname]] <- rr
  }
}
if (length(kegg_res) > 0) {
  kegg_all <- do.call(rbind, kegg_res)
  kegg_sig <- kegg_all[kegg_all$p.adjust < 0.05, ]
} else {
  kegg_all <- empty_gsea()
  kegg_sig <- empty_gsea()
}
write.table(kegg_all, file.path(out_dir, "gsea_kegg_all.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(kegg_sig, file.path(out_dir, "gsea_kegg_signif.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
plot_gsea_dot(kegg_sig, "kegg")

# --- 2b. GO BP (simplified) ---
go_res <- list()
for (cname in names(ranked_lists)) {
  gl <- ranked_lists[[cname]]
  res <- tryCatch({
    gseGO(geneList = gl, OrgDb = org.Hs.eg.db, ont = "BP",
          pvalueCutoff = 0.05, minGSSize = gsea_min, maxGSSize = gsea_max,
          eps = EPS, seed = 42)
  }, error = function(e) { NULL })
  if (!is.null(res) && nrow(res@result) > 0) {
    rr <- res@result
    rr$contrast <- cname
    go_res[[cname]] <- rr
  }
}
if (length(go_res) > 0) {
  go_all <- do.call(rbind, go_res)
  go_sig <- go_all[go_all$p.adjust < 0.05, ]
  # Simplify redundant GO terms
  go_sig <- tryCatch({
    simplify(go_sig, cutoff = 0.7, by = "p.adjust", select_fun = min)
  }, error = function(e) { go_sig })
} else {
  go_all <- empty_gsea()
  go_sig <- empty_gsea()
}
write.table(go_all, file.path(out_dir, "gsea_go_all.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(go_sig, file.path(out_dir, "gsea_go_signif.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
plot_gsea_dot(go_sig, "go")

# --- 2c. MSigDB Hallmarks ---
msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
term2gene_h <- msig_h[, c("gs_name", "ncbi_gene")]
names(term2gene_h) <- c("term", "gene")
term2gene_h$gene <- as.character(term2gene_h$gene)
hallmark <- run_gsea_collection(ranked_lists, term2gene_h, "hallmark")
plot_gsea_dot(hallmark$signif, "hallmark")

# --- 2d. Reactome ---
msig_r <- msigdbr(species = "Homo sapiens", collection = "C2",
                  subcollection = "CP:REACTOME")
term2gene_r <- msig_r[, c("gs_name", "ncbi_gene")]
names(term2gene_r) <- c("term", "gene")
term2gene_r$gene <- as.character(term2gene_r$gene)
reactome <- run_gsea_collection(ranked_lists, term2gene_r, "reactome")
plot_gsea_dot(reactome$signif, "reactome")

# --- 2e. Custom virus/innate-immune collection ---
custom_gmt <- file.path("resources", "gene_sets", "custom_virus_innate.gmt")
if (!file.exists(custom_gmt)) {
  custom_gmt <- file.path(dirname(sys.frame(1)$ofile), "..", "resources",
                          "gene_sets", "custom_virus_innate.gmt")
}
if (!file.exists(custom_gmt)) {
  custom_gmt <- normalizePath(file.path(out_dir, "..", "..", "pipeline",
                                        "resources", "gene_sets",
                                        "custom_virus_innate.gmt"),
                              mustWork = FALSE)
}
custom <- list(all = empty_gsea(), signif = empty_gsea())
if (file.exists(custom_gmt)) {
  term2gene_c <- read.gmt(custom_gmt)
  if (nrow(term2gene_c) > 0) {
    term2gene_c$gene <- as.character(term2gene_c$gene)
    custom <- run_gsea_collection(ranked_lists, term2gene_c, "custom")
    plot_gsea_dot(custom$signif, "custom")
  }
} else {
  cat("  custom GMT not found; skipping custom GSEA\n")
  write.table(empty_gsea(), file.path(out_dir, "gsea_custom_all.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(empty_gsea(), file.path(out_dir, "gsea_custom_signif.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# ==============================
# Part 3: Pathview Maps
# ==============================
cat("\n=== Part 2: Pathview ===\n")

pathview_contrasts <- intersect(c(
  grep("_vs_", contrasts, value = TRUE),
  grep("interaction_", contrasts, value = TRUE)
), contrasts)
if (length(pathview_contrasts) > 7) pathview_contrasts <- head(pathview_contrasts, 7)

if (nrow(kegg_sig) > 0) {
  top_kegg <- kegg_sig %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(min_padj = min(p.adjust), desc = dplyr::first(Description)) %>%
    dplyr::arrange(min_padj) %>%
    dplyr::slice_head(n = 20)

  top_kegg_ids <- top_kegg$ID
  pv_root <- normalizePath(file.path(out_dir, "pathview_output"), mustWork = FALSE)
  dir.create(pv_root, showWarnings = FALSE, recursive = TRUE)
  main_wd <- getwd()

  for (cid in top_kegg_ids) {
    cid_desc <- top_kegg$desc[top_kegg$ID == cid]
    cid_desc <- gsub("[^A-Za-z0-9_-]", "_", cid_desc)
    cat(sprintf("  Pathview: %s (%s)\n", cid, cid_desc))

    cid_dir <- file.path(pv_root, paste0(cid, "_", cid_desc))
    dir.create(cid_dir, showWarnings = FALSE, recursive = TRUE)

    for (pc in pathview_contrasts) {
      log2fc_col <- paste0(pc, "_log2FC")
      if (!log2fc_col %in% names(combined)) next

      fc_df <- combined[!is.na(combined$entrez) & !is.na(combined[[log2fc_col]]),
                        c("entrez", log2fc_col)]
      gene_fc <- setNames(fc_df[[log2fc_col]], as.character(fc_df$entrez))
      gene_fc <- gene_fc[!duplicated(names(gene_fc))]

      padj_col <- paste0(pc, "_padj")
      if (padj_col %in% names(combined)) {
        padj_vals <- combined[[padj_col]][match(names(gene_fc), as.character(combined$entrez))]
        gene_fc[is.na(padj_vals) | padj_vals >= 0.05] <- 0
      }

      if (sum(abs(gene_fc) > 0) < 1) next

      setwd(cid_dir)
      out_suffix <- gsub("_", ".", pc)

      tryCatch({
        pathview(gene.data = gene_fc, pathway.id = cid, species = "hsa",
                 gene.idtype = "ENTREZ", kegg.dir = pathview_dir,
                 out.suffix = out_suffix, limit = list(gene = 1, cpd = 1),
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

  stray_pngs <- list.files(main_wd, pattern = "^hsa.*\\.png$", full.names = TRUE)
  if (length(stray_pngs) > 0) file.remove(stray_pngs)
  stray_xmls <- list.files(main_wd, pattern = "^hsa.*\\.xml$", full.names = TRUE)
  if (length(stray_xmls) > 0) file.remove(stray_xmls)

  cat(sprintf("  Pathview maps saved to %s\n", pv_root))
} else {
  cat("  No significant KEGG enrichments for pathview\n")
}

# ==============================
# Part 4: GSVA Pathway Activity (visualization only)
# ==============================
cat("\n=== Part 3: GSVA (visualization only) ===\n")

metadata <- read.delim(file.path(results_dir, "tables", "metadata.tsv"),
                       stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample_id
metadata$cell_line  <- factor(metadata$cell_line)
metadata$time       <- factor(metadata$time, levels = unique(metadata$time))
metadata$treatment  <- factor(metadata$treatment)
metadata$group      <- factor(paste(metadata$cell_line, metadata$time,
                                    metadata$treatment, sep = "_"))

vst_counts <- as.matrix(read.table(file.path(results_dir, "tables",
                                              "vst_normalized_counts.tsv"),
                                   header = TRUE, row.names = 1, sep = "\t",
                                   check.names = FALSE))

# Hallmark gene sets
msigdb_h <- msigdbr(species = "Homo sapiens", collection = "H")
msigdb_h$entrez <- as.character(msigdb_h$ncbi_gene)
h_gene_sets <- split(msigdb_h$entrez, msigdb_h$gs_name)

ensg_to_entrez <- setNames(combined$entrez, combined$gene_id)
ensg_to_entrez <- ensg_to_entrez[!is.na(ensg_to_entrez)]
common_genes <- intersect(rownames(vst_counts), names(ensg_to_entrez))
vst_entrez <- vst_counts[common_genes, , drop = FALSE]
rownames(vst_entrez) <- ensg_to_entrez[common_genes]
dup_entrez <- duplicated(rownames(vst_entrez))
vst_entrez <- vst_entrez[!dup_entrez, , drop = FALSE]

cat(sprintf("  Expression matrix: %d genes x %d samples\n",
            nrow(vst_entrez), ncol(vst_entrez)))

cat("Running GSVA (Hallmark)...\n")
gsva_param <- gsvaParam(exprData = vst_entrez, geneSets = h_gene_sets,
                        minSize = gsea_min, maxSize = gsea_max, kcdf = "Gaussian")
gsva_res <- gsva(gsva_param)

gsva_df <- as.data.frame(gsva_res)
gsva_df$pathway <- rownames(gsva_df)
gsva_df <- gsva_df[, c("pathway", colnames(gsva_res))]
write.table(gsva_df, file = file.path(out_dir, "gsva_scores.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved gsva_scores.tsv\n")

# Sample-level heatmap of most variable Hallmark pathways
var_scores <- apply(gsva_res, 1, var, na.rm = TRUE)
top_pathways <- names(sort(var_scores, decreasing = TRUE))[seq_len(min(50, length(var_scores)))]
mat <- gsva_res[top_pathways, , drop = FALSE]
rownames(mat) <- gsub("HALLMARK_", "", rownames(mat))

annotation_col <- metadata[, c("cell_line", "time", "treatment"), drop = FALSE]
ann_colors <- list(
  cell_line = setNames(RColorBrewer::brewer.pal(max(3, length(levels(metadata$cell_line))), "Set1")[seq_along(levels(metadata$cell_line))], levels(metadata$cell_line)),
  treatment = setNames(RColorBrewer::brewer.pal(max(3, length(levels(metadata$treatment))), "Set2")[seq_along(levels(metadata$treatment))], levels(metadata$treatment)),
  time = setNames(RColorBrewer::brewer.pal(max(3, length(levels(metadata$time))), "Dark2")[seq_along(levels(metadata$time))], levels(metadata$time))
)

pdf(file.path(out_dir, "gsva_sample_heatmap.pdf"), width = 14,
    height = max(8, nrow(mat) * 0.3))
pheatmap(mat, annotation_col = annotation_col, annotation_colors = ann_colors,
         scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 8, fontsize_col = 8,
         main = "Top 50 variable Hallmark pathway scores per sample")
dev.off()
png(file.path(out_dir, "gsva_sample_heatmap.png"), width = 14,
    height = max(8, nrow(mat) * 0.3), units = "in", res = 150)
pheatmap(mat, annotation_col = annotation_col, annotation_colors = ann_colors,
         scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 8, fontsize_col = 8,
         main = "Top 50 variable Hallmark pathway scores per sample")
dev.off()
cat("  Saved gsva_sample_heatmap.pdf/png\n")

# Contrast mean-difference summary (no p-values — n=2 is underpowered)
cl_levs <- levels(metadata$cell_line)
tp_levs <- levels(metadata$time)
trt_levs <- levels(metadata$treatment)
ref_cl  <- design$factors$reference_cell_line %||% cl_levs[1]
ref_trt <- design$factors$reference_treatment %||% trt_levs[1]
nonref_trts <- setdiff(trt_levs, ref_trt)

gsva_summary <- list()
for (cl in cl_levs) {
  for (tp in tp_levs) {
    for (trt in nonref_trts) {
      tname <- paste0(cl, "_", tp, "_", trt, "_vs_", ref_trt)
      s_a <- metadata$sample_id[metadata$cell_line == cl &
                                metadata$time == tp &
                                metadata$treatment == trt]
      s_b <- metadata$sample_id[metadata$cell_line == cl &
                                metadata$time == tp &
                                metadata$treatment == ref_trt]
      if (length(s_a) == 0 || length(s_b) == 0) next
      for (pw in rownames(gsva_res)) {
        gsva_summary[[length(gsva_summary) + 1]] <- data.frame(
          contrast = tname,
          pathway = pw,
          mean_case = mean(gsva_res[pw, s_a], na.rm = TRUE),
          mean_control = mean(gsva_res[pw, s_b], na.rm = TRUE),
          mean_diff = mean(gsva_res[pw, s_a], na.rm = TRUE) - mean(gsva_res[pw, s_b], na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

if (length(gsva_summary) > 0) {
  gsva_summary_df <- do.call(rbind, gsva_summary)
  write.table(gsva_summary_df, file.path(out_dir, "gsva_contrast_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  Saved gsva_contrast_summary.tsv (%d rows)\n", nrow(gsva_summary_df)))
} else {
  write.table(data.frame(contrast = character(), pathway = character(),
                         mean_case = numeric(), mean_control = numeric(),
                         mean_diff = numeric(), stringsAsFactors = FALSE),
              file.path(out_dir, "gsva_contrast_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved gsva_contrast_summary.tsv (empty)\n")
}

# Remove old misleading differential GSVA file if it exists
old_gsva_diff <- file.path(out_dir, "gsva_diff_results.tsv")
if (file.exists(old_gsva_diff)) file.remove(old_gsva_diff)
old_gsva_dot <- file.path(out_dir, "gsva_diff_dotplot.pdf")
if (file.exists(old_gsva_dot)) file.remove(old_gsva_dot)

cat("\n=== Pathway Analysis Complete ===\n")
cat(sprintf("All outputs in: %s/\n", out_dir))
