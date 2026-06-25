#!/usr/bin/env Rscript
#
# Generate interactive HTML report for RNA-seq pathway analysis
# Combines GSEA results, DESeq2 stats, and pathview maps in a self-contained page
#

suppressPackageStartupMessages({
  library("htmltools")
  library("base64enc")
  library("org.Hs.eg.db")
  library("clusterProfiler")
  library("pathview")
  library("dplyr")
  library("ggplot2")
  library("jsonlite")
  library("yaml")
})

`%||%` <- function(a, b) if (is.null(a) || is.na(a) || length(a) == 0) b else a

tryRead <- function(path) {
  if (file.exists(path)) {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame()
  }
}

cat("=== Building Interactive HTML Report ===\n")

# ─── Snakemake integration ───
if (exists("snakemake")) {
  results_dir <- dirname(dirname(snakemake@output[["report"]]))
  out_dir     <- dirname(snakemake@output[["report"]])
} else {
  args <- commandArgs(trailingOnly = TRUE)
  results_dir <- if (length(args) >= 1) args[1] else "results"
  out_dir     <- if (length(args) >= 2) args[2] else file.path(results_dir, "pathway")
}
cat(sprintf("Results dir: %s\n", results_dir))
cat(sprintf("Pathway dir: %s\n", out_dir))

pv_cache  <- normalizePath(file.path(out_dir, "pathview_maps"), mustWork = FALSE)
pv_output <- normalizePath(file.path(out_dir, "pathview_output"), mustWork = FALSE)
dir.create(pv_cache, showWarnings = FALSE, recursive = TRUE)
dir.create(pv_output, showWarnings = FALSE, recursive = TRUE)

main_wd <- getwd()

# =============================================
# 1. Load data
# =============================================
cat("Loading data...\n")
combined <- read.delim(file.path(results_dir, "tables", "combined_results.tsv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
metadata <- read.delim(file.path(results_dir, "tables", "metadata.tsv"),
                       stringsAsFactors = FALSE)
gsea_kegg <- read.delim(file.path(out_dir, "gsea_kegg_signif.tsv"),
                        stringsAsFactors = FALSE)
gsea_hallmark <- read.delim(file.path(out_dir, "gsea_hallmark_signif.tsv"),
                            stringsAsFactors = FALSE)
gsea_reactome <- read.delim(file.path(out_dir, "gsea_reactome_signif.tsv"),
                            stringsAsFactors = FALSE)
gsea_custom   <- read.delim(file.path(out_dir, "gsea_custom_signif.tsv"),
                            stringsAsFactors = FALSE)
gsea_go       <- read.delim(file.path(out_dir, "gsea_go_signif.tsv"),
                            stringsAsFactors = FALSE)

# Helper to load a dotplot PNG as base64
load_png_b64 <- function(path) {
  if (file.exists(path)) {
    base64encode(readBin(path, "raw", file.info(path)$size))
  } else {
    NULL
  }
}

gsea_collection_pngs <- list(
  hallmark = load_png_b64(file.path(out_dir, "gsea_dotplot_hallmark.png")),
  reactome = load_png_b64(file.path(out_dir, "gsea_dotplot_reactome.png")),
  custom   = load_png_b64(file.path(out_dir, "gsea_dotplot_custom.png")),
  go       = load_png_b64(file.path(out_dir, "gsea_dotplot_go.png")),
  kegg     = load_png_b64(file.path(out_dir, "gsea_dotplot_kegg.png"))
)

gsva_heatmap_png <- load_png_b64(file.path(out_dir, "gsva_sample_heatmap.png"))

# --- Build gene-set enrichment tab content ---
build_enrichment_rows <- function(df, collection_label) {
  if (nrow(df) == 0) {
    return(sprintf('<tr><td colspan="5" class="text-muted">No significant %s terms</td></tr>', collection_label))
  }
  df <- df[order(df$p.adjust), ]
  df$Description <- substr(df$Description, 1, 100)
  rows <- character()
  for (i in seq_len(min(nrow(df), 500))) {
    rows <- c(rows, sprintf(
      '<tr><td>%s</td><td>%s</td><td>%.2f</td><td>%s</td><td>%s</td></tr>',
      df$contrast[i], df$Description[i], df$NES[i],
      format(df$p.adjust[i], digits = 2, scientific = TRUE),
      ifelse(df$NES[i] > 0, "Up", "Down")))
  }
  paste(rows, collapse = "\n")
}

enrichment_sections <- list()
collections <- list(
  Hallmark = list(df = gsea_hallmark, png = gsea_collection_pngs$hallmark),
  Reactome = list(df = gsea_reactome, png = gsea_collection_pngs$reactome),
  Custom   = list(df = gsea_custom,   png = gsea_collection_pngs$custom),
  GO_BP    = list(df = gsea_go,       png = gsea_collection_pngs$go),
  KEGG     = list(df = gsea_kegg,     png = gsea_collection_pngs$kegg)
)
for (cn in names(collections)) {
  item <- collections[[cn]]
  img_html <- if (!is.null(item$png)) {
    sprintf('<img src="data:image/png;base64,%s" class="img-fluid border mt-2" style="max-width:100%%" alt="%s dotplot">',
            item$png, cn)
  } else {
    '<p class="text-muted small">Dotplot not available</p>'
  }
  enrichment_sections[[cn]] <- sprintf('
<div class="tab-pane fade %s" id="enrich-%s">
  <h6>%s GSEA significant terms</h6>
  <div style="max-height:400px;overflow-y:auto">
    <table class="table table-sm table-striped small" id="enrichTable%s">
      <thead><tr><th>Contrast</th><th>Term</th><th>NES</th><th>padj</th><th>Dir</th></tr></thead>
      <tbody>%s</tbody>
    </table>
  </div>
  %s
</div>', ifelse(cn == "Hallmark", "show active", ""), tolower(cn), cn, cn,
  build_enrichment_rows(item$df, cn), img_html)
}

enrichment_nav <- paste(sapply(names(collections), function(cn) {
  sprintf('<li class="nav-item"><button class="nav-link %s" data-bs-toggle="tab" data-bs-target="#enrich-%s" type="button">%s</button></li>',
          ifelse(cn == "Hallmark", "active", ""), tolower(cn), cn)
}), collapse = "\n")

gsva_section_html <- if (!is.null(gsva_heatmap_png)) {
  sprintf('<div class="mt-4"><h6>GSVA Hallmark pathway scores per sample</h6>
    <img src="data:image/png;base64,%s" class="img-fluid border" style="max-width:100%%" alt="GSVA sample heatmap">
    <p class="text-muted small mt-1">Differential testing with only 2 replicates per group is not statistically viable, so no p-values are shown. Use this heatmap for visual exploration only.</p>
  </div>', gsva_heatmap_png)
} else {
  ""
}

# Build enrichment tab (with GSVA section inserted before closing tag)
enrichment_tab_html <- paste0('
<div class="tab-pane fade" id="enrichment">
  <p class="text-muted small">Hallmark, Reactome, Custom (innate immune / stress), GO BP, and KEGG gene-set enrichments. KEGG disease terms are retained only as an optional output; they are often driven by shared ribosomal, mitochondrial, or histone genes.</p>
  <ul class="nav nav-tabs mb-2" role="tablist">', enrichment_nav, '</ul>
  <div class="tab-content">', paste(unlist(enrichment_sections), collapse = "\n"), '</div>
  ', gsva_section_html, '
</div>')

# Load temporal analysis results
cluster_assign <- tryRead(file.path(results_dir, "cross_temporal", "cluster_assignments.tsv"))
cluster_prof  <- tryRead(file.path(results_dir, "cross_temporal", "cluster_mean_profiles.tsv"))
velocity      <- tryRead(file.path(results_dir, "cross_temporal", "velocity_summary.tsv"))
venn_genes    <- tryRead(file.path(results_dir, "cross_temporal", "venn_genelists.tsv"))
persist       <- tryRead(file.path(results_dir, "cross_temporal", "persistence_classes.tsv"))
gene_activity <- tryRead(file.path(results_dir, "cross_temporal", "gene_activity.tsv"))

# Cross-temporal analysis (Step G)
cross_tpersist <- tryRead(file.path(results_dir, "cross_cellline", "cross_temporal_persistence.tsv"))
cross_tvel      <- tryRead(file.path(results_dir, "cross_cellline", "cross_temporal_velocity.tsv"))
cross_tga       <- tryRead(file.path(results_dir, "cross_cellline", "cross_temporal_gene_activity.tsv"))
cross_tvenn     <- tryRead(file.path(results_dir, "cross_cellline", "cross_temporal_venn_genelists.tsv"))

# Map ENSG -> ENTREZ
ensg_ids <- combined$gene_id
id_map <- bitr(ensg_ids, fromType = "ENSEMBL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
combined$entrez <- id_map$ENTREZID[match(combined$gene_id, id_map$ENSEMBL)]

cat(sprintf("  %d genes, %d with ENTREZ\n", nrow(combined),
            sum(!is.na(combined$entrez))))

# =============================================
# 2. KEGG pathway gene membership
# =============================================
cat("Getting KEGG pathway gene memberships...\n")
kegg_membership <- as.list(org.Hs.egPATH2EG)

# =============================================
# 3. Define contrasts and their labels (auto from combined_results + design)
# =============================================
# ─── Read design ───
if (exists("snakemake")) {
  design <- jsonlite::fromJSON(snakemake@params$design)
} else {
  args <- commandArgs(trailingOnly = TRUE)
  design_file <- if (length(args) >= 3) args[3] else "../data/design.yaml"
  design <- yaml::read_yaml(design_file)
}

f <- design$factors
cl_list   <- f$cell_lines
tp_list   <- f$time_points
trt_list  <- f$treatments
cl_ids    <- sapply(cl_list, `[[`, "id")
tp_ids    <- sapply(tp_list, `[[`, "id")
trt_ids   <- sapply(trt_list, `[[`, "id")
cl_short  <- setNames(sapply(cl_list, `[[`, "short"), cl_ids)
tp_short  <- setNames(sapply(tp_list, `[[`, "short"), tp_ids)
trt_short <- setNames(sapply(trt_list, `[[`, "short"), trt_ids)
ref_cl    <- f$reference_cell_line %||% cl_ids[1]
ref_trt   <- f$reference_treatment %||% trt_ids[1]

# Extract contrast names from combined_results columns
all_contrasts <- unique(gsub("_log2FC$", "",
  grep("_log2FC$", names(combined), value = TRUE)))
cat(sprintf("  %d contrasts found\n", length(all_contrasts)))

# Load contrast info table (generated by 02_deseq2_analysis.R)
cinfo_path <- file.path(results_dir, "tables", "contrast_info.tsv")
if (file.exists(cinfo_path)) {
  cinfo <- read.delim(cinfo_path, stringsAsFactors = FALSE)
  rownames(cinfo) <- cinfo$contrast_name
  cat(sprintf("  Loaded contrast_info.tsv: %d contrasts\n", nrow(cinfo)))
} else {
  cinfo <- data.frame(contrast_name = all_contrasts, type = NA, cell_line = NA,
                       time = NA, treatment = NA, ref = NA, label = all_contrasts,
                       short = all_contrasts, stringsAsFactors = FALSE)
  rownames(cinfo) <- all_contrasts
}

# Build labels/shorts from contrast_info or fall back to old regex parsing
contrast_labels <- setNames(sapply(all_contrasts, function(cn) {
  if (cn %in% cinfo$contrast_name && !is.na(cinfo[cn, "label"])) return(cinfo[cn, "label"])
  gsub("_", " ", cn)
}), all_contrasts)

contrast_short <- setNames(sapply(all_contrasts, function(cn) {
  if (cn %in% cinfo$contrast_name && !is.na(cinfo[cn, "short"])) return(cinfo[cn, "short"])
  gsub("_", "", cn)
}), all_contrasts)

# Map contrast -> pathview file suffix (underscores to dots)
contrast_pv_suffix <- setNames(gsub("_", ".", all_contrasts), all_contrasts)

# =============================================
# 4. Compute per-pathway statistics
# =============================================
cat("Computing per-pathway statistics...\n")

unique_pwys <- unique(gsea_kegg$ID)
pwy_stats <- data.frame(
  ID = unique_pwys,
  Description = gsea_kegg$Description[match(unique_pwys, gsea_kegg$ID)],
  stringsAsFactors = FALSE
)

# Which contrasts are significant for each pathway?
# Build per-pathway summary by ID (name-matched, not positional)
pwy_summary <- do.call(rbind, lapply(split(gsea_kegg, gsea_kegg$ID), function(x) {
  best <- which.min(x$p.adjust)
  data.frame(
    ID = x$ID[1],
    n_sig_contrasts = nrow(x),
    contrast_list = paste(sort(unique(x$contrast)), collapse = ", "),
    best_padj = x$p.adjust[best],
    best_NES = x$NES[best],
    best_contrast = x$contrast[best],
    stringsAsFactors = FALSE
  )
}))

# Merge into pwy_stats by matching on ID (NOT by position)
idx <- match(pwy_stats$ID, pwy_summary$ID)
pwy_stats$n_sig_contrasts <- pwy_summary$n_sig_contrasts[idx]
pwy_stats$contrast_list   <- pwy_summary$contrast_list[idx]
pwy_stats$best_padj       <- pwy_summary$best_padj[idx]
pwy_stats$best_NES        <- pwy_summary$best_NES[idx]
pwy_stats$best_contrast   <- pwy_summary$best_contrast[idx]
pwy_stats$dominant_dir    <- ifelse(pwy_stats$best_NES > 0, "Up", "Down")

# Re-create sig_contrasts for pathview loop (named list: pathway -> vector of contrasts)
sig_contrasts <- split(gsea_kegg$contrast, gsea_kegg$ID)

# Per-pathway: count DE genes and baseline breakdown
pwy_stats$n_DE_total <- NA_integer_
pwy_stats$n_baseline_only <- NA_integer_
pwy_stats$n_treatment_only <- NA_integer_
pwy_stats$n_both <- NA_integer_
pwy_stats$n_btwn_only <- NA_integer_
pwy_stats$pct_up <- NA_integer_
pwy_stats$pct_down <- NA_integer_
pwy_stats$mean_abs_log2FC <- NA_real_

# Build stats directly keyed by ID (immune to reordering)
direct_lookup <- setNames(vector("list", length(unique_pwys)), unique_pwys)

for (i in seq_len(nrow(pwy_stats))) {
  pid <- pwy_stats$ID[i]
  pid_short <- sub("^hsa", "", pid)
  entrez_members <- kegg_membership[[pid_short]]
  if (is.null(entrez_members)) { next }
  members <- combined[combined$entrez %in% as.character(entrez_members), ]

  # Total DE genes in any contrast
  de_cols <- grep("_padj$", names(members), value = TRUE)
  is_de <- apply(members[, de_cols, drop = FALSE], 1, function(x) any(x < 0.05, na.rm = TRUE))

  # Treatment-vs-control DE (within any cell_line × treatment × timepoint)
  trt_vs_ctrl_contrasts <- cinfo$contrast_name[cinfo$type == "treatment_vs_control" & !is.na(cinfo$type)]
  trt_de_cols <- intersect(paste0(trt_vs_ctrl_contrasts, "_padj"), names(members))
  time_de <- if (length(trt_de_cols) > 0)
    apply(members[, trt_de_cols, drop = FALSE], 1, function(x) any(x < 0.05, na.rm = TRUE))
  else rep(FALSE, nrow(members))

  # Split into cell-line-specific for gene table
  cl1_trt_cols <- intersect(paste0(cinfo$contrast_name[cinfo$cell_line == ref_cl & cinfo$type == "treatment_vs_control" & !is.na(cinfo$type)], "_padj"), names(members))
  time_cl1_de <- if (length(cl1_trt_cols) > 0)
    apply(members[, cl1_trt_cols, drop = FALSE], 1, function(x) any(x < 0.05, na.rm = TRUE))
  else rep(FALSE, nrow(members))

  time_cl2_de <- if (length(cl_ids) >= 2) {
    nonref_cl_temp <- setdiff(cl_ids, ref_cl)[1]
    cl2_trt_cols <- intersect(paste0(cinfo$contrast_name[cinfo$cell_line == nonref_cl_temp & cinfo$type == "treatment_vs_control" & !is.na(cinfo$type)], "_padj"), names(members))
    if (length(cl2_trt_cols) > 0)
      apply(members[, cl2_trt_cols, drop = FALSE], 1, function(x) any(x < 0.05, na.rm = TRUE))
    else rep(FALSE, nrow(members))
  } else rep(FALSE, nrow(members))

  # Baseline DE (between cell lines at mock treatment or first timepoint)
  baseline_cols <- intersect(paste0(cinfo$contrast_name[cinfo$type == "between_cell_lines" & !is.na(cinfo$type)], "_padj"), names(members))
  mock_de <- if (length(baseline_cols) > 0)
    apply(members[, baseline_cols, drop = FALSE], 1, function(x) any(x < 0.05, na.rm = TRUE))
  else rep(FALSE, nrow(members))
  mock_de[is.na(mock_de)] <- FALSE

  # Between/Interaction DE
  btwn_de_cols <- intersect(paste0(cinfo$contrast_name[cinfo$type %in% c("between_cell_lines","interaction","between_treatments") & !is.na(cinfo$type)], "_padj"), names(members))
  btwn_de <- if (length(btwn_de_cols) > 0)
    apply(members[, btwn_de_cols, drop = FALSE], 1, function(x) any(x < 0.05, na.rm = TRUE))
  else rep(FALSE, nrow(members))

  # 4 mutually exclusive categories
  n_bsl_only  <- sum( mock_de & !time_de & !btwn_de, na.rm = TRUE)
  n_trt_only  <- sum(!mock_de &  time_de & !btwn_de, na.rm = TRUE)
  n_btwn_only <- sum(!mock_de & !time_de &  btwn_de, na.rm = TRUE)
  n_both      <- sum(mock_de + time_de + btwn_de >= 2, na.rm = TRUE)

  # n_DE_total = sum of all 4 (guaranteed to match)
  n_de <- n_bsl_only + n_trt_only + n_btwn_only + n_both
  pwy_stats$n_DE_total[i] <- n_de

  # Direction consistency among DE genes
  trt_fc_cols <- intersect(paste0(trt_vs_ctrl_contrasts, "_log2FC"), names(members))
  time_fcs <- unlist(sapply(trt_fc_cols, function(col) {
    if (col %in% names(members)) members[[col]][is_de] else numeric(0)
  }, simplify = FALSE))
  all_fcs <- as.numeric(time_fcs)
  all_fcs <- all_fcs[!is.na(all_fcs)]
  pct_up <- pct_dn <- 0
  if (length(all_fcs) > 0) {
    pct_up <- round(100 * sum(all_fcs > 0, na.rm = TRUE) / length(all_fcs), 0)
    pct_dn <- round(100 * sum(all_fcs < 0, na.rm = TRUE) / length(all_fcs), 0)
    pwy_stats$pct_up[i] <- pct_up
    pwy_stats$pct_down[i] <- pct_dn
  }

  # Mean abs log2FC from best contrast
  fc_col <- paste0(pwy_stats$best_contrast[i], "_log2FC")
  mean_lfc <- 0
  if (fc_col %in% names(members)) {
    mean_lfc <- mean(abs(members[[fc_col]]), na.rm = TRUE)
    pwy_stats$mean_abs_log2FC[i] <- mean_lfc
  }

  # Build DE gene table: which pathway genes are individually DE in each group?
  de_idx <- which(is_de)
  de_table_rows <- character()
  for (j in de_idx) {
    sym <- members$gene_symbol[j]
    if (is.na(sym) || sym == "--" || sym == "") sym <- members$gene_id[j]
    m <- if (mock_de[j]) "✓" else ""
    a <- if (time_cl1_de[j]) "✓" else ""
    e <- if (time_cl2_de[j]) "✓" else ""
    b <- if (btwn_de[j]) "✓" else ""
    n <- mock_de[j] + time_cl1_de[j] + time_cl2_de[j] + btwn_de[j]
    de_table_rows <- c(de_table_rows,
      sprintf('<tr><td><b>%s</b></td><td class="text-center">%s</td><td class="text-center">%s</td><td class="text-center">%s</td><td class="text-center">%s</td><td class="text-end">%d</td></tr>',
              sym, m, a, e, b, n))
  }
  de_gene_html <- paste(de_table_rows, collapse = "\n")

  # Store in direct lookup immediately (immune to reordering)
  direct_lookup[[pid]] <- list(
    n_DE_total       = n_de,
    n_baseline_only  = n_bsl_only,
    n_treatment_only = n_trt_only,
    n_btwn_only      = n_btwn_only,
    n_both           = n_both,
    pct_up           = pct_up,
    pct_down         = pct_dn,
    mean_abs_log2FC  = mean_lfc,
    de_gene_table    = de_gene_html
  )
}

# Rank: by best_padj (most significant first)
pwy_stats <- pwy_stats[order(pwy_stats$best_padj), ]
pwy_stats$rank <- seq_len(nrow(pwy_stats))

# Compute composite stat_score
pwy_stats$stat_score <- round(
  -log10(pwy_stats$best_padj) + (pwy_stats$n_sig_contrasts * 0.5) + (pwy_stats$n_DE_total / 10), 2)

cat(sprintf("  %d unique pathways analyzed\n", nrow(pwy_stats)))

# Build named lookup (immune to row reordering) merging reordered metadata + direct stats
pwy_lookup <- setNames(lapply(seq_len(nrow(pwy_stats)), function(i) {
  pid <- pwy_stats$ID[i]
  dl <- direct_lookup[[pid]]
  list(
    rank             = pwy_stats$rank[i],
    description      = pwy_stats$Description[i],
    best_padj        = pwy_stats$best_padj[i],
    best_NES         = pwy_stats$best_NES[i],
    best_contrast    = pwy_stats$best_contrast[i],
    dominant_dir     = pwy_stats$dominant_dir[i],
    n_sig_contrasts  = pwy_stats$n_sig_contrasts[i],
    contrast_list    = pwy_stats$contrast_list[i],
    n_DE_total       = if (!is.null(dl)) dl$n_DE_total else NA_integer_,
    n_baseline_only  = if (!is.null(dl)) dl$n_baseline_only else NA_integer_,
    n_treatment_only = if (!is.null(dl)) dl$n_treatment_only else NA_integer_,
    n_btwn_only      = if (!is.null(dl)) dl$n_btwn_only else NA_integer_,
    n_both           = if (!is.null(dl)) dl$n_both else NA_integer_,
    pct_up           = if (!is.null(dl)) dl$pct_up else 0L,
    pct_down         = if (!is.null(dl)) dl$pct_down else 0L,
    mean_abs_log2FC  = if (!is.null(dl)) dl$mean_abs_log2FC else 0,
    de_gene_table    = if (!is.null(dl)) dl$de_gene_table else ""
  )
}), pwy_stats$ID)

# Build core enrichment lookup: pathway/contrast -> list of ENTREZ IDs
core_lookup <- list()
for (j in seq_len(nrow(gsea_kegg))) {
  pid <- gsea_kegg$ID[j]
  ct  <- gsea_kegg$contrast[j]
  core_str <- gsea_kegg$core_enrichment[j]
  if (is.na(core_str) || core_str == "") next
  core_entrez <- strsplit(core_str, "/")[[1]]
  core_lookup[[paste(pid, ct, sep = "/")]] <- core_entrez
}

# =============================================
# 5. Generate missing Pathview maps
# =============================================
# Threshold summary for pathview maps:
#   Map coloring:   padj < 0.05   -> colored (red/blue by log2FC sign)
#                   padj >= 0.05  -> grey (log2FC forced to 0)
#   Map generation: at least 1 gene with padj < 0.05 and non-zero log2FC
#   Legend:         limit=list(gene=1) saturates |log2FC| >= 1 to max color
cat("Checking/generating pathview maps...\n")

# For each pathway, generate maps for contrasts where GSEA was significant
new_maps <- 0
for (i in seq_len(nrow(pwy_stats))) {
  pid  <- pwy_stats$ID[i]
  desc <- pwy_stats$Description[i]
  desc_safe <- gsub("[^A-Za-z0-9_-]", "_", desc)
  pv_dir <- file.path(pv_output, paste0(pid, "_", desc_safe))
  dir.create(pv_dir, showWarnings = FALSE, recursive = TRUE)

  # Contrasts where this pathway is significant
  sig_ct <- sig_contrasts[[pid]]

  for (ct in intersect(sig_ct, all_contrasts)) {
    suffix <- contrast_pv_suffix[ct]
    target_file <- file.path(pv_dir, paste0(pid, ".", suffix, ".png"))

    if (file.exists(target_file)) next  # already have it

    log2fc_col <- paste0(ct, "_log2FC")
    if (!log2fc_col %in% names(combined)) next

    fc_df <- combined[!is.na(combined$entrez) & !is.na(combined[[log2fc_col]]),
                      c("entrez", log2fc_col)]
    gene_fc <- setNames(fc_df[[log2fc_col]], as.character(fc_df$entrez))
    gene_fc <- gene_fc[!duplicated(names(gene_fc))]

    padj_col <- paste0(ct, "_padj")
    if (padj_col %in% names(combined)) {
      padj_vals <- combined[[padj_col]][match(names(gene_fc), as.character(combined$entrez))]
      gene_fc[is.na(padj_vals) | padj_vals >= 0.05] <- 0
    }

    if (sum(abs(gene_fc) > 0) < 1) next

    setwd(pv_dir)
    tryCatch({
      suppressMessages(
        pathview(gene.data = gene_fc, pathway.id = pid, species = "hsa",
                 gene.idtype = "ENTREZ", kegg.dir = pv_cache,
                  out.suffix = suffix, limit = list(gene = 1, cpd = 1),
                 low = list(gene = "blue", cpd = "blue"),
                 mid = list(gene = "grey", cpd = "grey"),
                 high = list(gene = "red", cpd = "yellow"),
                 kegg.native = TRUE, same.layer = FALSE)
      )
      new_maps <- new_maps + 1
    }, error = function(e) {
      cat(sprintf("    map failed: %s/%s: %s\n", pid, ct, e$message))
    })
    setwd(main_wd)
  }
}
cat(sprintf("  %d new maps generated\n", new_maps))

# =============================================
# 6. Base64-encode all available pathview maps
# =============================================
cat("Base64-encoding pathview maps...\n")

img_cache <- list()  # "pathway_id/contrast" -> base64 string

for (i in seq_len(nrow(pwy_stats))) {
  pid  <- pwy_stats$ID[i]
  desc <- pwy_stats$Description[i]
  desc_safe <- gsub("[^A-Za-z0-9_-]", "_", desc)
  pv_dir <- file.path(pv_output, paste0(pid, "_", desc_safe))

  if (!dir.exists(pv_dir)) next

  for (ct in all_contrasts) {
    suffix <- contrast_pv_suffix[ct]
    map_file <- file.path(pv_dir, paste0(pid, ".", suffix, ".png"))
    if (file.exists(map_file)) {
      b64 <- base64encode(readBin(map_file, "raw", file.info(map_file)$size))
      img_cache[[paste(pid, ct, sep = "/")]] <- b64
    }
  }
}
cat(sprintf("  %d images encoded\n", length(img_cache)))

# =============================================
# 7. Build GSEA dotplot as Plotly-ready data
# =============================================
cat("Preparing plot data...\n")

# Per contrast summary for bar chart
ct_summary <- as.data.frame(table(gsea_kegg$contrast))
names(ct_summary) <- c("contrast", "count")
ct_summary$label <- contrast_labels[ct_summary$contrast]

# =============================================
# 8. Build HTML document
# =============================================
cat("Building HTML...\n")

# --- Build ENTREZ -> gene info lookup for gene tables ---
entrez_to_info <- setNames(lapply(seq_len(nrow(combined)), function(i) list(
  gene_id   = combined$gene_id[i],
  symbol    = combined$gene_symbol[i]
)), combined$entrez)

# --- Helper: build leading-edge gene table for a contrast ---
# Gene table coloring:  log2FC > 0.3 -> red, < -0.3 -> blue, else grey
# Row bold:             padj < 0.05 in this contrast
# Gene source:          GSEA core_enrichment (max 20 genes)
make_gene_table <- function(pid, ct) {
  core_key <- paste(pid, ct, sep = "/")
  core_entrez <- core_lookup[[core_key]]
  if (is.null(core_entrez) || length(core_entrez) == 0) return("")

  log2fc_col <- paste0(ct, "_log2FC")
  padj_col   <- paste0(ct, "_padj")
  pval_col   <- paste0(ct, "_pvalue")
  if (!log2fc_col %in% names(combined)) return("")

  gene_rows <- c()
  for (ent in core_entrez) {
    info <- entrez_to_info[[ent]]
    # Fallback: match by gene symbol if ENTREZ lookup failed
    if (is.null(info)) {
      sym_rows <- which(combined$gene_symbol == ent)
      if (length(sym_rows) == 0) next
      r <- sym_rows[1]
    } else {
      rows <- which(combined$entrez == ent)
      if (length(rows) == 0) next
      r <- rows[1]
    }

    symbol  <- combined$gene_symbol[r]
    if (is.na(symbol) || symbol == "--" || symbol == "") symbol <- combined$gene_id[r]
    lfc     <- combined[[log2fc_col]][r]
    padj    <- if (padj_col %in% names(combined)) combined[[padj_col]][r] else NA
    pval    <- if (pval_col %in% names(combined)) combined[[pval_col]][r] else NA
    is_de   <- !is.na(padj) && padj < 0.05

    bold_span <- if (is_de) ' class="fw-bold"' else ""
    lfc_color <- if (!is.na(lfc)) {
      if (lfc > 0.3) "red" else if (lfc < -0.3) "blue" else "grey"
    } else "grey"
    padj_str <- if (is.na(padj)) "NA" else format(padj, digits = 2, scientific = TRUE)
    pval_str <- if (is.na(pval)) "NA" else format(pval, digits = 2, scientific = TRUE)
    lfc_str  <- if (is.na(lfc)) "NA" else sprintf("%+.2f", lfc)

    gene_rows <- c(gene_rows, sprintf(
      '<tr%s><td>%s</td><td style="color:%s">%s</td><td>%s</td><td>%s</td></tr>',
      bold_span, symbol, lfc_color, lfc_str, padj_str, pval_str))
    if (length(gene_rows) >= 20) break
  }
  if (length(gene_rows) == 0) return("")

  sprintf('<div class="mt-2"><small class="text-muted fw-bold">Leading-edge genes (DESeq2 padj < 0.05 shown bold):</small>
<table class="table table-sm table-borderless small mb-0" style="font-size:11px">
<thead><tr><th>Gene</th><th>log2FC</th><th>padj</th><th>pvalue</th></tr></thead><tbody>%s</tbody></table></div>',
    paste(gene_rows, collapse = "\n"))
}

# --- Helper: make pathway detail card ---
make_pwy_card <- function(pid, pwy_lookup, img_cache) {
  stats <- pwy_lookup[[pid]]
  if (is.null(stats)) return("")

  desc <- stats$description
  rank <- stats$rank

  # Subset GSEA rows for this pathway
  pwy_gsea <- gsea_kegg[gsea_kegg$ID == pid, ]
  pwy_gsea <- pwy_gsea[order(pwy_gsea$p.adjust), ]

  # Build per-contrast GSEA stats table
  gsea_rows <- ""
  for (j in seq_len(nrow(pwy_gsea))) {
    ct <- pwy_gsea$contrast[j]
    ct_label <- contrast_labels[ct]
    nes <- round(pwy_gsea$NES[j], 2)
    padj <- format(pwy_gsea$p.adjust[j], digits = 2, scientific = TRUE)
    gsea_rows <- paste0(gsea_rows,
      sprintf('<tr><td>%s</td><td>%.2f</td><td>%s</td></tr>', ct_label, nes, padj))
  }

  # Build image tabs with gene tables
  img_tabs <- ""
  img_panels <- ""
  first_tab <- TRUE

  for (ct in all_contrasts) {
    key <- paste(pid, ct, sep = "/")
    if (key %in% names(img_cache)) {
      active <- if (first_tab) "active" else ""
      show   <- if (first_tab) "show active" else ""
      ct_label <- contrast_short[ct]
      ct_full  <- contrast_labels[ct]
      img_data <- img_cache[[key]]
      gene_tbl <- make_gene_table(pid, ct)

      img_tabs <- paste0(img_tabs,
        sprintf('<button class="nav-link %s" id="tab-%s-%s" data-bs-toggle="tab"
                 data-bs-target="#panel-%s-%s" type="button" title="%s">%s</button>',
                 active, pid, ct, pid, ct, ct_full, ct_label))

      img_panels <- paste0(img_panels,
        sprintf('<div class="tab-pane fade %s" id="panel-%s-%s">
                 <div class="text-muted small mb-1">%s</div>
                 <img src="data:image/png;base64,%s" class="img-fluid border" alt="%s %s"
                      style="max-height:500px">
                 %s
                 </div>',
                 show, pid, ct, ct_full, img_data, pid, ct_full, gene_tbl))

      first_tab <- FALSE
    }
  }

  if (first_tab) {
    img_panels <- '<div class="alert alert-light small py-1 mb-0">No pathview map in any contrast &mdash; fewer than 1 gene with padj &lt; 0.05 and measurable fold change in this pathway.</div>'
  }

  baseline_html <- sprintf(
    '<span class="badge bg-primary me-1">Baseline only: %s</span>
     <span class="badge bg-success me-1">Treatment only: %s</span>
     <span class="badge bg-secondary me-1">Betw/Int only: %s</span>
     <span class="badge bg-warning text-dark me-1">Multiple: %s</span>',
    stats$n_baseline_only %||% 0,
    stats$n_treatment_only %||% 0,
    stats$n_btwn_only %||% 0,
    stats$n_both %||% 0
  )

  card <- tags$div(
    class = "card mb-3 pathway-card",
    id = paste0("pwy-", pid),
    `data-search` = paste(pid, desc, stats$contrast_list,
                          stats$best_contrast, sep = " "),
    `data-contrasts` = stats$contrast_list,
    `data-treatment-driven` = if (is.na(stats$n_baseline_only) && is.na(stats$n_treatment_only)) {
      "unknown"
    } else ifelse(((stats$n_treatment_only %||% 0) + (stats$n_btwn_only %||% 0)) >= (stats$n_baseline_only %||% 0), "true", "false"),
    `data-direction` = tolower(stats$dominant_dir),
    `data-best-padj` = stats$best_padj,
    tags$div(
      class = "card-header d-flex justify-content-between align-items-center",
      style = "cursor: pointer;",
      `data-bs-toggle` = "collapse",
      `data-bs-target` = paste0("#body-", pid),
      tags$strong(paste0("#", rank, " ", pid, " — ", desc)),
      tags$small(class = "text-muted ms-2",
                 sprintf("(best: %s)", contrast_labels[stats$best_contrast])),
      tags$div(
        class = "d-flex gap-2",
        tags$span(class = paste0("badge ", ifelse(stats$dominant_dir == "Up", "bg-danger", "bg-primary")),
                  sprintf("%s (%s)", ifelse(stats$dominant_dir == "Up", "↑ Up", "↓ Down"),
                          contrast_short[stats$best_contrast])),
        tags$span(class = "badge bg-secondary",
                  sprintf("NES %.2f", stats$best_NES)),
        tags$span(class = "badge bg-info",
                  sprintf("padj %s", format(stats$best_padj, digits = 2, scientific = TRUE))),
        tags$span(class = "badge bg-dark",
                  sprintf("%d contrasts", stats$n_sig_contrasts)),
        tags$span(class = "badge bg-light text-dark",
                  sprintf("%d sig genes", stats$n_DE_total %||% 0))
      )
    ),
    tags$div(
      id = paste0("body-", pid),
      class = "collapse",
      tags$div(
        class = "card-body",
        # Stats row
        tags$div(class = "row mb-3",
          tags$div(class = "col-md-6",
            tags$h6("Pathway Statistics"),
            HTML(sprintf(
              '<table class="table table-sm table-bordered">
               <tr><th>Best padj</th><td>%s</td><th>Best NES</th><td>%.2f</td></tr>
               <tr><th>Best contrast</th><td>%s</td><th># contrasts sig</th><td>%d</td></tr>
               <tr><th># DE genes</th><td>%d</td><th>Dominant dir</th><td>%s</td></tr>
               <tr><th>Mean |log2FC|</th><td>%.2f</td><th>%% Up / %% Down</th><td>%s%% / %s%%</td></tr>
               </table>',
               format(stats$best_padj, digits = 2, scientific = TRUE),
               stats$best_NES,
               contrast_labels[stats$best_contrast],
               stats$n_sig_contrasts,
               stats$n_DE_total %||% 0,
               stats$dominant_dir,
               stats$mean_abs_log2FC %||% 0,
               stats$pct_up %||% 0,
               stats$pct_down %||% 0
            )),
             tags$h6("Baseline DE Breakdown"),
             HTML(baseline_html),
             if (!is.null(stats$de_gene_table) && nchar(stats$de_gene_table) > 100) {
               HTML(sprintf(
                 '<details class="mt-2"><summary class="small fw-bold" style="cursor:pointer">%s %d DE Genes (click to expand)</summary>
                 <div style="max-height:300px;overflow-y:auto">
                 <table class="table table-sm table-borderless small mb-0" style="font-size:11px">
                 <thead><tr><th>Gene</th><th class="text-center">Base</th><th class="text-center">%s</th><th class="text-center">%s</th><th class="text-center">Betw/Int</th><th class="text-end">#</th></tr></thead>
                 <tbody>%s</tbody>
                 <tfoot class="fw-bold"><tr><td colspan="5">Total unique DE genes</td><td class="text-end">%d</td></tr></tfoot>
                 </table></div></details>',
                 "\u25BC", stats$n_DE_total %||% 0,
                 if (length(cl_ids) >= 1) cl_ids[1] else "CL1",
                 if (length(cl_ids) >= 2) cl_ids[2] else "CL2",
                 stats$de_gene_table, stats$n_DE_total %||% 0
               ))
             }
          ),
          tags$div(class = "col-md-6",
            tags$h6("GSEA per Contrast"),
            HTML(sprintf(
              '<table class="table table-sm table-striped"><thead><tr>
               <th>Contrast</th><th>NES</th><th>padj</th></tr></thead><tbody>%s</tbody></table>',
              gsea_rows
            ))
          )
        ),
        # Pathview maps with contrast tabs
        tags$h6("Pathview Maps"),
        tags$details(
          tags$summary("What do the tab labels mean? (click to expand)", style = "font-size:11px; color:#6c757d; cursor:pointer;"),
          tags$div(class = "ms-3 mt-1 small text-muted",
            HTML(paste(
              sapply(names(contrast_short), function(cn) {
                sprintf('<b>%s</b>=%s &nbsp;', contrast_short[cn], contrast_labels[cn])
              }), collapse = " "))
          )
        ),
        tags$details(
          tags$summary("How to read this pathway card", style = "font-size:11px; color:#6c757d; cursor:pointer;"),
          tags$div(class = "ms-3 mt-1 small text-muted",
            tags$ul(style = "padding-left: 1.2rem;",
              tags$li(HTML("<b>Map</b>: Only genes with padj &lt; 0.05 are colored (red = upregulated, blue = downregulated). Grey boxes = not significant.")),
              tags$li(HTML("<b>Leading-edge table</b> (below map): The subset of genes that drove the GSEA enrichment signal. Row bold = individually DE in this contrast.")),
              tags$li(HTML("<b>DE Genes</b> (expandable above left): All pathway member genes individually DE in at least one contrast category.")),
              tags$li(HTML("These three gene lists answer different questions and will <em>not</em> be identical. A pathway can be GSEA-significant even when few individual genes pass padj &lt; 0.05."))
            )
          )
        ),
        tags$ul(class = "nav nav-tabs mb-2", role = "tablist",
          HTML(img_tabs)
        ),
        tags$div(class = "tab-content", HTML(img_panels))
      )
    )
  )
  as.character(card)
}

# --- Build full page ---
cat("  Assembling page...\n")

# Build pathway cards
all_cards <- ""
for (i in seq_len(nrow(pwy_stats))) {
  pid <- pwy_stats$ID[i]
  all_cards <- paste0(all_cards, make_pwy_card(pid, pwy_lookup, img_cache), "\n")
}

# Build Plotly data for dotplot
plotly_dots <- list()
for (i in seq_len(nrow(gsea_kegg))) {
  plotly_dots[[i]] <- list(
    x = gsea_kegg$NES[i],
    y = gsea_kegg$Description[i],
    customdata = gsea_kegg$ID[i],
    contrast = gsea_kegg$contrast[i],
    padj = gsea_kegg$p.adjust[i],
    size = -log10(gsea_kegg$p.adjust[i])
  )
}
# Group by contrast for Plotly traces
dotplot_json <- list()
for (ct in unique(gsea_kegg$contrast)) {
  ct_data <- gsea_kegg[gsea_kegg$contrast == ct, ]
  dotplot_json[[ct]] <- list(
    x = ct_data$NES,
    y = ct_data$Description,
    customdata = ct_data$ID,
    text = paste0(ct_data$ID, "<br>", ct_data$Description, "<br>padj: ",
                  format(ct_data$p.adjust, scientific = TRUE, digits = 2)),
    marker = list(size = -log10(ct_data$p.adjust) * 3),
    type = "scatter",
    mode = "markers",
    name = contrast_labels[ct]
  )
}

# Build contrast summary data for bar chart
ct_summary_json <- list(
  x = as.list(ct_summary$count),
  y = as.list(ct_summary$label),
  type = "bar",
  orientation = "h",
  marker = list(color = "steelblue"),
  name = "Enriched pathways"
)

# =============================================
# 9. Prepare temporal analysis data for report
# =============================================
has_temporal <- nrow(cluster_prof) > 0 || nrow(velocity) > 0 || nrow(persist) > 0

# --- Clustering Plotly data ---
cluster_plotly <- "[]"
cluster_tbl_rows <- ""
cluster_select_opts <- '<option value="all" selected>All cell lines</option>'
if (nrow(cluster_prof) > 0) {
  tp_cols <- setdiff(names(cluster_prof), c("cell_line", "cluster", "n_genes"))
  traces <- list()
  cl_idx <- 0
  cl_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
                 "#F781BF", "#999999", "#E7298A", "#66A61E")
  for (cl in unique(cluster_prof$cell_line)) {
    cl_data <- cluster_prof[cluster_prof$cell_line == cl, ]
    for (i in seq_len(nrow(cl_data))) {
      cl_idx <- cl_idx + 1
      traces[[length(traces) + 1]] <- list(
        x = tp_cols,
        y = as.numeric(cl_data[i, tp_cols]),
        type = "scatter",
        mode = "lines+markers",
        name = paste0(cl, " ", cl_data$cluster[i], " (n=", cl_data$n_genes[i], ")"),
        line = list(color = cl_colors[((cl_idx - 1) %% length(cl_colors)) + 1], width = 2),
        marker = list(size = 8)
      )
    }
    cluster_select_opts <- paste0(cluster_select_opts,
      sprintf('<option value="%s">%s</option>', cl, cl))
  }
  cluster_plotly <- jsonlite::toJSON(traces, auto_unbox = TRUE, force = TRUE)

  if (nrow(cluster_assign) > 0) {
    for (i in seq_len(min(nrow(cluster_assign), 1000))) {
      row <- cluster_assign[i, ]
      memb_str <- paste(sprintf("%.2f", as.numeric(row[, grep("^C[0-9]+$", names(row))])), collapse = ", ")
      cluster_tbl_rows <- paste0(cluster_tbl_rows,
        sprintf('<tr><td>%s</td><td>%s</td><td>%s</td><td>%.3f</td><td>%s</td></tr>',
                row$cell_line, row$gene_id, row$cluster, row$membership_score, memb_str))
    }
  }
}

# --- Volcano Plotly data (for interactive viewer) ---
volcano_json <- "{}"
volcano_deg_tbl <- "{}"
volcano_select_opts <- ""
set.seed(42)
if (file.exists(cinfo_path) && nrow(cinfo) > 0) {
  volcano_data <- list()
  volcano_tables <- list()
  deg_cutoff <- 2000   # top N DEGs per contrast
  bg_n <- 500          # random non-DEGs for visual context

  for (cn in rownames(cinfo)) {
    log2fc_col <- paste0(cn, "_log2FC")
    padj_col   <- paste0(cn, "_padj")
    if (!all(c(log2fc_col, padj_col) %in% names(combined))) next

    df <- data.frame(
      gene_id = combined$gene_id,
      symbol = combined$gene_symbol,
      log2FC = combined[[log2fc_col]],
      padj = combined[[padj_col]],
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$padj), ]
    df$symbol[is.na(df$symbol) | df$symbol == "--" | df$symbol == ""] <- df$gene_id[is.na(df$symbol) | df$symbol == "--" | df$symbol == ""]
    df$is_de <- df$padj < 0.05
    df$neg_log10_padj <- pmin(-log10(df$padj), 50)

    de_df <- df[df$is_de, ]
    de_df <- de_df[order(-abs(de_df$log2FC)), ]
    if (nrow(de_df) > deg_cutoff) de_df <- de_df[1:deg_cutoff, ]
    bg_df <- df[!df$is_de, ]
    if (nrow(bg_df) > bg_n) bg_df <- bg_df[sample(nrow(bg_df), bg_n), ]
    plot_df <- rbind(de_df, bg_df)

    volcano_data[[cn]] <- list(
      x = plot_df$log2FC,
      y = plot_df$neg_log10_padj,
      text = paste0(plot_df$symbol, "<br>log2FC: ", round(plot_df$log2FC, 3),
                    "<br>padj: ", format(plot_df$padj, digits = 2, scientific = TRUE)),
      marker = list(
        color = ifelse(plot_df$is_de, ifelse(plot_df$log2FC > 0, "#E41A1C", "#377EB8"), "rgba(150,150,150,0.4)"),
        size = 4
      ),
      type = "scattergl",
      mode = "markers",
      name = cn
    )

    top_de <- head(de_df, 100)
    volcano_tables[[cn]] <- paste(
      sapply(seq_len(nrow(top_de)), function(j) {
        sprintf('<tr><td><b>%s</b></td><td>%+.3f</td><td>%s</td></tr>',
                top_de$symbol[j], top_de$log2FC[j],
                format(top_de$padj[j], digits = 2, scientific = TRUE))
      }), collapse = "\n")
  }
  volcano_json <- jsonlite::toJSON(volcano_data, auto_unbox = TRUE, force = TRUE)
  volcano_deg_tbl <- jsonlite::toJSON(volcano_tables, auto_unbox = TRUE, force = TRUE)
  volcano_select_opts <- paste(sapply(rownames(cinfo), function(cn) {
    sprintf('<option value="%s">%s</option>', cn, cinfo[cn, "label"])
  }), collapse = "\n")
  cat(sprintf("  Volcano data prepared for %d contrasts\n", length(volcano_data)))
}

# --- Velocity Plotly data ---
velocity_bar_json <- "{}"
velocity_box_json <- "{}"
velocity_tbl_html <- ""
if (nrow(velocity) > 0) {
  vel_traces <- list()
  cl_list <- unique(velocity$cell_line)
  up_colors <- c("#E41A1C", "#FF7F00", "#984EA3", "#A65628", "#F781BF", "#999999")
  dn_colors <- c("#377EB8", "#4DAF4A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  for (i in seq_along(cl_list)) {
    cl <- cl_list[i]
    cl_data <- velocity[velocity$cell_line == cl, ]
    up_c <- up_colors[((i - 1) %% length(up_colors)) + 1]
    dn_c <- dn_colors[((i - 1) %% length(dn_colors)) + 1]
    vel_traces[[length(vel_traces) + 1]] <- list(
      x = cl_data$timepoint, y = cl_data$n_up, type = "bar",
      name = paste0(cl, " Up"), marker = list(color = up_c))
    vel_traces[[length(vel_traces) + 1]] <- list(
      x = cl_data$timepoint, y = cl_data$n_down, type = "bar",
      name = paste0(cl, " Down"), marker = list(color = dn_c))
  }
  velocity_bar_json <- jsonlite::toJSON(vel_traces, auto_unbox = TRUE, force = TRUE)

  for (i in seq_len(nrow(velocity))) {
    row <- velocity[i, ]
    velocity_tbl_html <- paste0(velocity_tbl_html,
      sprintf('<tr><td>%s</td><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%.3f</td></tr>',
              row$cell_line, row$timepoint, row$n_up, row$n_down, row$n_total, row$mean_abs_log2FC))
  }
}

# --- Persistence data ---
gene_act_tbl_html <- ""

# --- Gene activity table (from gene_activity.tsv) ---
gene_act_select <- ""
if (nrow(gene_activity) > 0) {
  for (cl in unique(gene_activity$cell_line)) {
    gene_act_select <- paste0(gene_act_select,
      sprintf('<option value="%s"%s>%s</option>', cl,
              if (which(unique(gene_activity$cell_line) == cl) == 1) " selected" else "", cl))
  }
  for (i in seq_len(min(nrow(gene_activity), 2000))) {
    row <- gene_activity[i, ]
    sig_cols <- grep("^sig_", names(row), value = TRUE)
    lfc_cols <- grep("^log2FC_", names(row), value = TRUE)
    sig_html <- paste(sapply(sig_cols, function(sc) {
      sprintf('<td class="text-center">%s</td>', if (isTRUE(row[[sc]])) "&#10003;" else "")
    }), collapse = "")
    lfc_html <- paste(sapply(lfc_cols, function(lc) {
      v <- row[[lc]]
      sprintf('<td class="text-end">%s</td>', if (is.na(v)) "" else sprintf("%+.2f", v))
    }), collapse = "")
    gene_act_tbl_html <- paste0(gene_act_tbl_html,
      sprintf('<tr data-cl="%s"><td>%s</td><td>%s</td>%s%s<td>%s</td></tr>',
              row$cell_line, row$gene_symbol %||% row$gene_id, row$cell_line,
              sig_html, lfc_html, row$category))
  }
}

# --- Base64 encode Venn/UpSet and heatmap PNGs ---
# The DESeq2 script writes these per (cell_line x treatment) tag, e.g. A549_Ad26
venn_imgs <- list()
heat_imgs <- list()
nonref_trts <- setdiff(trt_ids, ref_trt)
for (cl in cl_ids) {
  for (trt in nonref_trts) {
    tag <- paste0(cl, "_", trt)
    # Try Venn first, then UpSet
    png_path <- file.path(results_dir, "cross_temporal", paste0("venn_plot_", tag, ".png"))
    if (!file.exists(png_path)) {
      png_path <- file.path(results_dir, "cross_temporal", paste0("upset_plot_", tag, ".png"))
    }
    if (file.exists(png_path)) {
      venn_imgs[[tag]] <- base64encode(readBin(png_path, "raw", file.info(png_path)$size))
    }
    # Heatmap
    heat_path <- file.path(results_dir, "cross_temporal", paste0("gene_activity_heatmap_", tag, ".png"))
    if (file.exists(heat_path)) {
      heat_imgs[[tag]] <- base64encode(readBin(heat_path, "raw", file.info(heat_path)$size))
    }
  }
}

# --- Build overlay images HTML ---
overlap_html <- ""
if (length(venn_imgs) > 0) {
  for (tag in names(venn_imgs)) {
    overlap_html <- paste0(overlap_html,
      sprintf('<div class="col-md-6 mb-3"><h6>%s</h6>
        <img src="data:image/png;base64,%s" class="img-fluid border" style="max-width:500px" alt="Venn %s"></div>',
        tag, venn_imgs[[tag]], tag))
  }
}
heat_html <- ""
if (length(heat_imgs) > 0) {
  for (tag in names(heat_imgs)) {
    heat_html <- paste0(heat_html,
      sprintf('<div class="col-md-6 mb-3"><h6>%s</h6>
        <img src="data:image/png;base64,%s" class="img-fluid border" style="max-width:500px" alt="Heatmap %s"></div>',
        tag, heat_imgs[[tag]], tag))
  }
}

# --- Gene activity table header (dynamic from data) ---
gene_act_header <- '<th>Gene</th><th>Cell Line</th>'
if (nrow(gene_activity) > 0) {
  tp_cols <- grep("^sig_", names(gene_activity), value = TRUE)
  for (tp_col in tp_cols) {
    tp_name <- sub("^sig_", "", tp_col)
    gene_act_header <- paste0(gene_act_header, sprintf('<th class="text-center">%s sig</th>', tp_name))
  }
  lfc_cols <- grep("^log2FC_", names(gene_activity), value = TRUE)
  for (lfc_col in lfc_cols) {
    tp_name <- sub("^log2FC_", "", lfc_col)
    gene_act_header <- paste0(gene_act_header, sprintf('<th class="text-end">%s log2FC</th>', tp_name))
  }
}
gene_act_header <- paste0(gene_act_header, '<th>Category</th>')

temporal_tab <- ''
if (has_temporal) {
  temporal_tab <- paste0('
<div class="tab-pane fade" id="temporal">
  <ul class="nav nav-tabs mb-3" role="tablist">
    <li class="nav-item"><button class="nav-link active" data-bs-toggle="tab" data-bs-target="#temp-clusters" type="button">Clusters</button></li>
    <li class="nav-item"><button class="nav-link" data-bs-toggle="tab" data-bs-target="#temp-velocity" type="button">Velocity</button></li>
    <li class="nav-item"><button class="nav-link" data-bs-toggle="tab" data-bs-target="#temp-persist" type="button">Persistence</button></li>
  </ul>
  <div class="tab-content">
    <!-- Clustering -->
    <div class="tab-pane fade show active" id="temp-clusters">
      <div class="row mb-3">
        <div class="col-md-3">
          <label class="form-label fw-bold">Cell Line:</label>
          <select id="clusterCL" class="form-select">', cluster_select_opts, '</select>
        </div>
      </div>
      <div id="clusterPlotly" style="height:450px"></div>
      <h6 class="mt-3">Cluster Assignments <small class="text-muted">(top 1000 genes)</small></h6>
      <div style="max-height:400px;overflow-y:auto">
        <table class="table table-sm table-striped" id="clusterTable">
          <thead><tr><th>Cell Line</th><th>Gene ID</th><th>Cluster</th><th>Score</th><th>Memberships</th></tr></thead>
          <tbody>', cluster_tbl_rows, '</tbody>
        </table>
      </div>
    </div>
    <!-- Velocity -->
    <div class="tab-pane fade" id="temp-velocity">
      <div class="row">
        <div class="col-md-6">
          <h6>DEG Count per Time Point</h6>
          <div id="velocityBar" style="height:400px"></div>
        </div>
        <div class="col-md-6">
          <h6>Summary Table</h6>
          <table class="table table-sm table-striped">
            <thead><tr><th>Cell Line</th><th>Time</th><th>Up</th><th>Down</th><th>Total</th><th>Mean |log2FC|</th></tr></thead>
            <tbody>', velocity_tbl_html, '</tbody>
          </table>
        </div>
      </div>
    </div>
    <!-- Persistence -->
    <div class="tab-pane fade" id="temp-persist">
      <h5>Gene Overlap (Venn/<i>UpSet</i>)</h5>
      <div class="row">', overlap_html, '</div>
      <h5 class="mt-3">Gene Activity Heatmaps <small class="text-muted">log2FC per timepoint, annotated by persistence category</small></h5>
      <div class="row">', heat_html, '</div>
      <h5 class="mt-3">Gene Activity Table <small class="text-muted">searchable, sortable — which gene is significant at which timepoint</small></h5>
      <div class="row mb-3">
        <div class="col-md-3">
          <label class="form-label fw-bold">Filter by Cell Line:</label>
          <select id="geneActCL" class="form-select">', gene_act_select, '</select>
        </div>
      </div>
      <div style="max-height:500px;overflow-y:auto">
        <table class="table table-sm table-striped" id="geneActTable">
          <thead><tr>', gene_act_header, '</tr></thead>
          <tbody>', gene_act_tbl_html, '</tbody>
        </table>
      </div>
    </div>
  </div>
</div>')
} else {
  temporal_tab <- '
<div class="tab-pane fade" id="temporal">
  <div class="alert alert-info">Temporal kinetics analysis was not run. Check that the temporal analysis outputs exist.</div>
</div>'
}

# ─── Cross-Temporal Divergence Tab ───
has_cross_temporal <- nrow(cross_tvel) > 0 || nrow(cross_tga) > 0

cross_divergence_tab <- ''

if (has_cross_temporal) {
  # Velocity bar chart data
  ct_velocity_bar_json <- "{}"
  ct_velocity_tbl_html <- ""
  if (nrow(cross_tvel) > 0) {
    ct_vel_traces <- list()
    pairs <- unique(cross_tvel$pair)
    colors <- c("#FF7F00", "#4DAF4A", "#377EB8", "#E41A1C", "#984EA3", "#A65628")
    for (i in seq_along(pairs)) {
      pr <- pairs[i]
      pr_data <- cross_tvel[cross_tvel$pair == pr, ]
      cl <- colors[((i - 1) %% length(colors)) + 1]
      ct_vel_traces[[length(ct_vel_traces) + 1]] <- list(
        x = pr_data$timepoint, y = pr_data$n_up, type = "bar",
        name = paste0(pr, " Up"), marker = list(color = cl))
      ct_vel_traces[[length(ct_vel_traces) + 1]] <- list(
        x = pr_data$timepoint, y = pr_data$n_down, type = "bar",
        name = paste0(pr, " Down"), marker = list(color = cl), opacity = 0.5)
    }
    ct_velocity_bar_json <- jsonlite::toJSON(ct_vel_traces, auto_unbox = TRUE, force = TRUE)
    for (i in seq_len(nrow(cross_tvel))) {
      row <- cross_tvel[i, ]
      ct_velocity_tbl_html <- paste0(ct_velocity_tbl_html,
        sprintf('<tr><td>%s</td><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%.3f</td></tr>',
                row$pair, row$timepoint, row$n_up, row$n_down, row$n_total, row$mean_abs_log2FC))
    }
  }

  # Base64 encode cross-temporal images
  ct_venn_imgs <- list()
  ct_heat_imgs <- list()
  for (pr in unique(c(cross_tvel$pair, cross_tga$pair))) {
    tag <- gsub("_vs_", "v", pr)
    # Venn
    venn_path <- file.path(results_dir, "cross_cellline", paste0("cross_temporal_venn_", tag, ".png"))
    if (file.exists(venn_path)) {
      ct_venn_imgs[[pr]] <- base64encode(readBin(venn_path, "raw", file.info(venn_path)$size))
    }
    # Heatmap
    heat_path <- file.path(results_dir, "cross_cellline", paste0("cross_temporal_activity_heatmap_", tag, ".png"))
    if (file.exists(heat_path)) {
      ct_heat_imgs[[pr]] <- base64encode(readBin(heat_path, "raw", file.info(heat_path)$size))
    }
  }

  ct_venn_html <- ""
  if (length(ct_venn_imgs) > 0) {
    for (pr in names(ct_venn_imgs)) {
      ct_venn_html <- paste0(ct_venn_html,
        sprintf('<div class="col-md-6 mb-3"><div class="card"><div class="card-header fw-bold">%s</div><div class="card-body text-center"><img src="data:image/png;base64,%s" class="img-fluid" style="max-height:500px"></div></div></div>',
                pr, ct_venn_imgs[[pr]]))
    }
  }

  ct_heat_html <- ""
  if (length(ct_heat_imgs) > 0) {
    for (pr in names(ct_heat_imgs)) {
      ct_heat_html <- paste0(ct_heat_html,
        sprintf('<div class="col-md-12 mb-3"><div class="card"><div class="card-header fw-bold">%s</div><div class="card-body text-center"><img src="data:image/png;base64,%s" class="img-fluid"></div></div></div>',
                pr, ct_heat_imgs[[pr]]))
    }
  }

  # Gene activity table
  ct_gene_act_tbl_html <- ""
  ct_gene_act_header <- '<th>Gene</th><th>Pair</th>'
  if (nrow(cross_tga) > 0) {
    sig_cols <- grep("^sig_", names(cross_tga), value = TRUE)
    for (sc in sig_cols) {
      ct_gene_act_header <- paste0(ct_gene_act_header,
        sprintf('<th class="text-center">%s sig</th>', sub("^sig_", "", sc)))
    }
    lfc_cols <- grep("^log2FC_", names(cross_tga), value = TRUE)
    for (lc in lfc_cols) {
      ct_gene_act_header <- paste0(ct_gene_act_header,
        sprintf('<th class="text-end">%s log2FC</th>', sub("^log2FC_", "", lc)))
    }
    ct_gene_act_header <- paste0(ct_gene_act_header, '<th>Divergence Category</th>')

    for (i in seq_len(min(nrow(cross_tga), 2000))) {
      row <- cross_tga[i, ]
      sig_html <- paste(sapply(sig_cols, function(sc) {
        sprintf('<td class="text-center">%s</td>', if (isTRUE(row[[sc]])) "&#10003;" else "")
      }), collapse = "")
      lfc_html <- paste(sapply(lfc_cols, function(lc) {
        v <- row[[lc]]
        sprintf('<td class="text-end">%s</td>', if (is.na(v)) "" else sprintf("%+.2f", v))
      }), collapse = "")
      ct_gene_act_tbl_html <- paste0(ct_gene_act_tbl_html,
        sprintf('<tr><td>%s</td><td>%s</td>%s%s<td>%s</td></tr>',
                row$gene_symbol %||% row$gene_id, row$pair, sig_html, lfc_html, row$category))
    }
  }

  # Full divergence tab
  cross_divergence_tab <- paste0('
<div class="tab-pane fade" id="divergence">
  <ul class="nav nav-tabs mb-3" role="tablist">
    <li class="nav-item"><button class="nav-link active" data-bs-toggle="tab" data-bs-target="#div-velocity" type="button">Divergence Velocity</button></li>
    <li class="nav-item"><button class="nav-link" data-bs-toggle="tab" data-bs-target="#div-patterns" type="button">Divergence Patterns</button></li>
    <li class="nav-item"><button class="nav-link" data-bs-toggle="tab" data-bs-target="#div-genes" type="button">Gene Table</button></li>
  </ul>
  <div class="tab-content">
    <!-- Velocity -->
    <div class="tab-pane fade show active" id="div-velocity">
      <div class="row">
        <div class="col-md-6">
          <h6>Between-Cell-Line DEG Count per Time Point</h6>
          <div id="ctVelocityBar" style="height:400px"></div>
        </div>
        <div class="col-md-6">
          <h6>Summary Table</h6>
          <table class="table table-sm table-striped">
            <thead><tr><th>Pair</th><th>Time</th><th>Up</th><th>Down</th><th>Total</th><th>Mean |log2FC|</th></tr></thead>
            <tbody>', ct_velocity_tbl_html, '</tbody>
          </table>
        </div>
      </div>
    </div>
    <!-- Patterns -->
    <div class="tab-pane fade" id="div-patterns">
      <h5>Cross-Timepoint Gene Overlap <small class="text-muted">Venn diagrams of between-cell-line DEGs across timepoints</small></h5>
      <div class="row">', ct_venn_html, '</div>
      <h5 class="mt-3">Cell Line Divergence Heatmaps <small class="text-muted">between-cell-line log2FC per timepoint, annotated by divergence category</small></h5>
      <div class="row">', ct_heat_html, '</div>
    </div>
    <!-- Gene Table -->
    <div class="tab-pane fade" id="div-genes">
      <h5>Cross-Temporal Gene Activity Table <small class="text-muted">searchable, sortable — which genes differ between cell lines at which timepoint</small></h5>
      <div style="max-height:500px;overflow-y:auto">
        <table class="table table-sm table-striped" id="ctGeneActTable">
          <thead><tr>', ct_gene_act_header, '</tr></thead>
          <tbody>', ct_gene_act_tbl_html, '</tbody>
        </table>
      </div>
    </div>
  </div>
</div>')
} else {
  cross_divergence_tab <- '
<div class="tab-pane fade" id="divergence">
  <div class="alert alert-info">Cross-cell-line temporal analysis was not run. Check that cross-temporal analysis outputs exist.</div>
</div>'
}

# --- Build head HTML manually (htmltools drops head tags) ---
head_html <- sprintf('
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>RNA-seq Pathway Analysis — Interactive Report</title>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
  <link rel="stylesheet" href="https://cdn.datatables.net/1.13.4/css/dataTables.bootstrap5.min.css">
  <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
  <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
  <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
  <script src="https://cdn.datatables.net/1.13.4/js/dataTables.bootstrap5.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
  <style>
    .pathway-card { transition: box-shadow 0.3s; }
    .pathway-card:hover { box-shadow: 0 4px 12px rgba(0,0,0,0.15); }
    .pathway-card.highlight { box-shadow: 0 0 0 3px #0d6efd; }
    .nav-tabs .nav-link { font-size: 11px; padding: 4px 8px; }
    .filter-bar { background: #f8f9fa; padding: 12px; border-radius: 6px; margin-bottom: 15px; }
    body { padding-bottom: 40px; }
  </style>
</head>')

# --- Build page as a single string ---
# tags$body works fine, get its HTML
body_html <- as.character(tags$body(
    tags$div(class = "container-fluid",

      # Header
      tags$div(class = "row mt-3 mb-4",
        tags$div(class = "col-12 text-center",
          tags$h3("RNA-seq Pathway Analysis — Interactive Report"),
          tags$p(class = "text-muted",
            sprintf("%s cell lines | treatments: %s | %d genes | %d enriched KEGG pathways",
                    paste(cl_ids, collapse=" & "),
                    paste(trt_ids, collapse=", "),
                    nrow(combined), nrow(pwy_stats)))
        )
      ),

      # Summary stats row
      tags$div(class = "row mb-3",
        tags$div(class = "col-md-2",
          tags$div(class = "card text-center bg-light", tags$div(class = "card-body",
            tags$h5(class = "text-primary", nrow(pwy_stats)), tags$small("Enriched pathways")))),
        tags$div(class = "col-md-2",
          tags$div(class = "card text-center bg-light", tags$div(class = "card-body",
            tags$h5(class = "text-success", format(nrow(gsea_kegg), big.mark = ",")),
            tags$small("Pathway × Contrast enrichments")))),
        tags$div(class = "col-md-2",
          tags$div(class = "card text-center bg-light", tags$div(class = "card-body",
            tags$h5(class = "text-danger", sum(pwy_stats$n_DE_total, na.rm = TRUE)),
            tags$small("Total DE gene-pathway hits")))),
        tags$div(class = "col-md-3",
          tags$div(class = "card text-center bg-light", tags$div(class = "card-body",
            tags$h5(class = "text-info", sum(pwy_stats$n_sig_contrasts >= 3, na.rm = TRUE)),
            tags$small("Pathways sig in 3+ contrasts (robust)")))),
        tags$div(class = "col-md-3",
          tags$div(class = "card text-center bg-light", tags$div(class = "card-body",
            tags$h5(class = "text-warning", length(img_cache)),
            tags$small("Pathview maps embedded"))))
      ),

      # Tab navigation
      tags$ul(class = "nav nav-tabs mb-3", id = "mainTabs", role = "tablist",
        tags$li(class = "nav-item",
          tags$button(class = "nav-link active", id = "tab-overview", `data-bs-toggle` = "tab",
                      `data-bs-target` = "#overview", type = "button", "Overview")),
        tags$li(class = "nav-item",
          tags$button(class = "nav-link", id = "tab-enrichment", `data-bs-toggle` = "tab",
                      `data-bs-target` = "#enrichment", type = "button", "Gene-set Enrichment")),
        tags$li(class = "nav-item",
          tags$button(class = "nav-link", id = "tab-pathways", `data-bs-toggle` = "tab",
                      `data-bs-target` = "#pathways", type = "button", "KEGG Pathview")),
        tags$li(class = "nav-item",
          tags$button(class = "nav-link", id = "tab-temporal", `data-bs-toggle` = "tab",
                      `data-bs-target` = "#temporal", type = "button", "Temporal Kinetics")),
        tags$li(class = "nav-item",
          tags$button(class = "nav-link", id = "tab-volcano", `data-bs-toggle` = "tab",
                      `data-bs-target` = "#volcano", type = "button", "Volcano Viewer")),
        if (has_cross_temporal) tags$li(class = "nav-item",
          tags$button(class = "nav-link", id = "tab-divergence", `data-bs-toggle` = "tab",
                      `data-bs-target` = "#divergence", type = "button", "Cell Line Divergence"))
        else HTML("")
      ),

      tags$div(class = "tab-content",
        # --- OVERVIEW TAB ---
        tags$div(class = "tab-pane fade show active", id = "overview",
          tags$div(class = "row",
            tags$div(class = "col-md-8",
              tags$h5("KEGG GSEA overview (click a dot to jump to KEGG pathway card)"),
              tags$div(id = "gseaDotplot", style = "height: 600px;")
            ),
            tags$div(class = "col-md-4",
              tags$h5("Enrichments per Contrast"),
              tags$div(id = "contrastBarplot", style = "height: 600px;")
            )
          )
        ),

        # --- GENE-SET ENRICHMENT TAB ---
        HTML(enrichment_tab_html),

        # --- KEGG PATHVIEW TAB ---
        tags$div(class = "tab-pane fade", id = "pathways",
          # Filter bar
          tags$div(class = "filter-bar",
            tags$div(class = "row align-items-center",
              tags$div(class = "col-md-3",
                tags$label(class = "form-label fw-bold", "Search:"),
                tags$input(type = "text", id = "pathwaySearch", class = "form-control",
                           placeholder = "Pathway name, gene, KEGG ID...")
              ),
              tags$div(class = "col-md-6",
                tags$label(class = "form-label fw-bold", "Filter:"),
                tags$div(class = "d-flex gap-3 flex-wrap",
                  tags$div(class = "form-check",
                    tags$input(class = "form-check-input filter-check", type = "checkbox", id = "filterTreatment"),
                    tags$label(class = "form-check-label", `for` = "filterTreatment", "Treatment-driven only")
                  ),
                  tags$div(class = "form-check",
                    tags$input(class = "form-check-input filter-check", type = "checkbox", id = "filterPad01"),
                    tags$label(class = "form-check-label", `for` = "filterPad01", "padj < 0.01 only")
                  ),
                  tags$div(class = "form-check",
                    tags$input(class = "form-check-input filter-check", type = "checkbox", id = "filterUp"),
                    tags$label(class = "form-check-label", `for` = "filterUp", "Up-regulated only")
                  ),
                  tags$div(class = "form-check",
                    tags$input(class = "form-check-input filter-check", type = "checkbox", id = "filterDown"),
                    tags$label(class = "form-check-label", `for` = "filterDown", "Down-regulated only")
                  )
                )
              ),
              tags$div(class = "col-md-3 text-end",
                tags$button(class = "btn btn-outline-secondary btn-sm", onclick = "expandAll()", "Expand All"),
                tags$button(class = "btn btn-outline-secondary btn-sm ms-1", onclick = "collapseAll()", "Collapse All")
              )
            )
          ),
          # Pathway cards
          HTML(all_cards)
        ),
        # Temporal Kinetics tab
        HTML(temporal_tab),
        # Volcano Viewer tab
        HTML(sprintf('<div class="tab-pane fade" id="volcano">
          <div class="row mb-2">
            <div class="col-md-4">
              <label class="form-label fw-bold">Select Contrast:</label>
              <select id="volcanoContrast" class="form-select">%s</select>
            </div>
            <div class="col-md-4">
              <small id="volcanoDEGcount" class="text-muted"></small>
            </div>
          </div>
          <div class="row">
            <div class="col-md-7">
              <div id="volcanoPlotly" style="height:550px"></div>
            </div>
            <div class="col-md-5">
              <h6>Top DEGs</h6>
              <div style="max-height:520px;overflow-y:auto">
                <table class="table table-sm table-hover small" style="font-size:11px">
                  <thead><tr><th>Gene</th><th>log2FC</th><th>padj</th></tr></thead>
                  <tbody id="volcanoDEGtable"></tbody>
                </table>
              </div>
            </div>
          </div>
        </div>', volcano_select_opts)),
        # Cell Line Divergence tab
        HTML(cross_divergence_tab)
      )
    ),

    # JavaScript
    HTML(sprintf('<script>
$(document).ready(function() {

  $("#pathwaySearch").on("keyup", updateVisibility);

  $(".filter-check").on("change", updateVisibility);

  function updateVisibility() {
    var searchVal = $("#pathwaySearch").val().toLowerCase();
    var treatmentOnly = $("#filterTreatment").is(":checked");
    var pad01 = $("#filterPad01").is(":checked");
    var onlyUp = $("#filterUp").is(":checked");
    var onlyDown = $("#filterDown").is(":checked");
    console.log("updateVisibility called, searchVal=" + searchVal);

    var totalCards = $(".pathway-card").length;
    console.log("Found " + totalCards + " pathway cards");

    var visibleCount = 0, hiddenCount = 0;
    $(".pathway-card").each(function() {
      var card = $(this);
      var show = true;

      if (searchVal) {
        var searchData = (card.attr("data-search") || "").toLowerCase();
        if (hiddenCount + visibleCount < 3) {
          console.log("card " + card.attr("id") + " data-search=" + searchData);
        }
        if (searchData.indexOf(searchVal) === -1) show = false;
      }
      if (treatmentOnly && card.attr("data-treatment-driven") !== "true") show = false;
      if (pad01) {
        var padj = parseFloat(card.attr("data-best-padj"));
        if (isNaN(padj) || padj >= 0.01) show = false;
      }
      if (onlyUp && card.attr("data-direction") !== "up") show = false;
      if (onlyDown && card.attr("data-direction") !== "down") show = false;

      card.toggle(show);
      if (show) visibleCount++; else hiddenCount++;
    });
    console.log("Visible: " + visibleCount + ", Hidden: " + hiddenCount);
  }

  window.expandAll = function() {
    $(".pathway-card .collapse").addClass("show");
  };
  window.collapseAll = function() {
    $(".pathway-card .collapse").removeClass("show");
  };
});

var gseaDotData = %s;
var gseaDotLayout = {
  title: "",
  xaxis: { title: "NES", zeroline: true },
  yaxis: { title: "", automargin: true, tickfont: { size: 9 } },
  margin: { l: 300, r: 20, t: 30, b: 50 },
  height: 600,
  showlegend: true,
  legend: { font: { size: 9 }, orientation: "h", y: 1.1 },
  hovermode: "closest"
};

Plotly.newPlot("gseaDotplot", Object.values(gseaDotData), gseaDotLayout);

var ctBarData = [%s];
var ctBarLayout = {
  title: "",
  xaxis: { title: "Number of enriched pathways" },
  yaxis: { title: "", automargin: true },
  margin: { l: 180, r: 20, t: 30, b: 30 },
  height: 600
};
Plotly.newPlot("contrastBarplot", ctBarData, ctBarLayout);

document.getElementById("gseaDotplot").on("plotly_click", function(data) {
  if (data.points && data.points.length > 0) {
    var pwyId = data.points[0].customdata;
    var pwTab = new bootstrap.Tab(document.getElementById("tab-pathways"));
    pwTab.show();
    setTimeout(function() {
      var card = document.getElementById("pwy-" + pwyId);
      if (card) {
        card.querySelector(".collapse").classList.add("show");
        card.scrollIntoView({ behavior: "smooth", block: "center" });
        $(card).addClass("highlight");
        setTimeout(function() { $(card).removeClass("highlight"); }, 2000);
      }
    }, 300);
  }
});

// Volcano viewer
var volcanoData = %s;
var volcanoTables = %s;
var volcanoRendered = false;

$("#volcanoContrast").on("change", function() {
  var cn = $(this).val();
  if (!cn || !volcanoData[cn]) return;
  var data = [volcanoData[cn]];
  var layout = {
    title: cn,
    xaxis: { title: "log2 Fold Change", zeroline: true, zerolinecolor: "#999" },
    yaxis: { title: "-log10(adj p-value)" },
    showlegend: false,
    margin: { l: 60, r: 20, t: 40, b: 50 },
    shapes: [
      { type: "line", x0: -1, x1: -1, y0: 0, y1: 1, yref: "paper",
        line: { color: "grey", dash: "dot", width: 1 } },
      { type: "line", x0: 1, x1: 1, y0: 0, y1: 1, yref: "paper",
        line: { color: "grey", dash: "dot", width: 1 } },
      { type: "line", x0: null, x1: null, y0: -Math.log10(0.05), y1: -Math.log10(0.05),
        xref: "paper", line: { color: "blue", dash: "dash", width: 1 } }
    ]
  };
  if (volcanoRendered) {
    Plotly.react("volcanoPlotly", data, layout);
  } else {
    volcanoRendered = true;
    Plotly.newPlot("volcanoPlotly", data, layout);
  }
  if (volcanoTables[cn]) {
    $("#volcanoDEGtable").html(volcanoTables[cn]);
  } else {
    $("#volcanoDEGtable").html("");
  }
});

document.getElementById("tab-volcano").addEventListener("shown.bs.tab", function() {
  if (!$("#volcanoContrast").val()) {
    var firstOpt = $("#volcanoContrast option").first().val();
    if (firstOpt) $("#volcanoContrast").val(firstOpt).trigger("change");
  } else {
    $("#volcanoContrast").trigger("change");
  }
});
</script>',
      jsonlite::toJSON(unname(dotplot_json), auto_unbox = TRUE, force = TRUE),
      jsonlite::toJSON(ct_summary_json, auto_unbox = TRUE, force = TRUE),
      volcano_json,
      volcano_deg_tbl
    )),
    # Temporal analysis JavaScript
    if (has_temporal) HTML(sprintf('<script>
$(document).ready(function() {
  var clusterRendered = false, velocityRendered = false;

  function renderCluster(filterCL) {
    if (!filterCL) filterCL = "all";
    var allClusterData = %s;
    var clusterLayout = {
      title: "", xaxis: { title: "Time Point" }, yaxis: { title: "Mean VST Expression" },
      height: 450, margin: { l: 60, r: 20, t: 30, b: 50 },
      legend: { font: { size: 10 }, orientation: "h", y: 1.15 }
    };
    var filtered = [];
    for (var i = 0; i < allClusterData.length; i++) {
      if (filterCL === "all" || allClusterData[i].name.indexOf(filterCL) === 0) {
        filtered.push(allClusterData[i]);
      }
    }
    if (clusterRendered) {
      Plotly.react("clusterPlotly", filtered, clusterLayout);
    } else {
      clusterRendered = true;
      Plotly.newPlot("clusterPlotly", filtered, clusterLayout);
    }
  }

  $("#clusterCL").on("change", function() {
    renderCluster($(this).val());
  });

  function renderVelocity() {
    if (velocityRendered) return;
    velocityRendered = true;
    var velocityData = %s;
    var velLayout = {
      barmode: "group",
      title: "", xaxis: { title: "Time Point" }, yaxis: { title: "Number of DEGs" },
      height: 400, margin: { l: 60, r: 20, t: 30, b: 50 },
      legend: { font: { size: 10 }, orientation: "h", y: 1.15 }
    };
    Plotly.newPlot("velocityBar", velocityData, velLayout);
  }

  // Render on tab activation for correct sizing
  document.getElementById("tab-temporal").addEventListener("shown.bs.tab", function(e) {
    renderCluster($("#clusterCL").val() || "all");
    var innerTab = document.querySelector("#temporal .nav-link.active");
    if (innerTab && innerTab.id) {
      if (innerTab.getAttribute("data-bs-target") === "#temp-velocity") renderVelocity();
    }
  });

  document.querySelectorAll("#temporal .nav-link").forEach(function(link) {
    link.addEventListener("shown.bs.tab", function(e) {
      if (e.target.getAttribute("data-bs-target") === "#temp-clusters") renderCluster($("#clusterCL").val() || "all");
      if (e.target.getAttribute("data-bs-target") === "#temp-velocity") renderVelocity();
    });
  });

  // DataTables
  if (document.getElementById("clusterTable")) {
    $("#clusterTable").DataTable({ pageLength: 25, searching: true, info: true,
      lengthChange: false, order: [[3, "desc"]] });
  }
  if (document.getElementById("geneActTable")) {
    var geneTable = $("#geneActTable").DataTable({
      pageLength: 25, searching: true, info: true, lengthChange: false,
      order: [] });
    $("#geneActCL").on("change", function() {
      var cl = $(this).val();
      geneTable.column(1).search(cl).draw();
    });
  }
});
</script>',
      cluster_plotly,
      velocity_bar_json
    )) else HTML(""),
    # Cross-temporal divergence JavaScript
    if (has_cross_temporal) HTML(sprintf('<script>
$(document).ready(function() {
  var ctVeloRendered = false;
  function renderCTVelocity() {
    if (ctVeloRendered) return;
    ctVeloRendered = true;
    var ctVeloData = %s;
    var ctVeloLayout = {
      title: "", xaxis: { title: "Time Point" }, yaxis: { title: "Cross-Cell-Line DEGs" },
      height: 400, margin: { l: 60, r: 20, t: 30, b: 50 },
      barmode: "group",
      legend: { font: { size: 10 }, orientation: "h", y: 1.15 }
    };
    Plotly.newPlot("ctVelocityBar", ctVeloData, ctVeloLayout);
  }
  document.getElementById("tab-divergence").addEventListener("shown.bs.tab", function(e) {
    renderCTVelocity();
  });
  if (document.getElementById("ctGeneActTable")) {
    $("#ctGeneActTable").DataTable({ pageLength: 25, searching: true, info: true, lengthChange: false, order: [] });
  }
});
</script>', ct_velocity_bar_json)) else HTML("")
  ))

# Combine into full HTML document
full_html <- sprintf('<!DOCTYPE html>\n<html lang="en">\n%s\n%s\n</html>',
                     head_html, body_html)

# Write output
out_file <- file.path(out_dir, "interactive_report.html")
cat(sprintf("  Writing %s ...\n", out_file))
writeLines(full_html, out_file)

file_size <- file.info(out_file)$size / 1e6
cat(sprintf("  Done. File size: %.1f MB\n", file_size))
cat(sprintf("  Output: %s\n", out_file))
cat("\n=== Interactive Report Complete ===\n")
