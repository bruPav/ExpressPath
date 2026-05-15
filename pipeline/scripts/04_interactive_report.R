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
combined <- read.delim(file.path(results_dir, "combined_results.tsv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
metadata <- read.delim(file.path(results_dir, "metadata.tsv"),
                       stringsAsFactors = FALSE)
gsea_kegg <- read.delim(file.path(out_dir, "gsea_kegg_signif.tsv"),
                        stringsAsFactors = FALSE)

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
cl_ids    <- sapply(cl_list, `[[`, "id")
tp_ids    <- sapply(tp_list, `[[`, "id")
cl_short  <- setNames(sapply(cl_list, `[[`, "short"), cl_ids)
tp_short  <- setNames(sapply(tp_list, `[[`, "short"), tp_ids)
ref_cl    <- f$reference_cell_line %||% cl_ids[1]
ref_tp    <- f$reference_time_point %||% tp_ids[1]

# Extract contrast names from combined_results columns
all_contrasts <- unique(gsub("_log2FC$", "",
  grep("_log2FC$", names(combined), value = TRUE)))
cat(sprintf("  %d contrasts found\n", length(all_contrasts)))

# Dynamic key contrast names (used in DE gene breakdown)
nonref_cl  <- setdiff(cl_ids, ref_cl)[1]
nonref_tps <- setdiff(tp_ids, ref_tp)
mock_btw_col <- paste0(nonref_cl, "_vs_", ref_cl, "_", ref_tp)
mock_padj_col <- paste0(mock_btw_col, "_padj")
mock_lfc_col  <- paste0(mock_btw_col, "_log2FC")

# Auto-generate labels and short names
contrast_labels <- setNames(sapply(all_contrasts, function(cn) {
  parts <- strsplit(cn, "_vs_")[[1]]
  if (length(parts) == 2) {
    left  <- parts[1]
    right <- parts[2]
    # Check if this is an interaction
    if (grepl("^interaction_", cn)) {
      tp <- sub("interaction_", "", cn)
      tp_label <- if (tp %in% names(tp_short)) tp else tp
      return(paste0("Interaction @", tp_label))
    }
    # Check if left is a cell_line_time and right is an interaction target
    left_parts <- strsplit(left, "_")[[1]]
    if (length(left_parts) >= 2 && left_parts[1] %in% cl_ids) {
      return(paste0(left_parts[1], " ", paste(left_parts[-1], collapse=" "),
                    " vs ", right))
    }
    return(paste0(gsub("_", " ", left), " vs ", gsub("_", " ", right)))
  }
  gsub("_", " ", cn)
}), all_contrasts)

contrast_short <- setNames(sapply(all_contrasts, function(cn) {
  if (grepl("^interaction_", cn)) {
    tp <- sub("interaction_", "", cn)
    return(paste0("I", tp_short[tp] %||% tp))
  }
  parts <- strsplit(cn, "_vs_")[[1]]
  if (length(parts) == 2) {
    left  <- parts[1]
    right <- parts[2]
    left_parts <- strsplit(left, "_")[[1]]
    if (length(left_parts) >= 2 && left_parts[1] %in% cl_ids) {
      cl <- left_parts[1]
      t2 <- paste(left_parts[-1], collapse="")
      if (right == "mock" || right == ref_tp) {
        return(paste0(cl_short[cl], tp_short[t2] %||% t2, "m"))
      }
      return(paste0(cl_short[cl], tp_short[t2] %||% t2, "v", tp_short[right] %||% right))
    }
    # Between cell lines
    if (left %in% cl_ids && right %in% cl_ids) {
      tp_label <- if (parts[2] == left) right else if (length(parts) > 2) parts[3] else "mock"
      return(paste0("B", tp_short[tp_label] %||% tp_label))
    }
  }
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

  # Baseline DE
  mock_de <- if (mock_padj_col %in% names(members))
    members[[mock_padj_col]] < 0.05 else rep(FALSE, nrow(members))
  mock_de[is.na(mock_de)] <- FALSE

  # Treatment DE (any within-cell-line time contrast)
  time_regex <- paste0("^(", paste(cl_ids, collapse = "|"), ")_(",
                      paste(nonref_tps, collapse = "|"), ")_vs_", ref_tp, "_padj$")
  time_de_cols <- grep(time_regex, names(members), value = TRUE)
  time_de <- if (length(time_de_cols) > 0)
    apply(members[, time_de_cols, drop = FALSE], 1, function(x) any(x < 0.05, na.rm = TRUE))
  else rep(FALSE, nrow(members))

  # Split into cell-line-specific for gene table
  time_cl1_cols <- grep(paste0("^", ref_cl, "_(", paste(nonref_tps, collapse="|"), ")_vs_"), names(members), value = TRUE)
  time_a549_de <- if (length(time_cl1_cols) > 0 && ref_cl %in% cl_ids)
    apply(members[, time_cl1_cols, drop = FALSE], 1, function(x) any(x < 0.05, na.rm = TRUE))
  else rep(FALSE, nrow(members))

  time_cl2_cols <- grep(paste0("^", nonref_cl, "_(", paste(nonref_tps, collapse="|"), ")_vs_"), names(members), value = TRUE)
  time_e6_de <- if (length(time_cl2_cols) > 0)
    apply(members[, time_cl2_cols, drop = FALSE], 1, function(x) any(x < 0.05, na.rm = TRUE))
  else rep(FALSE, nrow(members))

  # Between/Interaction DE
  btwn_regex <- paste0("^", nonref_cl, "_vs_", ref_cl, "_(", paste(nonref_tps, collapse = "|"),
                       ")_padj$|^interaction_(", paste(nonref_tps, collapse="|"), ")_padj$")
  btwn_de_cols <- grep(btwn_regex, names(members), value = TRUE)
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
  mock_fc <- if (mock_lfc_col %in% names(members)) members[[mock_lfc_col]][is_de] else numeric(0)
  time_fcs <- unlist(lapply(nonref_tps, function(tp) {
    col1 <- paste0(ref_cl, "_", tp, "_vs_", ref_tp, "_log2FC")
    col2 <- paste0(nonref_cl, "_", tp, "_vs_", ref_tp, "_log2FC")
    c(
      if (col1 %in% names(members)) members[[col1]][is_de] else numeric(0),
      if (col2 %in% names(members)) members[[col2]][is_de] else numeric(0)
    )
  }))
  all_fcs <- c(mock_fc, time_fcs)
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
    a <- if (time_a549_de[j]) "✓" else ""
    e <- if (time_e6_de[j]) "✓" else ""
    b <- if (btwn_de[j]) "✓" else ""
    n <- mock_de[j] + time_a549_de[j] + time_e6_de[j] + btwn_de[j]
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

    if (sum(abs(gene_fc) > 0.3) < 5) next

    setwd(pv_dir)
    tryCatch({
      suppressMessages(
        pathview(gene.data = gene_fc, pathway.id = pid, species = "hsa",
                 gene.idtype = "ENTREZ", kegg.dir = pv_cache,
                 out.suffix = suffix, limit = list(gene = 2, cpd = 1),
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
make_gene_table <- function(pid, ct) {
  core_key <- paste(pid, ct, sep = "/")
  core_entrez <- core_lookup[[core_key]]
  if (is.null(core_entrez) || length(core_entrez) == 0) return("")

  log2fc_col <- paste0(ct, "_log2FC")
  padj_col   <- paste0(ct, "_padj")
  if (!log2fc_col %in% names(combined)) return("")

  gene_rows <- c()
  for (ent in core_entrez) {
    info <- entrez_to_info[[ent]]
    if (is.null(info)) next
    rows <- which(combined$entrez == ent)
    if (length(rows) == 0) next
    r <- rows[1]

    symbol  <- combined$gene_symbol[r]
    if (is.na(symbol) || symbol == "--" || symbol == "") symbol <- combined$gene_id[r]
    lfc     <- combined[[log2fc_col]][r]
    padj    <- combined[[padj_col]][r]
    is_de   <- !is.na(padj) && padj < 0.05
    # Baseline DE status
    mock_col <- if (mock_padj_col %in% names(combined)) combined[[mock_padj_col]][r] else NA
    mock_de <- !is.na(mock_col) && mock_col < 0.05

    bold_span <- if (is_de) ' class="fw-bold"' else ""
    lfc_color <- if (!is.na(lfc)) {
      if (lfc > 0.3) "red" else if (lfc < -0.3) "blue" else "grey"
    } else "grey"
    padj_str <- if (is.na(padj)) "NA" else format(padj, digits = 2, scientific = TRUE)
    lfc_str  <- if (is.na(lfc)) "NA" else sprintf("%+.2f", lfc)
    mock_str <- if (mock_de) '<span class="text-primary">&#10003;</span>' else ""

    gene_rows <- c(gene_rows, sprintf(
      '<tr%s><td>%s</td><td style="color:%s">%s</td><td>%s</td><td>%s</td></tr>',
      bold_span, symbol, lfc_color, lfc_str, padj_str, mock_str))
    if (length(gene_rows) >= 20) break
  }
  if (length(gene_rows) == 0) return("")

  sprintf('<div class="mt-2"><small class="text-muted fw-bold">Leading-edge genes (DESeq2 padj < 0.05 shown bold, &#10003; = also DE at baseline):</small>
<table class="table table-sm table-borderless small mb-0" style="font-size:11px">
<thead><tr><th>Gene</th><th>log2FC</th><th>padj</th><th>BaseDE</th></tr></thead><tbody>%s</tbody></table></div>',
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
            sprintf("%s cell lines | %s treatment | %d genes | %d enriched KEGG pathways",
                    paste(cl_ids, collapse=" & "),
                    design$factors$treatment$label,
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
          tags$button(class = "nav-link", id = "tab-pathways", `data-bs-toggle` = "tab",
                      `data-bs-target` = "#pathways", type = "button", "Pathway Browser"))
      ),

      tags$div(class = "tab-content",
        # --- OVERVIEW TAB ---
        tags$div(class = "tab-pane fade show active", id = "overview",
          tags$div(class = "row",
            tags$div(class = "col-md-8",
              tags$h5("GSEA Enrichment (click a dot to jump to pathway)"),
              tags$div(id = "gseaDotplot", style = "height: 600px;")
            ),
            tags$div(class = "col-md-4",
              tags$h5("Enrichments per Contrast"),
              tags$div(id = "contrastBarplot", style = "height: 600px;")
            )
          )
        ),

        # --- PATHWAY BROWSER TAB ---
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
        )
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
</script>',
      jsonlite::toJSON(unname(dotplot_json), auto_unbox = TRUE, force = TRUE),
      jsonlite::toJSON(ct_summary_json, auto_unbox = TRUE, force = TRUE)
    ))
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
