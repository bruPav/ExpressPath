# RNA-Seq Time Course Analysis Pipeline

A pipeline for analyzing how gene expression changes over time in response to a treatment, across multiple cell lines.

**What it answers:**
1. Which genes respond to treatment in each cell line — and when?
2. What biological pathways are affected?
3. Do the two cell lines respond similarly or differently?

---

## How to run

```bash
cd pipeline
snakemake --use-conda -j4
```

Results appear in `results/<YYYYMMDD_HHMMSS>/`. Edit `config.yaml` to adjust thresholds.

---

## Quick start — open the interactive report

The file `results/<run>/pathway/interactive_report.html` is a self-contained web page you can open in any browser. It has three tabs:

| Tab | What it shows |
|-----|---------------|
| **Overview** | Which KEGG pathways are enriched in each comparison |
| **Pathway Browser** | Searchable pathway cards with gene-level detail, pathway maps, and stats |
| **Temporal Kinetics** | How gene expression changes over time — clusters, heatmaps, Venn diagrams |

---

## Key concepts explained (no jargon)

### What is a "differentially expressed gene" (DEG)?

A gene whose RNA level changes between two conditions in a way that is **unlikely to be noise**. The pipeline compares each treatment timepoint (1h, 3h) against the untreated mock to find these genes.

### What do the numbers mean?

| Value | What it tells you | How to read it |
|-------|-------------------|----------------|
| **log2FC** (log2 fold change) | How much the gene's RNA level changed | +1 = 2× higher, +2 = 4× higher, −1 = half as much, −2 = quarter |
| **padj** (adjusted p-value) | How confident we are the change is real | < 0.05 = 95%+ confident the gene truly changed |
| **baseMean** | Average RNA level across all samples | Higher = more abundant transcript |
| **NES** (normalized enrichment score) | For pathways: is the pathway up or down? | Positive = genes in this pathway tend to be upregulated; negative = downregulated |

### Why do some genes have huge log2FC values (e.g., +30)?

This happens when a gene produces **near-zero RNA in the mock** but measurable RNA after treatment. A gene going from "off" to "on" will mathematically produce a very large fold change — the value is inflated, but the direction and significance are real. These genes are often treatment-responsive non-coding RNAs and worth investigating.

---

## Understanding your results

Below is every file the pipeline produces, organized by what kind of question it answers.

---

## 1. Gene-level results — which genes respond?

These files are in the main results folder (e.g., `results/20260501_120000/`).

### `combined_results.tsv` — the master table

Every gene, every comparison, all statistics in one file. Open in Excel.

**Key columns:**
- `gene_id` / `gene_symbol` — gene identifier and name
- `lrt_padj` — overall: did this gene change at all across the experiment? (< 0.05 = yes)
- For each comparison (e.g., `E6_1h_vs_mock`):
  - `E6_1h_vs_mock_log2FC` — fold change of E6 at 1h relative to mock
  - `E6_1h_vs_mock_padj` — confidence that change is real
- `mock_is_DE` — TRUE if this gene differs between cell lines even before treatment (baseline difference)
- Annotation columns — functional info (GO terms, KEGG pathways, protein families, etc.)

**Use it when you want to:** Look up a specific gene, export DEG lists, or combine with other data.

### `signif_lrt.tsv`

Shorter version of combined_results: only genes that changed significantly overall (LRT padj < 0.05). This is usually the best starting point.

### `signif_lrt_foldchange1.tsv`

Filtered further: LRT-significant genes where at least one comparison has \|log2FC\| ≥ 1 (at least 2-fold change). Good for focusing on genes with clear biological effect sizes.

### `gene_activity.tsv` — the gene-by-timepoint table

One gene per row per cell line. Shows exactly which timepoints a gene is up, down, or not significant.

| Column | Meaning |
|--------|---------|
| `sig_1h` | TRUE if gene is significantly changed at 1h vs mock |
| `log2FC_1h` | Fold change at 1h (positive = up, negative = down) |
| `padj_1h` | Confidence at 1h |
| `category` | Temporal pattern: Transient, Sustained, etc. (see below) |

**Use it when you want to:** Search for genes that are up at 1h but gone by 3h, or genes that stay elevated across all timepoints.

---

## 2. Temporal patterns — when do genes respond?

### Persistence classification (`persistence_classes.tsv`)

Each DEG is classified into one of six patterns based on **which timepoints** it is significant at:

| Category | Pattern | Biological interpretation |
|----------|---------|--------------------------|
| **Transient** | Only significant at the earliest timepoint | Quick burst, then gone — early response genes |
| **Sustained** | Significant at first AND last timepoint | Long-lasting response — maintained over time |
| **Partially_Sustained** | Significant at earliest, but not latest | Fading response — starts but doesn't persist |
| **Secondary_Deferred** | Only significant at latest timepoint | Late response — needs time to kick in |
| **Transient_Mid** | Only significant at a middle timepoint | Brief mid-response (rare with only 2 timepoints) |
| **Complex** | Significant at multiple timepoints but not first→last | Mixed pattern |

These categories are visualized in the gene activity heatmaps (row-side color annotation).

### Response velocity (`velocity_summary.tsv`)

How many genes respond at each timepoint. Simple counts:

| Column | Meaning |
|--------|---------|
| `n_up` | Number of genes with higher RNA at this timepoint |
| `n_down` | Number of genes with lower RNA |
| `n_total` | Total DEGs at this timepoint |
| `mean_abs_log2FC` | Average fold-change magnitude |

**Use it when you want to:** Compare response speed between cell lines — does E6 respond faster than A549? Stronger?

### Venn diagrams (`venn_plot_<cell_line>.png`)

For each cell line: how many DEGs are shared vs unique across timepoints. Overlap = genes responding at multiple timepoints. Generated as both `.pdf` and `.png`.

The underlying data is in `venn_genelists.tsv` — a list of which gene is DE at which timepoint in which cell line.

### Gene activity heatmaps (`gene_activity_heatmap_<cell_line>.png`)

A color grid: rows = genes, columns = timepoints. Red = upregulated, blue = downregulated, white = no change. Genes are grouped by persistence category (row annotation). Color intensity is clamped to prevent a few extreme-outlier genes from washing out the scale for everything else.

**What to look for:** Blocks of red or blue across timepoints (sustained patterns), genes that switch from red to blue (complex regulation), and the overall balance of up vs down response.

---

## 3. Temporal clustering (Mfuzz)

### Mfuzz cluster profiles (`cluster_profiles_<cell_line>.pdf`)

Instead of hard categories, this groups genes by **expression shape over time** using soft clustering. Each gene belongs to every cluster with a membership score (a gene can have a "60% C1, 30% C2" pattern).

The PDF shows average expression profiles for each cluster. A cluster with a peak at 1h that drops at 3h = transient genes; a cluster that rises and stays high = sustained.

Raw data in `cluster_assignments.tsv` (gene → best cluster + membership scores) and `cluster_mean_profiles.tsv` (average shape of each cluster).

---

## 4. Pathway analysis — what biological processes are affected?

These files are in `results/<run>/pathway/`.

### GSEA — which KEGG pathways are enriched?

**Traditional analysis:** For each comparison, takes the full gene list (ranked by fold change), then asks: "are genes in this pathway clustered at the top of the list?" This detects subtle, coordinated changes even when individual genes have modest fold changes.

| File | Content |
|------|---------|
| `gsea_kegg_signif.tsv` | KEGG pathways with padj < 0.05 |
| `gsea_go_signif.tsv` | GO Biological Process terms with padj < 0.05 |

**Key columns in GSEA results:**
- `NES` — normalized enrichment score: positive = pathway is upregulated, negative = downregulated
- `padj` — confidence
- `core_enrichment` — the specific genes driving the enrichment
- `contrast` — which comparison this is for (e.g., `E6_1h_vs_mock`)

### GSVA — pathway activity per sample

**Different approach:** For each sample, calculates a single score representing the overall activity of a pathway (e.g., "How active is Apoptosis in E6 at 3h?"). Then tests whether pathway activity changes between conditions.

| File | Content |
|------|---------|
| `gsva_scores.tsv` | Pathway activity score for every pathway × every sample |
| `gsva_diff_results.tsv` | Which pathways have significantly different activity between conditions |

### Pathway visualizations

| File | What it shows |
|------|---------------|
| `gsea_dotplot_kegg.pdf` | Dot chart: one row per pathway, size = significance, color = direction, grouped by comparison |
| `gsva_diff_dotplot.pdf` | Same format but for GSVA pathway activity differences |
| `gsva_heatmap.pdf` | Color grid of pathway activity across all samples (rows = pathways, columns = samples) |
| `pathview_output/<pathway_id>/` | KEGG pathway maps with individual genes colored by fold change (red = up, blue = down) |

**How to read a pathview map:** Green boxes = gene products in the pathway. Red fill = gene is upregulated in this comparison, blue = downregulated, grey = no change/no data. White boxes are not measured genes. This shows you exactly where in the pathway your hits land.

---

## 5. Quality control plots — check these first

These help you verify the data quality before diving into results.

| File | What to look for |
|------|-----------------|
| `pca_plot.pdf` | Dots should cluster by cell line and timepoint, not by batch. If batch dominates, there may be a batch effect. |
| `pca_batch_plot.pdf` | Same PCA colored by batch — should see no strong separation by batch. |
| `sample_distance_heatmap.pdf` | Replicates of the same condition should form dark blocks (high similarity). Different conditions should show lighter colors (less similar). |
| `heatmap_top50.pdf` | Top 50 most-changed genes across all samples. Should show clear patterns, not random noise. |
| `volcano_<contrast>.pdf` | Each dot = a gene. X-axis = log fold change (how big the change is). Y-axis = significance (how confident we are). Red dots = significant genes. Look for roughly symmetric distribution (similar number of up and down-regulated genes). |

---

## 6. All output files — quick reference

### In `results/<run>/`

| File | What it contains | Format |
|------|-----------------|--------|
| `combined_results.tsv` | Master table: every gene × every comparison with stats and annotations | TSV |
| `signif_lrt.tsv` | Subset: only genes significantly changed overall (LRT padj < 0.05) | TSV |
| `signif_lrt_foldchange1.tsv` | Subset: LRT-sig genes with \|log2FC\| ≥ 1 in at least one comparison | TSV |
| `vst_normalized_counts.tsv` | Normalized expression values for every gene × sample (VST-transformed) | TSV |
| `gene_activity.tsv` | Per-gene, per-timepoint significance, fold change, and pattern category | TSV |
| `persistence_classes.tsv` | Temporal pattern classification for each DEG (Transient, Sustained, etc.) | TSV |
| `venn_genelists.tsv` | Which genes are DE at which timepoint (feeds Venn diagrams) | TSV |
| `velocity_summary.tsv` | DEG counts and mean fold change per timepoint per cell line | TSV |
| `cluster_assignments.tsv` | Mfuzz fuzzy cluster membership for each LRT-sig gene | TSV |
| `cluster_mean_profiles.tsv` | Average expression profile of each Mfuzz cluster | TSV |
| `counts_matrix.tsv` | Raw integer RNA counts (genes × samples) | TSV |
| `metadata.tsv` | Sample annotations (cell line, timepoint, replicate, batch) | TSV |
| `gene_annotations.tsv` | Functional annotations from eggNOG, GO, KEGG, Pfam, etc. | TSV |
| `pca_plot.pdf` | PCA colored by condition | PDF |
| `pca_batch_plot.pdf` | PCA colored by batch | PDF |
| `sample_distance_heatmap.pdf` | Sample-to-sample similarity matrix | PDF |
| `heatmap_top50.pdf` | Expression heatmap of top 50 DEGs | PDF |
| `volcano_<contrast>.pdf` | Volcano plot for key comparisons (up to 6) | PDF |
| `velocity_barplot.pdf` | Bar chart: DEG counts per timepoint | PDF |
| `velocity_fc_boxplot.pdf` | Boxplot: fold change distribution per timepoint | PDF |
| `cluster_profiles_<cell_line>.pdf` | Mfuzz cluster expression profiles (one per cell line) | PDF |
| `venn_plot_<cell_line>.pdf/.png` | Venn diagram: DEG overlap across timepoints | PDF+PNG |
| `gene_activity_heatmap_<cell_line>.pdf/.png` | Heatmap of log2FC by gene and timepoint, grouped by category | PDF+PNG |

### In `results/<run>/pathway/`

| File | What it contains | Format |
|------|-----------------|--------|
| `gsea_kegg_signif.tsv` | Significant KEGG pathway enrichments per contrast | TSV |
| `gsea_go_signif.tsv` | Significant GO Biological Process enrichments per contrast | TSV |
| `gsva_scores.tsv` | Hallmark pathway activity scores per sample | TSV |
| `gsva_diff_results.tsv` | Differential pathway activity between conditions | TSV |
| `gsea_dotplot_kegg.pdf` | Dot plot of GSEA KEGG results | PDF |
| `gsva_diff_dotplot.pdf` | Dot plot of GSVA differential results | PDF |
| `gsva_heatmap.pdf` | Heatmap of GSVA pathway activity | PDF |
| `interactive_report.html` | Self-contained interactive HTML report | HTML |
| `pathview_output/<pathway_id>/` | KEGG pathway maps with genes colored by fold change | PNG |

---

## 7. Key analysis decisions

| Decision | What we do | Why |
|----------|-----------|-----|
| Significance threshold | padj < 0.05 | Standard threshold: 95% confidence the change is real |
| Fold change estimation | `lfcShrink` with `ashr` method | Shrinks unreliable fold changes (low-count genes) toward zero while preserving strong signals |
| Gene filtering | Remove genes with < 10 total counts across all samples | Very low-count genes are unreliable |
| Temporal clustering | Mfuzz soft clustering on LRT-significant genes | Genes can belong to multiple patterns — reflects biology better than hard clusters |
| Persistence classification | Based on which timepoints a gene is significant at (padj < 0.05, no fold-change filter) | Tracks set membership, not magnitude — even a small but consistent response across timepoints counts |
| Heatmap color scale | Clamped at the 90th percentile of \|log2FC\| per cell line, with a floor of ±3 | Prevents a few extreme-outlier genes from washing out the color scale; the top 10% saturate at the color boundary |
| Replicates | Replicates handled by DESeq2's built-in dispersion model | Accounts for biological variability within each condition |

---

## 8. Configuration

Edit `pipeline/config.yaml` to adjust:

| Setting | Default | What it does |
|---------|---------|-------------|
| `pvalue_threshold` | 0.05 | Significance cutoff for DEG calling |
| `count_filter_threshold` | 10 | Minimum total reads across all samples to keep a gene |
| `temporal_n_clusters` | 6 | Number of Mfuzz clusters (or "auto") |
| `temporal_foldchange_threshold` | 0 | Additional \|log2FC\| filter for persistence/Venn (0 = disabled, use all significant genes) |
| `gsea_min_set_size` | 15 | Minimum genes a pathway must have to be tested |
| `gsea_max_set_size` | 500 | Maximum genes a pathway can have to be tested |

---

## 9. Troubleshooting

**"I changed config.yaml but nothing changed"** — Snakemake needs `configfile: "config.yaml"` in the Snakefile to read it (included by default). Run `snakemake --forceall --use-conda` to force a full re-run.

**"The heatmap looks all white except a few genes"** — This is clamped at P90 \|log2FC\| to prevent extreme outliers from compressing the color scale. The strongly-colored genes are the top 10% by fold change; the white ones still have real changes but at smaller magnitude. You can adjust the clamping in `02_deseq2_analysis.R` (lines 901–905).

**"Some genes have log2FC = 30"** — These are genes with near-zero RNA in the mock condition. The fold change number is mathematically huge but the biological reality is "gene turned on from off." Check their baseMean and raw counts to verify expression levels.

**"I want to add more cell lines or timepoints"** — Edit `data/design.yaml` and regenerate with `setup_design.html`. The pipeline auto-discovers contrasts from the design.
