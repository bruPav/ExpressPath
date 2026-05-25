# ExpressPath

RNA-seq time course analysis pipeline — DESeq2 differential expression →
temporal & cross-cell-line classification → GSEA pathway enrichment →
TF enrichment → Pathview KEGG maps → interactive HTML report.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/powered%20by-snakemake-green.svg)](https://snakemake.github.io)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-4.3-0092AC.svg)](https://bioconductor.org)

Metadata-driven. Single Snakemake command from raw counts to browsable results.

## Quick Start

1. Open `pipeline/setup_design.html` in your browser
2. Define cell lines, time points, treatment, replicates
3. Browse for your TSV file → map columns to samples
4. Download `design.yaml` → save to `data/`
5. Put your TSV file in `data/` (same name as shown in the GUI)
6. Run `./run.sh` (add `-j8` for more cores, `-n` for dry-run)
7. Open the report URL printed when the pipeline finishes

## Requirements

- conda (or mamba)
- snakemake (`conda install -c bioconda -c conda-forge snakemake`)

## Input

Tab-separated file with gene expression counts. Required columns:

| Column | Header | Example |
|--------|--------|---------|
| Gene ID | *(unnamed — first column)* | `ENSG00000000003` |
| Gene name | `gene_name` | `TSPAN6` |
| Count data | *any names* | `A549_mock_1_count`, `Sample_A`, … |

Your count columns can have any names — you map them to samples in the setup GUI.
Optional annotation columns (GO, KEGG, COG, etc.) are carried through to outputs
if present.

See `data/design.example.yaml` for the experiment configuration format.

## Output

All in `results/<timestamp>/`:

| File | Description |
|------|-------------|
| `combined_results.tsv` | All genes — LRT p-values + all pairwise log2FC and padj |
| `signif_lrt.tsv` | Genes with LRT padj < 0.05 |
| `counts_matrix.tsv` | Filtered count matrix (DESeq2 input) |
| `vst_normalized_counts.tsv` | VST-transformed counts |
| `cross_temporal/persistence_classes.tsv` | Per-cell-line temporal persistence categories |
| `cross_temporal/gene_activity.tsv` | Per-gene log2FC and significance at each timepoint |
| `cross_temporal/venn_genelists.tsv` | DEG sets per cell line × timepoint (for Venn/UpSet) |
| `cross_temporal/cross_cellline_shared.tsv` | Genes DE in both cell lines at each timepoint (concordance + magnitude divergence) |
| `cross_temporal/cross_cellline_specific.tsv` | Genes DE in only one cell line at each timepoint |
| `cross_cellline/cross_temporal_persistence.tsv` | Cross-cell-line temporal divergence categories |
| `cross_cellline/cross_temporal_gene_activity.tsv` | Per-gene between-cell-line log2FC at each timepoint |
| `cross_cellline/cross_temporal_shared.tsv` | Between-cell-line DEGs shared across timepoints |
| `cross_cellline/cross_temporal_specific.tsv` | Between-cell-line DEGs specific to one timepoint |
| `pathway/gsea_kegg_signif.tsv` | Enriched KEGG pathways (GSEA) |
| `pathway/gsea_go_signif.tsv` | Enriched GO terms (GSEA) |
| `pathway/gsva_scores.tsv` | Per-sample pathway activity scores |
| `pathway/pathview_output/` | KEGG pathway maps with log2FC overlay |
| `pathway/interactive_report.html` | Self-contained browsable report |
| `tf/tf_enrichment_results.tsv` | TF target enrichment (enrichR) — all contrasts |
| `tf/tf_enrichment_heatmap.pdf` | Heatmap of TF enrichment significance per analysis dimension |
| `tf/tf_regulatory_network_*.html` | Per-cell-line / shared / divergence TF–target regulatory networks |

## How It Works

```
data/design.yaml  +  data/your_data.tsv
        │
   [extract_counts]      Python — reads column_map from design
        │
   [deseq2_analysis]     R/DESeq2 — LRT + pairwise Wald contrasts
        │                                     │
        │                     ┌────────────────┘
        │                     ▼
        │              [temporal_mfuzz]  R/Mfuzz — clustering, persistence,
        │                                     cross-cell-line comparisons
        │
   [pathway_analysis]    R/clusterProfiler — GSEA + Pathview + GSVA
        │
   [tf_enrichment]       R/enrichR — TF target enrichment + regulatory networks
        │
   [interactive_report]  R/htmltools — self-contained HTML report
```

Contrasts are auto-generated from your experiment design — no hardcoded cell
line or time point names. Add more time points or rename cell lines in
`data/design.yaml` and everything adapts.

## Analysis Categories

The pipeline classifies DEGs into temporal activity categories at three levels:

### Per-Cell-Line Temporal Persistence

For each cell line, genes are classified by which treatment timepoints they are
DE at (vs mock). Output in `cross_temporal/persistence_classes.tsv`.

| Category | Meaning |
|----------|---------|
| `Transient` | DE only at the first treatment timepoint |
| `Transient_Mid` | DE at a single intermediate timepoint |
| `Secondary_Deferred` | DE only at the last treatment timepoint |
| `Sustained` | DE at the first AND last treatment timepoints, with contiguous significance across all intermediate timepoints |
| `Partially_Sustained` | DE contiguously from the first through an intermediate timepoint, but NOT at the last |
| `Intermittent` | DE at the first AND last treatment timepoints, but with gaps (non-contiguous) |
| `Complex` | Any other multi-timepoint pattern not fitting the above |

### Cross-Cell-Line at Each Timepoint

Compares DEG sets between two cell lines at each treatment timepoint.
Output in `cross_temporal/cross_cellline_shared.tsv` and
`cross_temporal/cross_cellline_specific.tsv`.

| Category | Meaning |
|----------|---------|
| `Concordant_Up` | Both cell lines upregulated (same direction) |
| `Concordant_Down` | Both cell lines downregulated (same direction) |
| `Discordant` | One up, one down (labeled as `{CL}_Up_{CL}_Down`) |
| `Magnitude Divergent` | Shared gene where \|log2FC ratio\| between cell lines > 2 |
| `Cell-line-specific` | DE in one cell line only (absent from the other) |

### Cross-Cell-Line Temporal Divergence (Part G)

Classifies how the *between-cell-line difference* evolves over time. Uses the
cell-line-vs-cell-line contrasts at each timepoint (e.g. E6 vs A549 at mock,
1h, 3h). Output in `cross_cellline/cross_temporal_persistence.tsv`.

| Category | Meaning |
|----------|---------|
| `Constitutive` | Between-cell-line difference significant at ALL timepoints |
| `Baseline_Only` | Difference only at the reference timepoint (pre-existing, disappears after treatment) |
| `Emergent_Early` | Difference appears only at the first treatment timepoint |
| `Emergent_Mid` | Difference appears at a single mid-treatment timepoint |
| `Emergent_Late` | Difference appears only at the last treatment timepoint |
| `Emergent_Sustained` | Difference at ALL treatment timepoints but NOT at baseline — treatment-induced and persistent |
| `Emergent_Complex` | Difference at multiple (but not all) treatment timepoints, not at baseline |
| `Convergent` | Difference present at baseline but ABSENT by the last timepoint (cell lines become more similar) |
| `Complex` | Any other multi-timepoint pattern |

### DESeq2 Contrasts

The statistical model `~ batch + cell_line + time + cell_line:time` produces four
types of pairwise Wald contrasts, each gated by flags in `data/design.yaml`:

| Comparison type | Flag | Example | What it tests |
|----------------|------|---------|---------------|
| Within cell line | `within_cell_line` | `A549_3h_vs_mock` | Treatment effect per cell line |
| Progression | `progression` | `A549_3h_vs_1h` | Response evolution between consecutive timepoints |
| Between cell lines | `between_cell_lines` | `E6_vs_A549_1h` | Cell line difference at each timepoint |
| Interaction | `interactions` | `interaction_3h` | Does the time effect differ between cell lines? |

An LRT (omnibus) test on `~ batch + cell_line + time + cell_line:time` vs.
`~ batch + cell_line` identifies genes that change in any way across the
experiment.

## Configuration

Edit `pipeline/config.yaml` for pipeline parameters (filters, thresholds).
Edit `data/design.yaml` for your experiment (cell lines, time points, column
mapping, batch labels). Both can be set with `setup_design.html`.

## License

MIT — see [LICENSE](LICENSE)

## Authors

Dr. Bruno Pavletić — [bruno.pavletic@irb.hr](mailto:bruno.pavletic@irb.hr)

Ruđer Bošković Institute, Zagreb, Croatia
