# ExpressPath

RNA-seq time course analysis pipeline — DESeq2 differential expression →
GSEA pathway enrichment → Pathview KEGG maps → interactive HTML report.

Metadata-driven. Single Snakemake command from raw counts to browsable results.

## Quick Start

1. Open `pipeline/setup_design.html` in your browser
2. Define cell lines, time points, treatment, replicates
3. Browse for your TSV file → map columns to samples
4. Download `design.yaml` → save to `data/`
5. Put your TSV file in `data/` (same name as shown in the GUI)
6. `cd pipeline && snakemake -j4 --use-conda`
7. Open `results/<timestamp>/pathway/interactive_report.html`

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
| `pathway/gsea_kegg_signif.tsv` | Enriched KEGG pathways (GSEA) |
| `pathway/gsea_go_signif.tsv` | Enriched GO terms (GSEA) |
| `pathway/gsva_scores.tsv` | Per-sample pathway activity scores |
| `pathway/pathview_output/` | KEGG pathway maps with log2FC overlay |
| `pathway/interactive_report.html` | Self-contained browsable report |

## How It Works

```
data/design.yaml  +  data/your_data.tsv
        │
   [extract_counts]      Python — reads column_map from design
        │
   [deseq2_analysis]     R/DESeq2 — LRT + pairwise Wald contrasts
        │
   [pathway_analysis]    R/clusterProfiler — GSEA + Pathview + GSVA
        │
   [interactive_report]  R/htmltools — self-contained HTML report
```

Contrasts are auto-generated from your experiment design — no hardcoded cell
line or time point names. Add more time points or rename cell lines in
`data/design.yaml` and everything adapts.

## Configuration

Edit `pipeline/config.yaml` for pipeline parameters (filters, thresholds).
Edit `data/design.yaml` for your experiment (cell lines, time points, column
mapping, batch labels). Both can be set with `setup_design.html`.

## License

MIT — see [LICENSE](LICENSE)

## Authors

Dr. Bruno Pavletić — [bruno.pavletic@irb.hr](mailto:bruno.pavletic@irb.hr)

Ruđer Bošković Institute, Zagreb, Croatia
