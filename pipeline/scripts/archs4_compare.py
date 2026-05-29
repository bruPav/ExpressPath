#!/usr/bin/env python3
"""
Public Transcriptome Comparison - Download relevant public RNA-seq datasets
from GEO and compare with your data via PCA.

Downloads supplementary count files from known GEO datasets (lightweight,
~1-5 MB each) and runs PCA to show where your transcriptomes fall relative
to public reference samples.

Usage:
    python pipeline/scripts/archs4_compare.py
    python pipeline/scripts/archs4_compare.py --input data/test_samplesheet.tsv
    python pipeline/scripts/archs4_compare.py --out-dir results/public_comparison
    python pipeline/scripts/archs4_compare.py --add-gse GSE12345,GSE67890
"""

import argparse
import csv
import gzip
import io
import os
import sys
import warnings
from collections import defaultdict

import numpy as np
import requests
import yaml

warnings.filterwarnings("ignore")

GEO_CACHE = "data/geo_cache"

DEFAULT_REF_DATASETS = [
    {
        "gse": "GTEx_v8",
        "label": "GTEx tissues (54)",
        "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
        "desc": "54 human tissues, median TPM across ~1000 donors",
        "format": "gct",
        "gene_col": 1,
        "skip_lines": 3,
        "data_type": "tpm",
    },
    {
        "gse": "GSE147507",
        "label": "A549 + viruses (SARS2/RSV/IAV/HPIV3)",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147507/suppl/GSE147507_RawReadCounts_Human.tsv.gz",
        "desc": "A549 & NHBE cells infected with SARS-CoV-2, RSV, IAV, HPIV3",
    },
    {
        "gse": "GSE173484",
        "label": "A549 treatments",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173484/suppl/GSE173484_allCounts_RNA-seq.txt.gz",
        "desc": "A549 cells with DOX-induced treatments + controls",
    },
]

FALLBACK_REF_DATASETS = [
    {
        "gse": "GSE162104",
        "label": "A549 SARS-CoV-2 timecourse",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162104/suppl/GSE162104_count.txt.gz",
        "desc": "A549 infected with SARS-CoV-2, time course",
    },
]


def parse_tsv(input_tsv, design):
    """Parse the TSV file into gene-level expression dict."""
    column_map = design.get("column_map", {})
    if not column_map:
        print("ERROR: design.yaml has no column_map section.", file=sys.stderr)
        sys.exit(1)

    count_columns = []
    for col_name, info in column_map.items():
        sid = info.get("sample_id", f"{info['cell_line']}_{info['time']}_{info['replicate']}")
        batch = info.get("batch", f"batch_{info['replicate']}")
        count_columns.append((col_name, sid, info["cell_line"], info["time"],
                              str(info["replicate"]), batch))

    col_names = [c[0] for c in count_columns]
    sample_ids = [c[1] for c in count_columns]

    with open(input_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = list(reader.fieldnames)

    gene_id_col = header[0]

    gene_to_symbol = {}
    counts = defaultdict(lambda: defaultdict(int))

    with open(input_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = row.get(gene_id_col, "").strip()
            symbol = row.get("gene_name", "").strip()
            if not gene_id or not symbol:
                continue
            gene_to_symbol[gene_id] = symbol
            for cn, sid in zip(col_names, sample_ids):
                val = row.get(cn, "0").strip()
                try:
                    counts[gene_id][sid] = int(float(val))
                except (ValueError, TypeError):
                    counts[gene_id][sid] = 0

    sample_meta = []
    for col_name, sid, cl, tp, rep, batch in count_columns:
        sample_meta.append({
            "sample_id": sid, "cell_line": cl,
            "time": tp, "replicate": rep, "batch": batch,
        })

    return counts, gene_to_symbol, sample_meta, sample_ids


def build_count_matrix(counts, genes_ordered, sample_ids):
    """Build gene x sample count matrix."""
    n_genes = len(genes_ordered)
    n_samples = len(sample_ids)
    matrix = np.zeros((n_genes, n_samples), dtype=np.float64)
    for i, gene_id in enumerate(genes_ordered):
        for j, sid in enumerate(sample_ids):
            matrix[i, j] = counts[gene_id].get(sid, 0)
    return matrix


def normalize_log2_cpm(matrix):
    """Compute log2(CPM + 1) normalization."""
    lib_sizes = matrix.sum(axis=0)
    lib_sizes[lib_sizes == 0] = 1
    cpm = matrix / lib_sizes * 1e6
    return np.log2(cpm + 1)


def download_suppl_counts(url, cache_dir=GEO_CACHE, gene_col=0, skip_lines=1):
    """Download and parse a supplementary count file.

    For standard TSV (gene_col=0, skip_lines=1): gene IDs in col 0,
    sample data in cols 1+.

    For GCT format (gene_col=1, skip_lines=3): GTEx-style, 3 metadata
    header rows, gene symbol in col 1, tissue data in cols 2+.

    Returns (matrix, gene_symbols, sample_ids) or None.
    """
    filename = os.path.basename(url).replace(".gz", "")
    cache_file = os.path.join(cache_dir, filename + ".npy")
    cache_meta = os.path.join(cache_dir, filename + "_meta.json")

    if os.path.exists(cache_file) and os.path.exists(cache_meta):
        import json
        matrix = np.load(cache_file)
        with open(cache_meta) as f:
            meta = json.load(f)
        return matrix, meta["genes"], meta["samples"]

    try:
        resp = requests.get(url, timeout=120)
        resp.raise_for_status()
    except Exception as e:
        print(f"    Download failed: {e}", file=sys.stderr)
        return None

    try:
        if url.endswith(".gz"):
            with gzip.GzipFile(fileobj=io.BytesIO(resp.content)) as f:
                text = f.read().decode("utf-8", errors="replace")
        else:
            text = resp.content.decode("utf-8", errors="replace")
    except Exception as e:
        print(f"    Decompress failed: {e}", file=sys.stderr)
        return None

    all_lines = text.split("\n")
    if len(all_lines) < skip_lines + 1:
        print(f"    File too short: {len(all_lines)} lines", file=sys.stderr)
        return None

    header_line_idx = skip_lines - 1
    header = all_lines[header_line_idx].strip().split("\t")

    if gene_col >= len(header):
        print(f"    gene_col={gene_col} out of range (header has {len(header)} cols)",
              file=sys.stderr)
        return None

    sample_names = header[gene_col + 1:]
    gene_symbols = []
    matrix_rows = []

    for line in all_lines[skip_lines:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) <= gene_col:
            continue
        gene = parts[gene_col].strip()
        if not gene or gene in ("", "N/A", "NA", "."):
            continue
        vals = []
        data_start = gene_col + 1
        for v in parts[data_start:]:
            try:
                vals.append(float(v))
            except (ValueError, TypeError):
                vals.append(0.0)
        while len(vals) < len(sample_names):
            vals.append(0.0)
        vals = vals[:len(sample_names)]
        gene_symbols.append(gene)
        matrix_rows.append(vals)

    if not gene_symbols:
        print(f"    No valid gene rows found", file=sys.stderr)
        return None

    matrix = np.array(matrix_rows, dtype=np.float64)

    os.makedirs(cache_dir, exist_ok=True)
    np.save(cache_file, matrix)
    import json
    with open(cache_meta, "w") as f:
        json.dump({"genes": gene_symbols, "samples": sample_names}, f)

    return matrix, gene_symbols, sample_names


def auto_discover_suppl(gse_id, cache_dir=GEO_CACHE):
    """Try to auto-discover supplementary count files from a GSE record."""
    try:
        import GEOparse
        gse = GEOparse.get_GEO(geo=gse_id, destdir=cache_dir, silent=True)
    except Exception as e:
        print(f"    GEO download failed for {gse_id}: {e}", file=sys.stderr)
        return None

    suppl_files = []
    for key in gse.metadata:
        if "suppl" in key.lower():
            files = gse.metadata[key]
            if isinstance(files, list):
                suppl_files.extend(files)
            else:
                suppl_files.append(str(files))

    count_keywords = ["count", "rawcount", "readcount", "expression", "rna-seq", "rnaseq"]
    count_files = []
    for f in suppl_files:
        fl = f.lower()
        if any(kw in fl for kw in count_keywords):
            count_files.append(f)

    score_file = None
    best_score = 0
    for cf in count_files:
        score = 0
        if "count" in cf.lower():
            score += 10
        if "raw" in cf.lower():
            score += 5
        if "human" in cf.lower():
            score += 3
        if score > best_score:
            best_score = score
            score_file = cf

    if score_file:
        url = score_file.replace("ftp://", "https://")
        print(f"    Found supplementary file: {os.path.basename(score_file)}")
        return download_suppl_counts(url, cache_dir)

    return None


def run_pca(user_matrix, user_genes, user_sample_ids, user_meta,
            ref_matrices, ref_genes_list, ref_sample_ids_list, ref_labels):
    """Run PCA on combined dataset with gene-symbol matching."""
    user_gene_to_idx = {}
    for i, g in enumerate(user_genes):
        sym = g.strip().upper()
        if sym:
            user_gene_to_idx[sym] = i

    ref_gene_to_idx_list = []
    ref_gene_sets = []
    for ref_genes in ref_genes_list:
        g_to_i = {}
        for i, g in enumerate(ref_genes):
            sym = g.strip().upper()
            if sym:
                g_to_i[sym] = i
        ref_gene_to_idx_list.append(g_to_i)
        ref_gene_sets.append(set(g_to_i.keys()))

    common_genes = set(user_gene_to_idx.keys())
    for gs in ref_gene_sets:
        common_genes &= gs

    if not common_genes:
        print("ERROR: No common genes found!", file=sys.stderr)
        return None

    gene_list = sorted(common_genes)
    n_genes = len(gene_list)
    print(f"  Common genes: {n_genes}")

    combined = np.zeros((n_genes, user_matrix.shape[1]), dtype=np.float64)
    for i, gene in enumerate(gene_list):
        combined[i, :] = user_matrix[user_gene_to_idx[gene], :]

    all_sample_ids = list(user_sample_ids)
    all_sample_labels = ["your_data"] * len(user_sample_ids)
    all_sample_groups = [
        f"{sm['cell_line']}_{sm['time']}" for sm in user_meta
    ]
    all_dataset_labels = ["your_data"] * len(user_sample_ids)

    for k, (ref_mat_raw, ref_genes, ref_sids, ref_label) in enumerate(
        zip(ref_matrices, ref_genes_list, ref_sample_ids_list, ref_labels)
    ):
        g2i = ref_gene_to_idx_list[k]
        n_ref = len(ref_sids)
        ref_block_unnorm = np.zeros((n_genes, n_ref), dtype=np.float64)
        for i, gene in enumerate(gene_list):
            if gene in g2i:
                ref_block_unnorm[i, :] = ref_mat_raw[g2i[gene], :]
        ref_block = normalize_log2_cpm(ref_block_unnorm)
        combined = np.hstack([combined, ref_block])
        all_sample_ids.extend(ref_sids)
        all_sample_labels.extend(["public_data"] * n_ref)
        all_sample_groups.extend([ref_label] * n_ref)
        all_dataset_labels.extend([ref_label] * n_ref)

    # Z-score normalize per gene
    row_means = combined.mean(axis=1, keepdims=True)
    row_stds = combined.std(axis=1, keepdims=True)
    row_stds[row_stds == 0] = 1
    combined_norm = (combined - row_means) / row_stds

    centered = combined_norm - combined_norm.mean(axis=1, keepdims=True)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    n_samps = combined.shape[1]
    pc1 = Vt[0, :n_samps] * S[0]
    pc2 = Vt[1, :n_samps] * S[1]
    var_pc1 = S[0]**2 / (S**2).sum() * 100
    var_pc2 = S[1]**2 / (S**2).sum() * 100

    return pc1, pc2, var_pc1, var_pc2, all_sample_labels, all_sample_groups, all_dataset_labels


def plot_pca(pc1, pc2, var1, var2, labels, groups, dataset_labels, out_path):
    """Generate PCA scatter plot with cell-line colors and timepoint intensity."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(14, 10))

    is_your = np.array([l == "your_data" for l in labels])
    is_public = ~is_your

    # --- Your samples: one color per cell line, alpha per timepoint ---
    cell_line_colors = {"A549": "#2166ac", "E6": "#d6604d"}
    time_alpha = {"mock": 0.30, "1h": 0.60, "3h": 0.95}
    time_order = {"mock": 1, "1h": 2, "3h": 3}

    your_groups_raw = [g for g, y in zip(groups, is_your) if y]
    your_labels_raw = [l for l, y in zip(labels, is_your) if y]
    your_pc1 = pc1[is_your]
    your_pc2 = pc2[is_your]

    for cl in ["A549", "E6"]:
        for tp in ["mock", "1h", "3h"]:
            mask = np.array([
                g.startswith(cl) and tp in g and lbl == "your_data"
                for g, lbl in zip(groups, labels)
            ])
            if not mask.any():
                continue
            color = cell_line_colors.get(cl, "#888888")
            alpha = time_alpha.get(tp, 0.5)
            ax.scatter(pc1[mask], pc2[mask], c=color, alpha=alpha,
                       s=60, edgecolors="black", linewidths=1.0, zorder=6)

    # Build custom legend entries for your samples
    legend_entries = []
    for cl in ["A549", "E6"]:
        color = cell_line_colors.get(cl, "#888888")
        for tp in ["mock", "1h", "3h"]:
            alpha = time_alpha.get(tp, 0.5)
            label = f"{cl} {tp}"
            legend_entries.append(
                Line2D([0], [0], marker="o", color="w", markerfacecolor=color,
                       markersize=8, alpha=alpha, markeredgecolor="black",
                       markeredgewidth=0.8, label=label)
            )

    # --- Divider ---
    legend_entries.append(
        Line2D([0], [0], marker="", color="w", label="─" * 20)
    )

    # --- Public datasets ---
    public_datasets = sorted(set(d for d, p in zip(dataset_labels, is_public) if p))
    public_colors = [
        "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf", "#999999"
    ]
    for i, ds in enumerate(public_datasets):
        mask = np.array([d == ds and l == "public_data" for d, l in zip(dataset_labels, labels)])
        color = public_colors[i % len(public_colors)]
        ax.scatter(pc1[mask], pc2[mask], c=color, label=ds,
                   s=35, alpha=0.60, edgecolors="white", linewidths=0.3, zorder=4)
        legend_entries.append(
            Line2D([0], [0], marker="o", color="w", markerfacecolor=color,
                   markersize=8, alpha=0.60, markeredgecolor="white",
                   markeredgewidth=0.3, label=ds)
        )

    ax.set_xlabel(f"PC1 ({var1:.1f}%)", fontsize=14)
    ax.set_ylabel(f"PC2 ({var2:.1f}%)", fontsize=14)
    ax.set_title("PCA: Your Transcriptomes vs Public Reference Datasets",
                 fontsize=16, fontweight="bold")
    ax.legend(handles=legend_entries, bbox_to_anchor=(1.02, 1),
              loc="upper left", fontsize=9, title="Samples", title_fontsize=10,
              framealpha=0.9)
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.5, alpha=0.4)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.5, alpha=0.4)
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  PCA plot saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare your transcriptomes to public GEO RNA-seq data via PCA"
    )
    parser.add_argument("--input", default="data/test_samplesheet.tsv",
                        help="Path to input TSV count file")
    parser.add_argument("--design", default="data/design.yaml",
                        help="Path to design.yaml")
    parser.add_argument("--out-dir", default="results/public_comparison",
                        help="Output directory")
    parser.add_argument("--add-gse", default="",
                        help="Comma-separated extra GSE IDs to try")
    parser.add_argument("--no-download", action="store_true",
                        help="Skip downloading, use cached data only")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: {args.input} not found", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.design):
        print(f"ERROR: {args.design} not found", file=sys.stderr)
        sys.exit(1)

    with open(args.design) as f:
        design = yaml.safe_load(f)

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(GEO_CACHE, exist_ok=True)

    print("=" * 80)
    print("  Public Transcriptome Comparison (GEO)")
    print("=" * 80)
    print(f"  Input:    {args.input}")
    print(f"  Design:   {args.design}")
    print(f"  Output:   {args.out_dir}")

    # 1. Parse user data
    print("\n[1/4] Parsing your count data...")
    counts, gene_to_symbol, sample_meta, sample_ids = parse_tsv(args.input, design)
    genes_ordered = list(gene_to_symbol.keys())
    user_symbols = [gene_to_symbol[g] for g in genes_ordered]
    user_matrix = build_count_matrix(counts, genes_ordered, sample_ids)
    print(f"  Genes: {len(genes_ordered)}, Samples: {len(sample_ids)}")

    # 2. Normalize user data
    print("[2/4] Normalizing your data (log2 CPM)...")
    user_norm = normalize_log2_cpm(user_matrix)

    # 3. Download reference datasets
    ref_datasets = list(DEFAULT_REF_DATASETS)
    if args.add_gse:
        for extra_gse in args.add_gse.split(","):
            extra_gse = extra_gse.strip()
            if extra_gse:
                ref_datasets.append({
                    "gse": extra_gse, "label": extra_gse,
                    "url": None, "desc": "User-specified dataset",
                })

    ref_matrices = []
    ref_gene_symbols = []
    ref_sample_ids_list = []
    ref_labels = []

    print("[3/4] Downloading public reference data...")
    n_loaded = 0

    for ds in ref_datasets:
        print(f"  [{ds['gse']}] {ds['label']}")
        print(f"          {ds['desc']}")

        result = None
        if ds.get("url"):
            result = download_suppl_counts(
                ds["url"],
                gene_col=ds.get("gene_col", 0),
                skip_lines=ds.get("skip_lines", 1),
            )
        else:
            if not args.no_download:
                result = auto_discover_suppl(ds["gse"])

        if result is None:
            print(f"    Failed. Trying auto-discovery...")
            result = auto_discover_suppl(ds["gse"])

        if result is not None:
            matrix, genes, sids = result
            ref_matrices.append(matrix)
            ref_gene_symbols.append(genes)
            ref_sample_ids_list.append(sids)
            ref_labels.append(ds["label"])
            n_loaded += 1
            print(f"    Loaded: {len(genes)} genes x {len(sids)} samples")
        else:
            print(f"    Skipped: could not download/parse")

    if n_loaded == 0:
        print("\n  ERROR: No reference datasets loaded. Check internet connection", file=sys.stderr)
        print("  or try different GSE IDs with --add-gse.", file=sys.stderr)
        return

    # 4. PCA
    print(f"\n[4/4] Running PCA ({n_loaded} reference datasets)...")
    pca_result = run_pca(
        user_norm, user_symbols, sample_ids, sample_meta,
        ref_matrices, ref_gene_symbols, ref_sample_ids_list, ref_labels,
    )

    if pca_result is None:
        print("PCA failed: no common genes between datasets.", file=sys.stderr)
        return

    pc1, pc2, var1, var2, all_labels, all_groups, all_datasets = pca_result

    out_plot = os.path.join(args.out_dir, "pca_public_comparison.png")
    plot_pca(pc1, pc2, var1, var2, all_labels, all_groups, all_datasets, out_plot)

    # Summary
    summary_path = os.path.join(args.out_dir, "comparison_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Public Transcriptome Comparison Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Your data: {len(sample_ids)} samples\n")
        f.write(f"Reference datasets: {n_loaded}\n")
        for lbl in ref_labels:
            f.write(f"  - {lbl}\n")
        f.write(f"\nPCA: PC1={var1:.1f}%, PC2={var2:.1f}%\n")
    print(f"  Summary: {summary_path}")

    print(f"\n{'=' * 80}")
    print("  Done! PCA plot ready.")
    print(f"  Your samples (colored by cell_line + timepoint) are shown")
    print(f"  alongside public reference datasets.")
    print(f"\n  Output: {out_plot}")

    # Interpretation hints
    print(f"\n  Interpretation:")
    print(f"  - If your samples cluster WITH the public samples of a category,")
    print(f"    your transcriptome resembles that tissue/treatment.")
    print(f"  - If your samples cluster SEPARATELY, your transcriptome may have")
    print(f"    unique features or novel response patterns.")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()
