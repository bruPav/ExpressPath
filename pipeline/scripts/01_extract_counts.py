#!/usr/bin/env python3
"""
Extract count matrix, metadata, and gene annotations from an RNA-Seq samples TSV.
Handles naming inconsistencies via configurable column renames.
Supports optional batch metadata file for sample_id → batch mapping.
"""

import csv
import os
import sys
import json

# ─── Snakemake integration ───
try:
    from snakemake.script import snakemake
    input_tsv  = snakemake.input["tsv"]
    out_dir    = os.path.dirname(snakemake.output["counts"])
    count_filter = snakemake.params.get("count_filter", 10)
    design = json.loads(snakemake.params.get("design", "{}"))
except ImportError:
    # Standalone mode — read from command line or use defaults
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument("input_tsv")
    parser.add_argument("--out-dir", default="results")
    parser.add_argument("--count-filter", type=int, default=10)
    parser.add_argument("--design", default="../data/design.yaml")
    args = parser.parse_args()
    input_tsv   = args.input_tsv
    out_dir     = args.out_dir
    count_filter = args.count_filter
    with open(args.design) as f:
        design = yaml.safe_load(f)

os.makedirs(out_dir, exist_ok=True)

# ─── Read column_map from design ───
column_map = design.get("column_map", {})
if not column_map:
    print("ERROR: design.yaml has no column_map section. Run setup_design.html first.",
          file=sys.stderr)
    sys.exit(1)

# Build COUNT_COLUMNS from explicit map
COUNT_COLUMNS = []
for col_name, info in column_map.items():
    sid   = info.get("sample_id", f"{info['cell_line']}_{info['time']}_{info['replicate']}")
    batch = info.get("batch", f"batch_{info['replicate']}")
    COUNT_COLUMNS.append((
        col_name,
        sid,
        info["cell_line"],
        info["time"],
        str(info["replicate"]),
        batch
    ))

sample_ids = [c[1] for c in COUNT_COLUMNS]

# Annotation columns to preserve
ANNOTATION_COLUMNS = [
    "COG_class", "COG_class_annotation", "GO_annotation", "KEGG_annotation",
    "KEGG_pathway_annotation", "KOG_class", "KOG_class_annotation",
    "Pfam_annotation", "Swiss_Prot_annotation", "eggNOG_class",
    "eggNOG_class_annotation", "NR_annotation", "TF_family",
    "GO_second_level_annotation",
]


def main():
    # Read header
    with open(input_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = list(reader.fieldnames)

    # Verify all count columns from column_map exist in TSV
    count_col_names = [c[0] for c in COUNT_COLUMNS]
    missing = [c for c in count_col_names if c not in header]
    if missing:
        print(f"ERROR: Missing count columns in TSV: {missing}", file=sys.stderr)
        sys.exit(1)

    # Build metadata rows
    metadata_rows = []
    for col_name, sample_id, cell_line, time_point, rep, batch in COUNT_COLUMNS:
        metadata_rows.append({
            "sample_id": sample_id,
            "cell_line": cell_line,
            "time": time_point,
            "replicate": rep,
            "batch": batch,
        })

    # Verify annotation columns
    missing_annot = [c for c in ANNOTATION_COLUMNS if c not in header]
    if missing_annot:
        print(f"WARNING: Missing annotation columns: {missing_annot}", file=sys.stderr)
    annot_cols = [c for c in ANNOTATION_COLUMNS if c in header]

    # Parse all rows
    counts_data = {}
    gene_annot  = {}
    gene_symbol = {}
    gene_order  = []

    with open(input_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = row.get("", "").strip()
            symbol  = row.get("gene_name", "").strip()

            if not gene_id:
                continue

            gene_order.append(gene_id)
            gene_symbol[gene_id] = symbol

            counts_data[gene_id] = {}
            gene_annot[gene_id] = {}
            for cn in annot_cols:
                gene_annot[gene_id][cn] = row.get(cn, "").strip()

            for col_name, sid in zip(count_col_names, sample_ids):
                val = row.get(col_name, "0").strip()
                try:
                    count = int(float(val))
                except (ValueError, TypeError):
                    count = 0
                counts_data[gene_id][sid] = count

    print(f"Parsed {len(gene_order)} genes")
    print(f"Found {len(sample_ids)} samples")

    # Filter low-count genes
    n_filtered = 0
    kept_genes = []
    for gene_id in gene_order:
        total = sum(counts_data[gene_id].get(sid, 0) for sid in sample_ids)
        if total >= count_filter:
            kept_genes.append(gene_id)
        else:
            n_filtered += 1

    print(f"Filtered {n_filtered} low-count genes (total count < {count_filter})")
    print(f"Kept {len(kept_genes)} genes")

    # Write counts matrix
    counts_path = os.path.join(out_dir, "counts_matrix.tsv")
    with open(counts_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gene_id"] + sample_ids)
        for gene_id in kept_genes:
            row = [gene_id] + [counts_data[gene_id].get(sid, 0) for sid in sample_ids]
            writer.writerow(row)
    print(f"Wrote counts matrix -> {counts_path} ({len(kept_genes)} genes x {len(sample_ids)} samples)")

    # Write metadata
    meta_path = os.path.join(out_dir, "metadata.tsv")
    with open(meta_path, "w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t",
                                fieldnames=["sample_id", "cell_line", "time", "replicate", "batch"])
        writer.writeheader()
        writer.writerows(metadata_rows)
    print(f"Wrote metadata -> {meta_path}")

    # Write gene annotations
    annot_path = os.path.join(out_dir, "gene_annotations.tsv")
    with open(annot_path, "w", newline="") as f:
        fields = ["gene_id", "gene_symbol"] + annot_cols
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for gene_id in kept_genes:
            row = {"gene_id": gene_id, "gene_symbol": gene_symbol.get(gene_id, "")}
            row.update(gene_annot.get(gene_id, {}))
            writer.writerow(row)
    print(f"Wrote gene annotations -> {annot_path}")

    print("Done.")


if __name__ == "__main__":
    main()
