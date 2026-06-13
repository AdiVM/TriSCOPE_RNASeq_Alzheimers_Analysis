"""
SEAAD_extract_matrix_only.py

Extracts and normalizes the gene expression matrix from the SEA-AD h5ad,
writing one parquet file per cell type. Each cell type is processed in
sub-chunks of SUBCHUNK_SIZE cells to stay under 200GB RAM.

Metadata parquet is assumed to already exist from the previous run.

Output files (one per cell type):
  SEAAD_Matrix_Ex.parquet
  SEAAD_Matrix_Inh.parquet
  SEAAD_Matrix_Ast.parquet
  SEAAD_Matrix_Oli.parquet
  SEAAD_Matrix_Mic.parquet
  SEAAD_Matrix_Opc.parquet

Normalization: CPM -> log2(x+1)  (matches original Morabito pipeline)

Usage:
    python SEAAD_extract_matrix_only.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import pyarrow as pa
import pyarrow.parquet as pq
import os
import gc

# ── Paths ──────────────────────────────────────────────────────────────────────
H5AD_PATH  = '/n/scratch/users/a/adm808/SEA_AD/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad'
OUT_DIR    = '/n/groups/patel/adithya/SEAAD_Outputs/'
META_PATH  = os.path.join(OUT_DIR, 'SEAAD_CellMetadata.parquet')

os.makedirs(OUT_DIR, exist_ok=True)

CELL_TYPES_TO_KEEP = ['Ex', 'Inh', 'Ast', 'Oli', 'Mic', 'Opc']

# 100K cells x 30K genes x 4 bytes = ~12GB per sub-chunk — well within 200GB
SUBCHUNK_SIZE = 100_000


def normalize_cpm_log2(counts_matrix):
    """
    CPM -> log2(x + 1) normalization.
    Input:  scipy sparse or dense numpy array, shape (cells, genes)
    Output: dense numpy float32 array, same shape
    """
    if sp.issparse(counts_matrix):
        counts_matrix = counts_matrix.toarray()
    counts_matrix = counts_matrix.astype(np.float32)

    total_counts = counts_matrix.sum(axis=1, keepdims=True)
    total_counts = np.where(total_counts == 0, 1, total_counts)

    cpm = counts_matrix / total_counts * 1e6
    log2cpm = np.log2(cpm + 1)

    return log2cpm.astype(np.float32)


def main():
    # ── Load metadata (already saved) ─────────────────────────────────────────
    print("Loading metadata parquet...")
    meta_df = pd.read_parquet(META_PATH)
    print(f"  Metadata shape: {meta_df.shape}")

    # ── Load h5ad in backed mode ───────────────────────────────────────────────
    print("\nLoading h5ad (backed mode)...")
    adata = sc.read_h5ad(H5AD_PATH, backed='r')
    gene_names = adata.var_names.tolist()
    print(f"  Total cells in h5ad: {adata.n_obs:,}")
    print(f"  Total genes:         {len(gene_names):,}")

    # ── Pre-compute barcode position index once ────────────────────────────────
    print("\nPre-computing barcode position index...")
    original_obs_arr = np.array(adata.obs_names)
    sort_order = np.argsort(original_obs_arr)
    sorted_obs = original_obs_arr[sort_order]
    print(f"  Index built over {len(original_obs_arr):,} barcodes.")

    # ── Process one cell type at a time ───────────────────────────────────────
    for ct in CELL_TYPES_TO_KEEP:
        out_path = os.path.join(OUT_DIR, f'SEAAD_Matrix_{ct}.parquet')

        # Skip if already done (safe to re-run)
        if os.path.exists(out_path):
            print(f"\n  [{ct}] Already exists, skipping: {out_path}")
            continue

        print(f"\n{'='*60}")
        print(f"  Processing cell type: {ct}")

        # All barcodes for this cell type
        ct_barcodes = meta_df[meta_df['broad_cell_type'] == ct].index.tolist()
        n_cells = len(ct_barcodes)
        print(f"    Total cells: {n_cells:,}")
        print(f"    Sub-chunk size: {SUBCHUNK_SIZE:,}")
        print(f"    Number of sub-chunks: {int(np.ceil(n_cells / SUBCHUNK_SIZE))}")

        # Fast positional lookup for ALL barcodes of this cell type
        ct_barcodes_arr = np.array(ct_barcodes)
        insert_positions = np.searchsorted(sorted_obs, ct_barcodes_arr)
        barcode_positions = sort_order[insert_positions]

        # PyArrow writer — appends each sub-chunk to the same parquet file
        writer = None
        schema = None

        n_subchunks = int(np.ceil(n_cells / SUBCHUNK_SIZE))

        for i in range(n_subchunks):
            start = i * SUBCHUNK_SIZE
            end   = min(start + SUBCHUNK_SIZE, n_cells)
            print(f"    Sub-chunk {i+1}/{n_subchunks}: cells {start:,} - {end:,}")

            sub_barcodes  = ct_barcodes[start:end]
            sub_positions = barcode_positions[start:end]

            # Sort for sequential HDF5 disk access
            disk_order            = np.argsort(sub_positions)
            sub_positions_sorted  = sub_positions[disk_order]
            restore_order         = np.argsort(disk_order)

            # Pull from disk
            counts_slice = adata.layers['UMIs'][sub_positions_sorted, :]

            # Restore original barcode order
            counts_slice = counts_slice[restore_order, :]

            # Normalize
            norm_slice = normalize_cpm_log2(counts_slice)
            del counts_slice
            gc.collect()

            # Build DataFrame
            chunk_df = pd.DataFrame(
                norm_slice,
                index=sub_barcodes,
                columns=gene_names,
                dtype=np.float32
            )
            chunk_df.index.name = 'TAG'
            del norm_slice
            gc.collect()

            # Convert to Arrow table and append to parquet file
            table = pa.Table.from_pandas(chunk_df)
            del chunk_df
            gc.collect()

            if writer is None:
                schema = table.schema
                writer = pq.ParquetWriter(out_path, schema, compression='snappy')

            writer.write_table(table)
            del table
            gc.collect()

            print(f"      Written sub-chunk {i+1}/{n_subchunks}")

        if writer is not None:
            writer.close()

        size_gb = os.path.getsize(out_path) / 1e9
        print(f"    Saved: {out_path}  ({size_gb:.2f} GB)")

    # ── Final check ───────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("EXTRACTION COMPLETE — files written:")
    for ct in CELL_TYPES_TO_KEEP:
        out_path = os.path.join(OUT_DIR, f'SEAAD_Matrix_{ct}.parquet')
        if os.path.exists(out_path):
            size_gb = os.path.getsize(out_path) / 1e9
            print(f"  {out_path}  ({size_gb:.2f} GB)")
        else:
            print(f"  MISSING: {out_path}")
    print("=" * 60)


if __name__ == '__main__':
    main()