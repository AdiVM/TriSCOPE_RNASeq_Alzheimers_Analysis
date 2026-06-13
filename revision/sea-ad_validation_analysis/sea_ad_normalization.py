"""
sea_ad_normalization.py

Produces two parquet files from the SEA-AD h5ad:
  1. SEAAD_CellMetadata.parquet   — one row per cell, key columns only
  2. SEAAD_NormalizedCellMatrix.parquet — cells x genes, CPM -> log2(x+1)
     processed in 6 chunks (one per broad cell type) to keep RAM manageable.

Normalization matches original Morabito pipeline:
  CPM = counts / total_counts_per_cell * 1e6
  X   = log2(CPM + 1)

Run on cluster with at least 64GB RAM recommended.
Usage:
    python sea_ad_normalization.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os
import gc

# ── Paths ──────────────────────────────────────────────────────────────────────
H5AD_PATH    = '/n/scratch/users/a/adm808/SEA_AD/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad'
OUT_DIR      = '/n/groups/patel/adithya/SEAAD_Outputs/'
META_OUT     = os.path.join(OUT_DIR, 'SEAAD_CellMetadata.parquet')
MATRIX_OUT   = os.path.join(OUT_DIR, 'SEAAD_NormalizedCellMatrix.parquet')

os.makedirs(OUT_DIR, exist_ok=True)

# ── Cell type mapping: Subclass -> broad_cell_type ─────────────────────────────
CELLTYPE_MAP = {
    # Excitatory neurons
    'L2/3 IT'   : 'Ex',
    'L4 IT'     : 'Ex',
    'L5 IT'     : 'Ex',
    'L5 ET'     : 'Ex',
    'L6 IT'     : 'Ex',
    'L6 IT Car3': 'Ex',
    'L6 CT'     : 'Ex',
    'L6b'       : 'Ex',
    'L5/6 NP'   : 'Ex',
    # Inhibitory neurons
    'Pvalb'     : 'Inh',
    'Sst'       : 'Inh',
    'Vip'       : 'Inh',
    'Lamp5'     : 'Inh',
    'Sncg'      : 'Inh',
    'Lamp5 Lhx6': 'Inh',
    'Sst Chodl' : 'Inh',
    'Chandelier' : 'Inh',
    'Pax6'      : 'Inh',
    # Glia / non-neuronal
    'Astrocyte'     : 'Ast',
    'Oligodendrocyte': 'Oli',
    'Microglia-PVM' : 'Mic',
    'OPC'           : 'Opc',
    # Excluded — not mapped
    'Endothelial': None,
    'VLMC'       : None,
}

# Broad cell types we will process (one chunk each)
CELL_TYPES_TO_KEEP = ['Ex', 'Inh', 'Ast', 'Oli', 'Mic', 'Opc']

# Metadata columns to retain in the metadata parquet
META_COLS = [
    'Donor ID',
    'Sex',
    'Age at Death',
    'PMI',
    'APOE Genotype',
    'Cognitive Status',
    'Overall AD neuropathological Change',
    'Subclass',           # kept for reference / QC
    'Neurotypical reference',
]


def normalize_cpm_log2(counts_matrix):
    """
    CPM -> log2(x + 1) normalization.
    Input:  scipy sparse or dense numpy array, shape (cells, genes)
    Output: dense numpy float32 array, same shape
    """
    if sp.issparse(counts_matrix):
        counts_matrix = counts_matrix.toarray()
    counts_matrix = counts_matrix.astype(np.float32)

    # Per-cell total counts
    total_counts = counts_matrix.sum(axis=1, keepdims=True)

    # Avoid divide-by-zero for cells with 0 total counts (shouldn't happen post-QC)
    total_counts = np.where(total_counts == 0, 1, total_counts)

    # CPM
    cpm = counts_matrix / total_counts * 1e6

    # log2(CPM + 1)
    log2cpm = np.log2(cpm + 1)

    return log2cpm.astype(np.float32)


def main():
    print("=" * 60)
    print("Loading h5ad (backed mode — expression matrix stays on disk)")
    print("=" * 60)
    adata = sc.read_h5ad(H5AD_PATH, backed='r')
    obs_df = adata.obs.copy()
    gene_names = adata.var_names.tolist()

    print(f"Total cells loaded: {len(obs_df):,}")
    print(f"Total genes:        {len(gene_names):,}")

    # ── Step 1: Build broad_cell_type column ──────────────────────────────────
    print("\nBuilding broad_cell_type from Subclass...")
    obs_df['broad_cell_type'] = obs_df['Subclass'].map(CELLTYPE_MAP)

    # ── Step 2: Drop Reference donors and unmapped cell types ─────────────────
    print("Filtering out Reference donors and unmapped cell types...")
    n_before = len(obs_df)

    obs_df = obs_df[obs_df['Cognitive Status'] != 'Reference'].copy()
    obs_df = obs_df[obs_df['broad_cell_type'].isin(CELL_TYPES_TO_KEEP)].copy()

    n_after = len(obs_df)
    print(f"  Cells before filter: {n_before:,}")
    print(f"  Cells after filter:  {n_after:,}")

    # ── Step 3: Build binary label ────────────────────────────────────────────
    obs_df['alzheimers_or_control'] = (
        obs_df['Cognitive Status'].eq('Dementia').astype(int)
    )

    # ── Step 4: Build and save metadata parquet ───────────────────────────────
    print("\nSaving metadata parquet...")

    # TAG = cell barcode (obs_df index) — this is the join key with the matrix
    meta_out_df = obs_df[META_COLS + ['broad_cell_type', 'alzheimers_or_control']].copy()
    meta_out_df.index.name = 'TAG'   # align with original pipeline convention

    # Verify no duplicate TAGs
    assert meta_out_df.index.nunique() == len(meta_out_df), \
        "ERROR: Duplicate cell barcodes found in metadata index!"

    meta_out_df.to_parquet(META_OUT)
    print(f"  Saved: {META_OUT}")
    print(f"  Shape: {meta_out_df.shape}")
    print(f"\n  Diagnosis breakdown (donor-level):")
    donor_diag = meta_out_df.groupby('Donor ID')['alzheimers_or_control'].first()
    print(f"    AD (1):      {donor_diag.sum()} donors")
    print(f"    Control (0): {(donor_diag == 0).sum()} donors")
    print(f"\n  Cell type breakdown:")
    print(meta_out_df['broad_cell_type'].value_counts().to_string())

    # ── Step 5: Normalize and save expression matrix, chunk by cell type ──────
    print("\n" + "=" * 60)
    print("Extracting and normalizing expression matrix (UMIs layer)")
    print("Chunking by cell type to manage RAM")
    print("=" * 60)

    # Pre-compute sorted index for fast barcode position lookup (done once, reused per chunk)
    print("Pre-computing barcode position index...")
    original_obs_arr = np.array(adata.obs_names)
    sort_order = np.argsort(original_obs_arr)
    sorted_obs = original_obs_arr[sort_order]
    print(f"  Index built over {len(original_obs_arr):,} barcodes.")

    all_chunks = []

    for ct in CELL_TYPES_TO_KEEP:
        print(f"\n  Processing cell type: {ct}")

        # Get barcodes for this cell type
        ct_barcodes = meta_out_df[meta_out_df['broad_cell_type'] == ct].index.tolist()
        print(f"    Cells: {len(ct_barcodes):,}")

        # Get integer positions using pre-computed sorted index
        ct_barcodes_arr = np.array(ct_barcodes)
        insert_positions = np.searchsorted(sorted_obs, ct_barcodes_arr)
        barcode_positions = sort_order[insert_positions]

        # Sort positions before disk read — HDF5 backed mode is much faster
        # with sequential (sorted) integer indices; random seeks are very slow.
        # Track the reorder so we can restore original barcode order after.
        disk_order = np.argsort(barcode_positions)
        barcode_positions_sorted = barcode_positions[disk_order]
        restore_order = np.argsort(disk_order)  # inverts disk_order

        # Slice UMIs layer — sequential access on disk
        print(f"    Loading UMIs from disk (sorted sequential access)...")
        counts_slice = adata.layers['UMIs'][barcode_positions_sorted, :]

        # Restore original barcode order so index aligns with ct_barcodes
        counts_slice = counts_slice[restore_order, :]

        # Normalize
        print(f"    Normalizing (CPM -> log2(x+1))...")
        norm_slice = normalize_cpm_log2(counts_slice)

        # Build DataFrame
        chunk_df = pd.DataFrame(
            norm_slice,
            index=ct_barcodes,
            columns=gene_names,
            dtype=np.float32
        )
        chunk_df.index.name = 'TAG'

        all_chunks.append(chunk_df)
        print(f"    Done. Shape: {chunk_df.shape}")

        # Free memory before next chunk
        del counts_slice, norm_slice, chunk_df
        gc.collect()

    # ── Step 6: Concatenate all chunks and save ───────────────────────────────
    print("\nConcatenating all cell type chunks...")
    full_matrix = pd.concat(all_chunks, axis=0)
    del all_chunks
    gc.collect()

    print(f"Final matrix shape: {full_matrix.shape}")

    # Verify TAG alignment with metadata
    assert set(full_matrix.index) == set(meta_out_df.index), \
        "ERROR: Matrix TAGs and metadata TAGs do not match!"

    print(f"\nSaving matrix parquet to: {MATRIX_OUT}")
    full_matrix.to_parquet(MATRIX_OUT, compression='snappy')
    print("Done.")

    print("\n" + "=" * 60)
    print("EXTRACTION COMPLETE")
    print(f"  Metadata : {META_OUT}")
    print(f"  Matrix   : {MATRIX_OUT}")
    print("=" * 60)


if __name__ == '__main__':
    main()