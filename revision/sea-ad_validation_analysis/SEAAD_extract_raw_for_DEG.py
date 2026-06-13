"""
SEAAD_extract_raw_for_DEG.py

Extracts RAW UMI counts per cell type from SEA-AD A9 h5ad, writes sparse MTX
files + per-cell-type metadata CSVs for consumption by R DEG pipeline.

One cell type at a time, streamed subchunk-by-subchunk directly to disk
to stay under RAM limits (no in-memory accumulation of the full matrix).

Output (per cell type):
  {cell_type}_counts.mtx       — sparse raw UMI counts, cells x genes
  {cell_type}_barcodes.csv     — cell barcodes (rows of mtx)
  {cell_type}_genes.csv        — gene names (cols of mtx)
  {cell_type}_metadata.csv     — per-cell metadata for DE model

Usage:
    python SEAAD_extract_raw_for_DEG.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os
import gc

# ── Paths ──────────────────────────────────────────────────────────────────────
H5AD_PATH  = '/n/scratch/users/a/adm808/SEA_AD/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad'
META_PATH  = '/n/groups/patel/adithya/SEAAD_Outputs/SEAAD_CellMetadata.parquet'
OUT_DIR    = '/n/scratch/users/a/adm808/SEAAD_DEG_input/'

os.makedirs(OUT_DIR, exist_ok=True)

CELL_TYPES_TO_KEEP = ['Ast', 'Oli', 'Mic', 'Inh', 'Ex']
# 'Ex', 'Inh', 'Ast', 'Oli', 'Mic', 'Opc'
SUBCHUNK_SIZE = 100_000

# Metadata columns needed by the R DEG pipeline
# (match the glmmTMB model: counts ~ AD + sex + age + pmi + (1|donor) + offset(log(libsize)))
META_COLS_FOR_DEG = [
    'Donor ID',
    'Sex',
    'Age at Death',
    'PMI',
    'APOE Genotype',
    'alzheimers_or_control',
    'broad_cell_type',
]


def main():
    # ── Load metadata ─────────────────────────────────────────────────────────
    print("Loading metadata parquet...", flush=True)
    meta_df = pd.read_parquet(META_PATH)
    print(f"  Metadata shape: {meta_df.shape}", flush=True)

    # ── Load h5ad in backed mode ───────────────────────────────────────────────
    print("\nLoading h5ad (backed mode)...", flush=True)
    adata = sc.read_h5ad(H5AD_PATH, backed='r')
    gene_names = adata.var_names.tolist()
    n_genes = len(gene_names)
    print(f"  Total cells in h5ad: {adata.n_obs:,}", flush=True)
    print(f"  Total genes:         {n_genes:,}", flush=True)

    # ── Pre-compute barcode position index ─────────────────────────────────────
    print("\nPre-computing barcode position index...", flush=True)
    original_obs_arr = np.array(adata.obs_names)
    sort_order = np.argsort(original_obs_arr)
    sorted_obs = original_obs_arr[sort_order]

    # Save gene names ONCE (same for every cell type)
    genes_path = os.path.join(OUT_DIR, 'genes.csv')
    pd.DataFrame({'gene': gene_names}).to_csv(genes_path, index=False)
    print(f"  Saved gene list: {genes_path}", flush=True)

    # ── Process one cell type at a time ────────────────────────────────────────
    for ct in CELL_TYPES_TO_KEEP:
        mtx_path      = os.path.join(OUT_DIR, f'{ct}_counts.mtx')
        barcodes_path = os.path.join(OUT_DIR, f'{ct}_barcodes.csv')
        metadata_path = os.path.join(OUT_DIR, f'{ct}_metadata.csv')

        if os.path.exists(mtx_path) and os.path.getsize(mtx_path) > 0:
            print(f"\n  [{ct}] Already exists, skipping: {mtx_path}", flush=True)
            continue

        # Clean up any 0-byte partial from a previous crash
        if os.path.exists(mtx_path) and os.path.getsize(mtx_path) == 0:
            print(f"\n  [{ct}] Found 0-byte partial, removing: {mtx_path}", flush=True)
            os.remove(mtx_path)

        print(f"\n{'='*60}", flush=True)
        print(f"  Processing cell type: {ct}", flush=True)

        ct_meta = meta_df[meta_df['broad_cell_type'] == ct].copy()
        ct_barcodes = ct_meta.index.tolist()
        n_cells = len(ct_barcodes)
        print(f"    Total cells: {n_cells:,}", flush=True)

        ct_barcodes_arr = np.array(ct_barcodes)
        insert_positions = np.searchsorted(sorted_obs, ct_barcodes_arr)
        barcode_positions = sort_order[insert_positions]

        n_subchunks = int(np.ceil(n_cells / SUBCHUNK_SIZE))

        # ── Stream-write MTX directly to disk ──────────────────────────────────
        # Never accumulate full matrix in RAM. Reserve header space, write
        # coordinate triplets per subchunk, rewrite header with true nnz at end.
        tmp_path = mtx_path + '.tmp'
        total_nnz = 0
        row_offset = 0

        with open(tmp_path, 'wb') as f:
            f.write(b'%%MatrixMarket matrix coordinate integer general\n')
            # Reserve 63 bytes for "n_rows n_cols n_nnz\n"; rewrite once known
            header_pos = f.tell()
            f.write(b' ' * 63 + b'\n')

            for i in range(n_subchunks):
                start = i * SUBCHUNK_SIZE
                end   = min(start + SUBCHUNK_SIZE, n_cells)
                print(f"    Sub-chunk {i+1}/{n_subchunks}: cells {start:,} - {end:,}", flush=True)

                sub_positions = barcode_positions[start:end]
                disk_order           = np.argsort(sub_positions)
                sub_positions_sorted = sub_positions[disk_order]
                restore_order        = np.argsort(disk_order)

                # Pull raw UMI counts from disk (sequential access)
                counts_slice = adata.layers['UMIs'][sub_positions_sorted, :]

                # Restore original barcode order, ensure sparse
                if sp.issparse(counts_slice):
                    counts_slice = counts_slice[restore_order, :]
                else:
                    counts_slice = sp.csr_matrix(counts_slice[restore_order, :])

                # Convert to COO and stream triplets to disk (mtx is 1-indexed)
                coo = counts_slice.tocoo()
                rows = (coo.row + row_offset + 1).astype(np.int64)
                cols = (coo.col + 1).astype(np.int64)
                data = coo.data.astype(np.int64)

                np.savetxt(f, np.column_stack([rows, cols, data]), fmt='%d')

                total_nnz += coo.nnz
                row_offset += (end - start)

                del counts_slice, coo, rows, cols, data
                gc.collect()

            # Rewrite header with actual dims (space-padded to fill reserved area)
            header = f'{n_cells} {n_genes} {total_nnz}'
            header = header.ljust(63)
            f.seek(header_pos)
            f.write(header.encode('ascii'))

        # Atomic rename — skip logic won't see partial files
        os.rename(tmp_path, mtx_path)
        print(f"    MTX written: ({n_cells:,} x {n_genes:,}, {total_nnz:,} non-zeros)", flush=True)

        # Barcodes aligned to MTX rows
        pd.DataFrame({'barcode': ct_barcodes}).to_csv(barcodes_path, index=False)

        # Per-cell metadata aligned to MTX rows (same order as ct_barcodes)
        ct_meta_out = ct_meta.loc[ct_barcodes, META_COLS_FOR_DEG].copy()
        ct_meta_out.index.name = 'barcode'
        ct_meta_out.reset_index().to_csv(metadata_path, index=False)

        print(f"    Saved: {mtx_path}", flush=True)
        print(f"    Saved: {barcodes_path}", flush=True)
        print(f"    Saved: {metadata_path}", flush=True)

    print("\n" + "=" * 60, flush=True)
    print("DONE. Files ready for R DEG pipeline.", flush=True)
    print("=" * 60, flush=True)


if __name__ == '__main__':
    main()