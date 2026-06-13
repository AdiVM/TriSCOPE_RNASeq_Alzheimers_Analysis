"""
SEAAD_extract_Ex_split.py

Splits SEA-AD Excitatory neurons into 3 deterministic donor-grouped MTX sets
to mirror the 427 pipeline's Ex_set1/Ex_set2/Ex_set3 structure.

Why split by donor (not cell index):
  * RUVr fits W factors on donor-level pseudobulk per set.
  * A donor straddling sets would split that donor's contribution to its
    own pseudobulk row, corrupting the RUVr factor estimation.

Why split at all:
  * Single Ex MTX is 4.4B non-zeros, exceeds R's 32-bit sparse index limit
    (~2.15B). 3 sets keeps each well below the cap.

Determinism:
  * Donor IDs are sorted lexicographically, then partitioned into 3 contiguous
    groups (sizes 27/27/26 for 80 donors).
  * A manifest is written listing every donor's set assignment.

Output (per set s in {1,2,3}):
  Ex_set{s}_counts.mtx       — raw UMI counts, cells x genes
  Ex_set{s}_barcodes.csv     — cell barcodes (row order of mtx)
  Ex_set{s}_metadata.csv     — per-cell metadata aligned to mtx rows
Plus:
  Ex_split_manifest.csv      — donor -> set_id mapping

Usage:
    python SEAAD_extract_Ex_split.py
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

CELL_TYPE     = 'Ex'
N_SETS        = 3
SUBCHUNK_SIZE = 100_000

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
    print("Loading metadata parquet...", flush=True)
    meta_df = pd.read_parquet(META_PATH)
    meta_ex = meta_df[meta_df['broad_cell_type'] == CELL_TYPE].copy()
    print(f"  Total Ex cells: {len(meta_ex):,}", flush=True)

    # ── Deterministic donor partition (sorted alphabetically -> 3 contiguous groups) ──
    donors_sorted = sorted(meta_ex['Donor ID'].unique())
    n_donors = len(donors_sorted)
    print(f"  Total Ex donors: {n_donors}", flush=True)

    # Split sizes: ceil divisions front-loaded so set1 >= set2 >= set3
    base = n_donors // N_SETS
    rem  = n_donors % N_SETS
    sizes = [base + (1 if i < rem else 0) for i in range(N_SETS)]
    print(f"  Donor split sizes: {sizes}", flush=True)

    donor_to_set = {}
    start = 0
    for s_idx, size in enumerate(sizes, start=1):
        for d in donors_sorted[start:start + size]:
            donor_to_set[d] = s_idx
        start += size

    # Write manifest BEFORE any extraction so split is reproducible / auditable.
    manifest_path = os.path.join(OUT_DIR, 'Ex_split_manifest.csv')
    pd.DataFrame({
        'donor_id': list(donor_to_set.keys()),
        'set_id':   list(donor_to_set.values()),
    }).sort_values(['set_id', 'donor_id']).to_csv(manifest_path, index=False)
    print(f"  Manifest written: {manifest_path}", flush=True)

    # ── Load h5ad backed ────────────────────────────────────────────────────────
    print("\nLoading h5ad (backed mode)...", flush=True)
    adata = sc.read_h5ad(H5AD_PATH, backed='r')
    gene_names = adata.var_names.tolist()
    n_genes = len(gene_names)
    print(f"  Genes: {n_genes:,}", flush=True)

    print("Pre-computing barcode position index...", flush=True)
    original_obs_arr = np.array(adata.obs_names)
    sort_order = np.argsort(original_obs_arr)
    sorted_obs = original_obs_arr[sort_order]

    # ── Process each Ex set ─────────────────────────────────────────────────────
    for s in range(1, N_SETS + 1):
        mtx_path      = os.path.join(OUT_DIR, f'Ex_set{s}_counts.mtx')
        barcodes_path = os.path.join(OUT_DIR, f'Ex_set{s}_barcodes.csv')
        metadata_path = os.path.join(OUT_DIR, f'Ex_set{s}_metadata.csv')

        if os.path.exists(mtx_path) and os.path.getsize(mtx_path) > 0:
            print(f"\n  [Ex set{s}] Already exists, skipping: {mtx_path}", flush=True)
            continue

        if os.path.exists(mtx_path) and os.path.getsize(mtx_path) == 0:
            print(f"\n  [Ex set{s}] Found 0-byte partial, removing.", flush=True)
            os.remove(mtx_path)

        print(f"\n{'='*60}\n  Processing Ex set {s} of {N_SETS}", flush=True)

        set_donors = [d for d, sid in donor_to_set.items() if sid == s]
        set_meta = meta_ex[meta_ex['Donor ID'].isin(set_donors)].copy()
        set_barcodes = set_meta.index.tolist()
        n_cells = len(set_barcodes)
        print(f"    Donors: {len(set_donors)}  Cells: {n_cells:,}", flush=True)

        set_barcodes_arr = np.array(set_barcodes)
        insert_positions = np.searchsorted(sorted_obs, set_barcodes_arr)
        barcode_positions = sort_order[insert_positions]

        n_subchunks = int(np.ceil(n_cells / SUBCHUNK_SIZE))
        tmp_path = mtx_path + '.tmp'
        total_nnz = 0
        row_offset = 0

        with open(tmp_path, 'wb') as f:
            f.write(b'%%MatrixMarket matrix coordinate integer general\n')
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

                counts_slice = adata.layers['UMIs'][sub_positions_sorted, :]

                if sp.issparse(counts_slice):
                    counts_slice = counts_slice[restore_order, :]
                else:
                    counts_slice = sp.csr_matrix(counts_slice[restore_order, :])

                coo = counts_slice.tocoo()
                rows = (coo.row + row_offset + 1).astype(np.int64)
                cols = (coo.col + 1).astype(np.int64)
                data = coo.data.astype(np.int64)

                np.savetxt(f, np.column_stack([rows, cols, data]), fmt='%d')

                total_nnz += coo.nnz
                row_offset += (end - start)

                del counts_slice, coo, rows, cols, data
                gc.collect()

            header = f'{n_cells} {n_genes} {total_nnz}'.ljust(63)
            f.seek(header_pos)
            f.write(header.encode('ascii'))

        os.rename(tmp_path, mtx_path)
        print(f"    MTX written: ({n_cells:,} x {n_genes:,}, {total_nnz:,} non-zeros)", flush=True)
        # sanity: confirm under R's 32-bit signed int sparse index limit
        if total_nnz >= 2**31:
            print(f"    WARNING: nnz={total_nnz} still >= 2^31; R readMM may fail.", flush=True)

        pd.DataFrame({'barcode': set_barcodes}).to_csv(barcodes_path, index=False)
        set_meta_out = set_meta.loc[set_barcodes, META_COLS_FOR_DEG].copy()
        set_meta_out.index.name = 'barcode'
        set_meta_out.reset_index().to_csv(metadata_path, index=False)

        print(f"    Saved: {mtx_path}", flush=True)
        print(f"    Saved: {barcodes_path}", flush=True)
        print(f"    Saved: {metadata_path}", flush=True)

    print("\n" + "=" * 60, flush=True)
    print("DONE. Ex split into 3 sets.", flush=True)
    print("=" * 60, flush=True)


if __name__ == '__main__':
    main()
