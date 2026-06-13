"""
label_seaad_subpops_centroid.py

Stage A of the SEA-AD -> Mathys subpopulation labeling pipeline, centroid
correlation variant.

This is a sibling of label_seaad_subpops.py with one and only one methodological
change: per-Mathys-subpop centroid (mean log-normalized expression across the
cells defining that subpop in Mathys) is computed from the Mathys reference
matrix, and each SEA-AD cell is assigned to the broad-matched subpop whose
centroid has the highest Pearson correlation with the cell's expression
vector over the genes shared between datasets.

Motivation: see the RESULTS LOG block in label_seaad_subpops.py. The
sc.tl.score_genes implementation was degenerate for Ex (99.2% of cells went
to Ex2 because Ex2's 52 markers are pan-expressed glutamatergic-synaptic
genes in SEA-AD even though they are uniquely-Ex2 within Mathys). Pearson
correlation is scale-invariant and uses the entire shared transcriptome,
so cells are matched on overall expression *profile* rather than mean
expression of a small marker set.

Everything else is intentionally identical to label_seaad_subpops.py:
  - SEA-AD matrix loader (pyarrow streaming into pre-allocated float32 array)
  - SEA-AD metadata join on TAG with strict reindex check
  - top1 / top2 / margin computation via argpartition + stable sort
  - Output parquet schema (Donor ID, Subclass, broad_cell_type,
    alzheimers_or_control, score_<subpop>..., top1, top1_score, top2,
    top2_score, margin)
  - SLURM submission pattern (same per-broad mem/partition budgets)
  - CLI argument names where they overlap

The only output difference is the filename: results land in
  {output_dir}/seaad_subpop_centroid_{broad}.parquet
so the score_genes parquets and the centroid parquets can coexist.

Usage:
    python label_seaad_subpops_centroid.py --broad_cell_type Ex

Inputs (defaults match the project layout):
  --matrix_dir     dir containing SEAAD_Matrix_{broad}.parquet (cells x genes, log2(CPM+1))
  --meta_path      SEAAD_CellMetadata.parquet (TAG index, has broad_cell_type and Subclass)
  --mathys_matrix  Mathys log-normalized expression parquet (genes x cells)
  --mathys_meta    Mathys cell metadata parquet (has TAG and Subcluster columns)
"""

import argparse
import gc
import os
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


# SEA-AD broad_cell_type -> Mathys subpop name prefix.
# SEA-AD uses "Inh" for inhibitory but Mathys uses "In*" (In0..In11).
BROAD_TO_PREFIX = {
    'Ex':  'Ex',
    'Inh': 'In',
    'Ast': 'Ast',
    'Oli': 'Oli',
    'Mic': 'Mic',
    'Opc': 'Opc',
}

DEFAULT_MATRIX_DIR    = '/n/groups/patel/adithya/SEAAD_Outputs'
DEFAULT_META_PATH     = '/n/groups/patel/adithya/SEAAD_Outputs/SEAAD_CellMetadata.parquet'
DEFAULT_MATHYS_MATRIX = '/home/adm808/NormalizedCellMatrixSyn18485175.parquet'
DEFAULT_MATHYS_META   = '/home/adm808/New_CellMetadataSyn1848517.parquet'
# Distinct from the score_genes labels/ directory so the two methods' outputs
# can be inspected side-by-side without overwriting each other.
DEFAULT_OUTPUT_DIR    = '/n/groups/patel/adithya/Alz_Outputs/Final_Outputs/seaad_pseudotime/labels_centroid'

# ---- Reproducibility ---------------------------------------------------------
# The centroid pipeline is fully deterministic by construction:
#   - No random sampling of control genes (cf. sc.tl.score_genes' ctrl_size).
#   - Subpops are iterated in sorted() order.
#   - Gene order is taken from the parquet schema (deterministic on disk).
#   - Mean, matmul, and L2-norm are all reduction operations whose output
#     depends only on inputs (modulo last-bit BLAS thread-order effects,
#     which do not change top1 except for cells whose margin is at the
#     fp32 noise floor -- those cells are exactly what the downstream
#     margin QC is designed to flag).
# RANDOM_SEED is set defensively even though no numpy RNG is used, so that
# any future addition that does use one inherits a fixed seed automatically.
# The same RANDOM_SEED is used for every broad cell type (every SLURM job),
# making the per-broad runs identical in their (non-)randomness.
RANDOM_SEED = 0


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--broad_cell_type', required=True, choices=list(BROAD_TO_PREFIX))
    p.add_argument('--matrix_dir',    default=DEFAULT_MATRIX_DIR)
    p.add_argument('--meta_path',     default=DEFAULT_META_PATH)
    p.add_argument('--mathys_matrix', default=DEFAULT_MATHYS_MATRIX)
    p.add_argument('--mathys_meta',   default=DEFAULT_MATHYS_META)
    p.add_argument('--output_dir',    default=DEFAULT_OUTPUT_DIR)
    p.add_argument('--force', action='store_true',
                   help='Overwrite output if it already exists.')
    return p.parse_args()


def main():
    args = parse_args()
    broad  = args.broad_cell_type
    prefix = BROAD_TO_PREFIX[broad]

    # Defensive seed: no RNG is invoked by this script, but if future edits
    # ever add one, this guarantees the same per-broad behavior.
    np.random.seed(RANDOM_SEED)

    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, f'seaad_subpop_centroid_{broad}.parquet')
    if os.path.exists(out_path) and not args.force:
        print(f'[skip] {out_path} already exists; pass --force to overwrite.')
        return 0

    t0 = time.time()

    # ---------- [1/6] Load Mathys metadata & matrix, build per-subpop centroids ----------
    print(f'[1/6] Loading Mathys metadata {args.mathys_meta}', flush=True)
    mref_meta = pd.read_parquet(args.mathys_meta)
    for col in ('TAG', 'Subcluster'):
        if col not in mref_meta.columns:
            sys.exit(f'Mathys metadata missing required column: {col!r}')
    mref_meta = mref_meta[['TAG', 'Subcluster']].dropna(subset=['Subcluster']).copy()
    mref_meta['TAG'] = mref_meta['TAG'].astype(str)
    mref_meta['Subcluster'] = mref_meta['Subcluster'].astype(str)

    sub_rx = rf'^{prefix}\d+$'
    mref_meta_b = mref_meta.loc[mref_meta['Subcluster'].str.match(sub_rx)]
    subpops = sorted(mref_meta_b['Subcluster'].unique())
    if not subpops:
        sys.exit(f'No Mathys subpops found for broad={broad} (prefix={prefix})')
    counts = mref_meta_b['Subcluster'].value_counts().reindex(subpops)
    print(f'        {len(subpops)} broad-matched subpops; Mathys cells per subpop:', flush=True)
    for sp in subpops:
        print(f'          {sp}: {int(counts[sp]):,}', flush=True)

    print(f'        loading Mathys matrix {args.mathys_matrix}', flush=True)
    mref = pd.read_parquet(args.mathys_matrix)  # genes x cells, float64, log-normalized
    # The Mathys matrix is genes (rows) x cells (cols); guard against an unexpected orientation.
    if mref.shape[0] > mref.shape[1]:
        sys.exit(f'Mathys matrix shape={mref.shape} looks transposed (expected genes < cells)')
    mathys_genes = np.asarray(mref.index.astype(str))
    mathys_cells_set = set(mref.columns.astype(str))
    miss = [c for c in mref_meta_b['TAG'].tolist() if c not in mathys_cells_set]
    if miss:
        sys.exit(f'{len(miss)} Mathys cells in metadata are missing from matrix columns '
                 f'(first 3: {miss[:3]})')
    # Cast Mathys values to float32 once; release the float64 DataFrame.
    mref_values = mref.to_numpy(dtype=np.float32)  # (n_genes_m, n_cells_m)
    mref_col_idx = {c: i for i, c in enumerate(mref.columns.astype(str))}
    del mref
    gc.collect()

    n_genes_m = mref_values.shape[0]
    C_mref = np.empty((len(subpops), n_genes_m), dtype=np.float32)
    for i, sp in enumerate(subpops):
        sp_tags = mref_meta_b.loc[mref_meta_b['Subcluster'] == sp, 'TAG'].tolist()
        col_idxs = np.array([mref_col_idx[c] for c in sp_tags], dtype=np.int64)
        if col_idxs.size == 0:
            sys.exit(f'No Mathys cells found for subpop {sp}')
        # Mean across this subpop's cells -> one centroid vector of length n_genes_m
        C_mref[i] = mref_values[:, col_idxs].mean(axis=1)
    del mref_values
    gc.collect()
    print(f'        built {len(subpops)} centroids of length {n_genes_m}', flush=True)

    # ---------- [2/6] Load SEA-AD matrix (streaming, identical to label_seaad_subpops.py) ----------
    matrix_path = os.path.join(args.matrix_dir, f'SEAAD_Matrix_{broad}.parquet')
    print(f'[2/6] Loading SEA-AD matrix {matrix_path} (streaming via pyarrow)', flush=True)
    pf = pq.ParquetFile(matrix_path)
    field_names = pf.schema_arrow.names
    if 'TAG' not in field_names:
        sys.exit(f'Expected "TAG" field in parquet schema; got first fields {field_names[:3]}')
    tag_field_idx = field_names.index('TAG')
    gene_cols = [n for n in field_names if n != 'TAG']
    n_rows = pf.metadata.num_rows
    X = np.empty((n_rows, len(gene_cols)), dtype=np.float32)
    tags = np.empty(n_rows, dtype=object)
    row = 0
    for batch in pf.iter_batches(batch_size=50_000):
        n = batch.num_rows
        tags[row:row + n] = batch.column(tag_field_idx).to_numpy(zero_copy_only=False)
        j = 0
        for k in range(batch.num_columns):
            if k == tag_field_idx:
                continue
            X[row:row + n, j] = batch.column(k).to_numpy(zero_copy_only=False)
            j += 1
        row += n
        del batch
    if row != n_rows:
        sys.exit(f'Streaming load row mismatch: expected {n_rows}, got {row}')
    mat = pd.DataFrame(X, index=pd.Index(tags, name='TAG'), columns=gene_cols, copy=False)
    if mat.index.name != 'TAG':
        sys.exit(f'Expected matrix index named "TAG"; got {mat.index.name!r}')
    print(f'        matrix shape={mat.shape}  dtype={mat.dtypes.iloc[0]}', flush=True)

    # ---------- [3/6] Load SEA-AD metadata (identical to label_seaad_subpops.py) ----------
    print(f'[3/6] Loading metadata {args.meta_path}', flush=True)
    meta = pd.read_parquet(args.meta_path)
    if meta.index.name != 'TAG':
        sys.exit(f'Expected metadata index named "TAG"; got {meta.index.name!r}')
    meta_broad = meta.loc[meta['broad_cell_type'] == broad]
    missing = mat.index.difference(meta_broad.index)
    if len(missing):
        sys.exit(f'{len(missing)} matrix TAGs not found in metadata for broad={broad} '
                 f'(first 3: {list(missing[:3])})')
    obs = meta_broad.reindex(mat.index).copy()
    if not obs.index.equals(mat.index):
        sys.exit('TAG alignment failed between matrix and metadata after reindex.')

    # ---------- [4/6] Intersect genes (SEA-AD x Mathys) and subset both matrices ----------
    print(f'[4/6] Intersecting SEA-AD x Mathys gene panels', flush=True)
    mathys_gene_set = set(mathys_genes.tolist())
    sa_gene_arr = np.asarray(gene_cols)
    keep_mask = np.array([g in mathys_gene_set for g in sa_gene_arr], dtype=bool)
    shared_genes = sa_gene_arr[keep_mask]
    if shared_genes.size == 0:
        sys.exit('No shared genes between SEA-AD and Mathys gene panels.')
    print(f'        shared: {shared_genes.size:,} / {sa_gene_arr.size:,} SEA-AD / '
          f'{n_genes_m:,} Mathys', flush=True)

    # Reorder Mathys centroids to match the SEA-AD shared-gene order.
    mathys_pos = {g: i for i, g in enumerate(mathys_genes.tolist())}
    m_indices = np.array([mathys_pos[g] for g in shared_genes], dtype=np.int64)
    C_shared = np.ascontiguousarray(C_mref[:, m_indices], dtype=np.float32)  # (n_subpops, n_shared)
    del C_mref
    gc.collect()

    # Subset SEA-AD X to shared genes. This is the peak-memory step: briefly
    # holds X (n_cells x n_genes_all) and X_shared (n_cells x n_shared) together.
    sa_indices = np.where(keep_mask)[0]
    X_shared = np.ascontiguousarray(X[:, sa_indices], dtype=np.float32)
    del mat, X
    gc.collect()
    print(f'        X_shared: {X_shared.shape}  C_shared: {C_shared.shape}', flush=True)

    # ---------- [5/6] Pearson correlation: each SEA-AD cell vs each Mathys centroid ----------
    print(f'[5/6] Computing Pearson centroid correlations', flush=True)
    # Mean-center along the gene axis (axis=1) for both matrices.
    X_shared -= X_shared.mean(axis=1, keepdims=True)
    C_shared = C_shared - C_shared.mean(axis=1, keepdims=True)
    X_norm = np.linalg.norm(X_shared, axis=1)  # (n_cells,)
    C_norm = np.linalg.norm(C_shared, axis=1)  # (n_subpops,)
    # Guard against zero-norm rows (constant expression -> correlation undefined).
    n_zero_cells = int((X_norm == 0).sum())
    n_zero_subpops = int((C_norm == 0).sum())
    X_norm_safe = np.where(X_norm == 0, 1.0, X_norm).astype(np.float32)
    C_norm_safe = np.where(C_norm == 0, 1.0, C_norm).astype(np.float32)
    # (n_cells, n_shared) @ (n_shared, n_subpops) -> (n_cells, n_subpops)
    scores = (X_shared @ C_shared.T) / (X_norm_safe[:, None] * C_norm_safe[None, :])
    scores = scores.astype(np.float32)
    if n_zero_cells:
        print(f'        [warn] {n_zero_cells} SEA-AD cells had zero-norm expression -> NaN row',
              flush=True)
        scores[X_norm == 0, :] = np.nan
    if n_zero_subpops:
        bad = [sp for sp, n in zip(subpops, C_norm) if n == 0]
        print(f'        [warn] zero-norm Mathys centroids (constant across genes): {bad}',
              flush=True)
        for j, sp in enumerate(subpops):
            if C_norm[j] == 0:
                scores[:, j] = np.nan
    del X_shared, C_shared
    gc.collect()

    # ---------- [6/6] Top1/top2/margin and write parquet (identical to label_seaad_subpops.py) ----------
    print(f'[6/6] Computing top1/top2/margin and writing parquet', flush=True)
    # If any subpop column is all-NaN (zero-norm centroid), drop from argmax search.
    col_all_nan = np.isnan(scores).all(axis=0)
    if col_all_nan.any():
        dropped = [sp for sp, bad in zip(subpops, col_all_nan) if bad]
        print(f'        dropping all-NaN score columns from argmax: {dropped}', flush=True)
    valid_subpops_list = [sp for sp, bad in zip(subpops, col_all_nan) if not bad]
    valid_subpops = np.array(valid_subpops_list)
    s = scores[:, ~col_all_nan]

    n_cells, n_valid = s.shape
    if n_valid == 0:
        sys.exit('All score columns were NaN; nothing to rank.')

    if n_valid == 1:
        top1_idx = np.zeros(n_cells, dtype=np.int64)
        top2_idx = np.zeros(n_cells, dtype=np.int64)
        top2_valid = False
    else:
        part = np.argpartition(s, kth=[n_valid - 2, n_valid - 1], axis=1)
        cand = part[:, -2:]
        rows = np.arange(n_cells)
        cand_scores = s[rows[:, None], cand]
        order_within = np.argsort(-cand_scores, axis=1, kind='stable')
        top1_idx = cand[rows, order_within[:, 0]]
        top2_idx = cand[rows, order_within[:, 1]]
        top2_valid = True

    rows = np.arange(n_cells)
    top1_score = s[rows, top1_idx]
    if top2_valid:
        top2_score = s[rows, top2_idx]
    else:
        top2_score = np.full(n_cells, np.nan, dtype=np.float32)
    top1 = valid_subpops[top1_idx]
    top2 = (valid_subpops[top2_idx] if top2_valid
            else np.array([''] * n_cells, dtype=object))
    margin = top1_score - top2_score

    out = pd.DataFrame(index=obs.index.copy())
    out.index.name = 'TAG'
    for c in ['Donor ID', 'Subclass', 'broad_cell_type', 'alzheimers_or_control']:
        if c in obs.columns:
            out[c] = obs[c].to_numpy()
    for i, sp in enumerate(subpops):
        out[f'score_{sp}'] = scores[:, i]
    out['top1']       = top1
    out['top1_score'] = top1_score
    out['top2']       = top2
    out['top2_score'] = top2_score
    out['margin']     = margin

    out.to_parquet(out_path)
    elapsed = time.time() - t0
    print(f'[done] wrote {out_path} '
          f'({out.shape[0]:,} cells x {out.shape[1]} cols) in {elapsed/60:.1f} min',
          flush=True)

    print('\nTop1 assignment counts:')
    print(out['top1'].value_counts().to_string())
    print('\nMargin summary:')
    print(out['margin'].describe().to_string())
    print('\nSubclass x top1 crosstab (first 20 rows):')
    ct = pd.crosstab(out['Subclass'], out['top1'])
    print(ct.iloc[:20].to_string())

    return 0


if __name__ == '__main__':
    sys.exit(main())
