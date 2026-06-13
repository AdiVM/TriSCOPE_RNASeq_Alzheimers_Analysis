"""
label_seaad_subpops.py

Stage A of the SEA-AD -> Mathys subpopulation labeling pipeline.

For one SEA-AD broad cell type (Ex/Inh/Ast/Oli/Mic/Opc), scores every cell
against the Mathys subpopulations matching that broad type using
sc.tl.score_genes, then records the argmax (top1), second-best (top2),
and margin (top1 - top2) per cell -- alongside the full per-subpop score
matrix.

Confidence filtering on the margin is deliberately NOT done here; it is a
downstream decision made after inspecting the score distribution in the
Stage A QC notebook.

================================================================================
RESULTS LOG -- 2026-05-20 first end-to-end run on all 6 broads
================================================================================
This sc.tl.score_genes-based script produced defensible labels for 5/6 broads
but was degenerate for Ex. See label_seaad_subpops_centroid.py for the fix.

Per-broad outcomes (job 40818944-40818948, smoke test 40817460):

  broad | n_cells | elapsed | MaxRSS  | ReqMem | top1 distribution
  ------+---------+---------+---------+--------+-----------------------------
  Opc   |  26,654 |   1:13  |  10.4 G |  116G  | OK: 3 subpops, healthy spread
  Mic   |  41,187 |   1:10  |  20.9 G |  124G  | OK: 4 subpops, healthy spread
  Ast   |  83,983 |   1:59  |  29.0 G |  140G  | OK: 4 subpops, healthy spread
  Oli   | 142,064 |   3:29  |  51.8 G |  160G  | OK: 5 subpops (Oli1 only 3 markers)
  Inh   | 393,273 |  23:29  | 141.7 G |  160G  | OK: 12 subpops, clean Subclass map
  Ex    | 610,593 |  34:00  | 219.6 G |  249G  | DEGENERATE: 99.2% -> Ex2

Ex degeneracy root cause: Mathys Ex2's 52 markers (GRIA4, GRIN1, SYN3, KCNIP4,
CSMD1, EPHA6, CBLN2, CNTNAP5, LRRTM4, TENM2, FLRT2, ...) are synaptic /
cell-adhesion genes that are uniquely-in-Mathys-Ex2 within the Mathys dataset
but pan-expressed across all SEA-AD cortical Ex cells. Marker-name-only scoring
(score_genes) cannot distinguish them because the gene set is essentially a
generic glutamatergic-synaptic signature once you leave the Mathys reference
context. Overlap with other Mathys Ex subpops was low (max 25% with Ex0, mostly
<10%), so this is NOT a marker-leakage problem -- it's that Mathys defined Ex2
as the "vanilla glutamatergic" cluster *relative to its more specialized
neighbors*, and that contrast does not transfer to SEA-AD.

Resolution: rerun the labeling using Pearson centroid correlation against the
Mathys reference matrix instead of marker-set scoring. That is implemented in
label_seaad_subpops_centroid.py, which keeps everything in this script
identical (loader, metadata join, top1/top2/margin computation, parquet output
schema, SLURM submission pattern) and only swaps the inner scoring loop.

Per-broad SLURM budgets observed here remain usable as guidance for any future
score_genes-style rerun on this dataset:
  short partition, mem >= 1.5 * MaxRSS, time as listed above.
================================================================================

Usage:
    python label_seaad_subpops.py --broad_cell_type Ex

Inputs (defaults match the project layout):
  --markers    Mathys subpop marker table (Excel, sheet "Subpopulation_markers", header=1)
  --matrix_dir Directory containing SEAAD_Matrix_{broad}.parquet (cells x genes, log2(CPM+1))
  --meta_path  SEAAD_CellMetadata.parquet (TAG index, has broad_cell_type and Subclass)

Output:
  {output_dir}/seaad_subpop_scores_{broad}.parquet
    Columns: per-cell metadata (Donor ID, Subclass, broad_cell_type, alzheimers_or_control)
             + score_<subpop> for every broad-matched Mathys subpop
             + top1, top1_score, top2, top2_score, margin
"""

import argparse
import gc
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scanpy as sc


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

DEFAULT_MATRIX_DIR = '/n/groups/patel/adithya/SEAAD_Outputs'
DEFAULT_META_PATH  = '/n/groups/patel/adithya/SEAAD_Outputs/SEAAD_CellMetadata.parquet'
DEFAULT_MARKERS    = '/n/groups/patel/adithya/Alz_Outputs/Final_Outputs/seaad_pseudotime/mathys_markers.xlsx'
DEFAULT_OUTPUT_DIR = '/n/groups/patel/adithya/Alz_Outputs/Final_Outputs/seaad_pseudotime/labels'


def _subpops_for_prefix(markers_df, prefix):
    rx = rf'^{prefix}\d+$'
    return sorted(markers_df.loc[markers_df['subpopulation'].astype(str).str.match(rx),
                                 'subpopulation'].unique())


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--broad_cell_type', required=True, choices=list(BROAD_TO_PREFIX))
    p.add_argument('--matrix_dir', default=DEFAULT_MATRIX_DIR)
    p.add_argument('--meta_path',  default=DEFAULT_META_PATH)
    p.add_argument('--markers',    default=DEFAULT_MARKERS)
    p.add_argument('--output_dir', default=DEFAULT_OUTPUT_DIR)
    p.add_argument('--ctrl_size',    type=int, default=50)
    p.add_argument('--n_bins',       type=int, default=25)
    p.add_argument('--random_state', type=int, default=0)
    p.add_argument('--force', action='store_true',
                   help='Overwrite output if it already exists.')
    return p.parse_args()


def main():
    args = parse_args()
    broad  = args.broad_cell_type
    prefix = BROAD_TO_PREFIX[broad]

    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, f'seaad_subpop_scores_{broad}.parquet')
    if os.path.exists(out_path) and not args.force:
        print(f'[skip] {out_path} already exists; pass --force to overwrite.')
        return 0

    t0 = time.time()

    print(f'[1/5] Loading markers from {args.markers}', flush=True)
    mk = pd.read_excel(args.markers, sheet_name='Subpopulation_markers', header=1)
    mk = mk[['gene.name', 'subpopulation']].dropna(subset=['subpopulation', 'gene.name'])
    subpops = _subpops_for_prefix(mk, prefix)
    if not subpops:
        sys.exit(f'No Mathys subpops found for broad={broad} (prefix={prefix})')
    markers = {sp: mk.loc[mk['subpopulation'] == sp, 'gene.name'].astype(str).tolist()
               for sp in subpops}
    print(f'        {len(subpops)} broad-matched subpops: {subpops}', flush=True)

    matrix_path = os.path.join(args.matrix_dir, f'SEAAD_Matrix_{broad}.parquet')
    print(f'[2/5] Loading SEA-AD matrix {matrix_path} (streaming via pyarrow)', flush=True)
    # Stream parquet row-batches into a pre-allocated float32 array to avoid the
    # ~2x peak that pd.read_parquet incurs (arrow buffers + pandas frame + to_numpy copy).
    # Peak here is X (n_cells x n_genes x 4B) plus one batch's worth (~7GB at batch_size=50k for Ex).
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
    # Wrap as DataFrame without copying X so downstream code is unchanged.
    mat = pd.DataFrame(X, index=pd.Index(tags, name='TAG'), columns=gene_cols, copy=False)
    if mat.index.name != 'TAG':
        sys.exit(f'Expected matrix index named "TAG"; got {mat.index.name!r}')
    print(f'        matrix shape={mat.shape}  dtype={mat.dtypes.iloc[0]}', flush=True)

    print(f'[3/5] Loading metadata {args.meta_path}', flush=True)
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

    print(f'[4/5] Building AnnData and scoring {len(subpops)} subpops', flush=True)
    adata = ad.AnnData(
        X=mat.to_numpy(),
        obs=obs,
        var=pd.DataFrame(index=mat.columns.astype(str)),
    )
    del mat
    gc.collect()
    print(f'        AnnData: n_obs={adata.n_obs:,}  n_vars={adata.n_vars:,}', flush=True)

    var_set = set(adata.var_names)
    score_cols = []
    for sp in subpops:
        gl_full = markers[sp]
        gl = [g for g in gl_full if g in var_set]
        col = f'score_{sp}'
        if len(gl) == 0:
            print(f'        [warn] {sp}: 0/{len(gl_full)} markers in SEA-AD panel -> NaN score',
                  flush=True)
            adata.obs[col] = np.nan
        else:
            sc.tl.score_genes(
                adata,
                gene_list=gl,
                ctrl_size=args.ctrl_size,
                n_bins=args.n_bins,
                score_name=col,
                random_state=args.random_state,
                use_raw=False,
                copy=False,
            )
            print(f'        scored {sp}: used {len(gl)}/{len(gl_full)} markers',
                  flush=True)
        score_cols.append(col)

    print(f'[5/5] Computing top1/top2/margin and writing parquet', flush=True)
    scores = adata.obs[score_cols].to_numpy(dtype=np.float32)  # (cells, n_subpops)

    # If any score column is all-NaN (no markers overlapped), drop it from the
    # argmax search so it doesn't poison every cell's top1.
    col_all_nan = np.isnan(scores).all(axis=0)
    if col_all_nan.any():
        dropped = [c for c, bad in zip(score_cols, col_all_nan) if bad]
        print(f'        dropping all-NaN score columns from argmax: {dropped}',
              flush=True)
    valid_cols = [c for c, bad in zip(score_cols, col_all_nan) if not bad]
    valid_subpops = np.array([c[len('score_'):] for c in valid_cols])
    s = scores[:, ~col_all_nan]

    n_cells, n_valid = s.shape
    if n_valid == 0:
        sys.exit('All score columns were NaN; nothing to rank.')

    # Argpartition for top-2 (faster than full sort)
    if n_valid == 1:
        top1_idx = np.zeros(n_cells, dtype=np.int64)
        top2_idx = np.zeros(n_cells, dtype=np.int64)
        top2_valid = False
    else:
        # argpartition with kth=[-1, -2] puts top-2 at the last two positions (unordered)
        part = np.argpartition(s, kth=[n_valid - 2, n_valid - 1], axis=1)
        cand = part[:, -2:]   # two largest per row, order within pair not guaranteed
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

    out = pd.DataFrame(index=adata.obs_names.copy())
    out.index.name = 'TAG'
    for c in ['Donor ID', 'Subclass', 'broad_cell_type', 'alzheimers_or_control']:
        if c in adata.obs.columns:
            out[c] = adata.obs[c].to_numpy()
    for c in score_cols:
        out[c] = adata.obs[c].to_numpy()
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
