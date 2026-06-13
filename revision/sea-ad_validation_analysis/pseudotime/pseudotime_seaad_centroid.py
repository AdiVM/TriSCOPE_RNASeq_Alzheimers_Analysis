"""
pseudotime_seaad_centroid.py

Per-subcluster pseudotime pipeline on SEA-AD, using centroid-labeled
subpopulations (top1 from label_seaad_subpops_centroid.py outputs).

Script body is a copy of the Mathys per-subcluster pseudotime analysis from
/home/adm808/Analysis_Scripts/Plotting_Figures/pseudotime_hvg.ipynb (cell 0)
with the following and ONLY the following changes:

  1. Data loaded from SEA-AD parquets instead of the Mathys CellMatrix:
       expression  -> /n/groups/patel/adithya/SEAAD_Outputs/SEAAD_Matrix_{broad}.parquet
                      (cells x genes, float32, already CPM->log2(x+1) normalized)
       metadata    -> /n/groups/patel/adithya/SEAAD_Outputs/SEAAD_CellMetadata.parquet
       subcluster  -> top1 from .../labels_centroid/seaad_subpop_centroid_{broad}.parquet
                      (Mathys subpop name assigned per SEA-AD cell by Pearson centroid corr)

  2. Two scanpy normalization lines REMOVED because the SEA-AD matrix is
     already log-normalized (see sea_ad_normalization_matrix.py):
       sc.pp.normalize_total(adata_template, target_sum=1e4)
       sc.pp.log1p(adata_template)

  3. Mathys clinical column names replaced with the SEA-AD column names that
     carry the same information (no derived/injected columns):
       dcfdx_lv in {4,5}  ->  alzheimers_or_control == 1
       ceradsc            ->  'Overall AD neuropathological Change'
                              (ordered Categorical: High < Intermediate < Low < Not AD,
                               so .max() = 'Not AD' = healthiest, matching ceradsc=4)
       cogdx              ->  'Cognitive Status'
                              (ordered Categorical: No dementia < Dementia,
                               so .min() = 'No dementia' = healthiest, matching cogdx=1)
       projid             ->  'Donor ID'
       Subcluster         ->  centroid top1

     The CERAD classifier binary outcome 'ceradsc in {1,2}' (definite/probable AD)
     becomes 'Overall AD neuropathological Change' in {'High','Intermediate'}
     (high/moderate AD pathology) -- same clinical definition.

     The pd.to_numeric(...) coercions that the Mathys script applied to ceradsc
     and cogdx were no-ops on already-typed numeric columns; with SEA-AD's
     ordered Categorical encoding those coercions are removed (they would error
     on string categories). The .max()/.min()/.isin() and dropna() that consume
     those columns are otherwise unchanged.

  4. The outer for-loop runs for the single subcluster passed on the CLI; the
     per-subcluster validation row and ROC interpolation arrays are written
     to their own files so a separate aggregator can combine them into the
     final publication CSV + ROC plots after all jobs complete.

Usage:
    python pseudotime_seaad_centroid.py --subcluster Ex5
"""

import argparse
import os
import re
import shutil

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, roc_curve
from scipy.stats import combine_pvalues, spearmanr
from statsmodels.stats.multitest import multipletests
from statsmodels.nonparametric.smoothers_lowess import lowess
import seaborn as sns


# Mathys subcluster prefix -> SEA-AD broad_cell_type (matches the centroid pipeline).
SUBCLUSTER_PREFIX_TO_BROAD = {
    'Ex': 'Ex', 'In': 'Inh', 'Ast': 'Ast',
    'Oli': 'Oli', 'Mic': 'Mic', 'Opc': 'Opc',
}

# Ordered Categoricals so .max()/.min() pick the healthiest level on each axis:
#   Mathys ceradsc: 1 (def AD) < 2 < 3 < 4 (no AD)        -> max = healthiest
#   SEA-AD pathology: High < Intermediate < Low < Not AD  -> max = 'Not AD'    (healthiest)
#   Mathys cogdx: 1 (NCI) < 2 < ... (higher = more dementia) -> min = healthiest
#   SEA-AD cognition: No dementia < Dementia              -> min = 'No dementia' (healthiest)
PATHOLOGY_CATEGORIES = ['High', 'Intermediate', 'Low', 'Not AD']
COGNITION_CATEGORIES = ['No dementia', 'Dementia']


def subcluster_to_broad(sub):
    m = re.match(r'^([A-Za-z]+)\d+$', sub)
    if not m:
        raise ValueError(f'Unrecognized subcluster name: {sub!r}')
    prefix = m.group(1)
    if prefix not in SUBCLUSTER_PREFIX_TO_BROAD:
        raise ValueError(f'Unknown subcluster prefix {prefix!r} for {sub!r}')
    return SUBCLUSTER_PREFIX_TO_BROAD[prefix]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--subcluster', required=True,
                   help='Mathys subpop name (e.g. Ex5, In3, Opc2) to process.')
    p.add_argument('--matrix_dir',
                   default='/n/groups/patel/adithya/SEAAD_Outputs')
    p.add_argument('--meta_path',
                   default='/n/groups/patel/adithya/SEAAD_Outputs/SEAAD_CellMetadata.parquet')
    p.add_argument('--labels_dir',
                   default='/n/groups/patel/adithya/Alz_Outputs/Final_Outputs/seaad_pseudotime/labels_centroid')
    p.add_argument('--output_dir',
                   default='/n/groups/patel/adithya/Alz_Outputs/Final_Outputs/seaad_pseudotime/pseudotime_centroid')
    return p.parse_args()


def main():
    args = parse_args()
    broad = subcluster_to_broad(args.subcluster)

    # --- Set up Paths ---
    expr_path   = os.path.join(args.matrix_dir, f'SEAAD_Matrix_{broad}.parquet')
    meta_path   = args.meta_path
    labels_path = os.path.join(args.labels_dir, f'seaad_subpop_centroid_{broad}.parquet')
    output_path = args.output_dir + '/'
    os.makedirs(output_path, exist_ok=True)

    # --- Load Data ---
    # SEA-AD per-cell metadata (TAG index).
    meta_df = pd.read_parquet(meta_path)

    # Centroid labels for this broad (TAG index, has 'top1' column with assigned Mathys subpop).
    labels_df = pd.read_parquet(labels_path)
    # 'Subcluster' column = centroid-assigned Mathys subpop (replaces Mathys's own Subcluster).
    meta_df = meta_df.join(labels_df['top1'].rename('Subcluster'), how='inner')

    # SEA-AD pathology and cognition columns -> ordered Categoricals so .max()/.min()
    # preserve clinical semantics directly (no derived numeric columns).
    meta_df['Overall AD neuropathological Change'] = pd.Categorical(
        meta_df['Overall AD neuropathological Change'],
        categories=PATHOLOGY_CATEGORIES, ordered=True,
    )
    meta_df['Cognitive Status'] = pd.Categorical(
        meta_df['Cognitive Status'],
        categories=COGNITION_CATEGORIES, ordered=True,
    )

    # AD/Control label (direct replacement for Mathys's `dcfdx_lv in {4,5}` derivation).
    meta_df['AD'] = meta_df['alzheimers_or_control'].map({1: 'AD', 0: 'Control'}).astype('category')
    print("AD status defined based on 'alzheimers_or_control':")
    print(meta_df['AD'].value_counts())

    # Mathys metadata had TAG as a column; restore that so the script body's
    # meta_sub['TAG'].tolist() / meta_sub.set_index('TAG') work unchanged.
    meta_df = meta_df.reset_index()

    # SEA-AD expression matrix (cells on rows; no transpose needed -- the Mathys
    # script transposed because its source was genes-on-rows).
    expr_df = pd.read_parquet(expr_path)
    expr_df_transposed = expr_df  # name kept identical to the Mathys script body

    # --- Analysis Settings ---
    all_subclusters = [args.subcluster]  # one per SLURM job

    random_seeds = [41, 42, 43]
    # 57, 62, 72, 73

    validation_results = []
    summary_roc_data = []
    summary_roc_data_cerad = []
    roc_plot_data_cerad = []

    # --- Main Loop (single iteration; loop kept for byte-identical script body) ---
    for current_subcluster in all_subclusters:
        print(f"\n--- Processing Subcluster: {current_subcluster} ---")

        meta_sub = meta_df[meta_df['Subcluster'] == current_subcluster]

        if meta_sub.shape[0] < 50 or len(meta_sub['AD'].unique()) < 2:
            print(f"Skipping {current_subcluster}: <50 cells or only one AD/Control group.")
            continue

        safe_subcluster_name = re.sub(r'[^a-zA-Z0-9_.-]', '_', current_subcluster)
        subcluster_cells = meta_sub['TAG'].tolist()
        expr_sub = expr_df_transposed.loc[subcluster_cells]

        # --- Create AnnData ---
        adata_template = sc.AnnData(X=expr_sub)
        adata_template.obs = meta_sub.set_index('TAG').loc[adata_template.obs.index]

        # --- NON-CIRCULAR: Use HVGs instead of DEGs ---
        # Two normalization lines from the Mathys script are REMOVED here because
        # the SEA-AD matrix is already CPM -> log2(x+1) normalized upstream
        # (see sea_ad_normalization_matrix.py):
        #     sc.pp.normalize_total(adata_template, target_sum=1e4)
        #     sc.pp.log1p(adata_template)
        print("Selecting HVGs for pseudotime (non-circular approach)...")

        # Select top 2000 HVGs (or all genes if fewer available)
        n_top = min(2000, adata_template.n_vars - 1)
        sc.pp.highly_variable_genes(adata_template, n_top_genes=n_top, flavor='seurat')
        hvg_genes = adata_template.var[adata_template.var['highly_variable']].index.tolist()

        if len(hvg_genes) < 50:
            print(f"Skipping {current_subcluster}: only {len(hvg_genes)} HVGs found.")
            continue

        adata_template = adata_template[:, hvg_genes].copy()
        print(f"Using {len(hvg_genes)} HVGs to build trajectories.")

        # --- PCA (no Harmony - not needed within subclusters) ---
        n_pcs = min(30, len(hvg_genes) - 1)
        sc.pp.pca(adata_template, n_comps=n_pcs)
        print(f"PCA complete: {adata_template.obsm['X_pca'].shape}")

        # --- FAST neighbor graph ---
        print(f"Building neighbor graph for {adata_template.n_obs} cells...")
        sc.pp.neighbors(adata_template, n_neighbors=10, n_pcs=n_pcs, method='umap')
        print("Neighbor graph complete!")

        # --- Downstream: UMAP, clustering, diffusion map ---
        sc.tl.umap(adata_template)
        sc.tl.leiden(adata_template, resolution=0.5)
        sc.tl.diffmap(adata_template)

        # --- Loop over random seeds to test trajectory robustness ---
        subcluster_aucs = []
        roc_plot_data = []
        top_trajectory_genes_from_gt = []  # will be filled after seed loop
        top_50_trajectory_genes_from_gt = []  # defensive init (used in final row regardless)

        subcluster_pvals, subcluster_odds_ratios = [], []
        all_fprs, all_tprs = [], []  # For mean ROC over seeds

        cerad_aucs = []
        cerad_odds_ratios = []
        cerad_pvals = []
        cerad_or_cis = []

        # Store pseudotime per seed (for later averaging)
        pseudotime_per_seed = []

        for seed in random_seeds:
            print(f"\n  -- Running analysis for seed: {seed} --")

            # Create a fresh copy for each run
            adata = adata_template.copy()

            # --- Select Root Cell Based on pathology and cognition (with current seed) ---
            # SEA-AD analogs of Mathys (max ceradsc, then min cogdx):
            #   max('Overall AD neuropathological Change') = 'Not AD'  (healthiest pathology)
            #   min('Cognitive Status')                    = 'No dementia' (healthiest cognition)
            subcluster_meta_dpt = adata.obs.copy()
            subcluster_meta_dpt.dropna(
                subset=['Overall AD neuropathological Change', 'Cognitive Status'],
                inplace=True,
            )

            if subcluster_meta_dpt.empty:
                print("  Skipping trajectory ordering: No cells with valid pathology/cognition values.")
                break

            max_pathology = subcluster_meta_dpt['Overall AD neuropathological Change'].max()
            root_candidates = subcluster_meta_dpt[
                subcluster_meta_dpt['Overall AD neuropathological Change'] == max_pathology
            ]
            min_cognition = root_candidates['Cognitive Status'].min()
            root_candidates = root_candidates[
                root_candidates['Cognitive Status'] == min_cognition
            ]

            root_cell_id = root_candidates.sample(n=1, random_state=seed).index[0]
            adata.uns['iroot'] = np.where(adata.obs_names == root_cell_id)[0][0]

            # --- Calculate Pseudotime (DPT) ---
            sc.tl.dpt(adata)

            # Store pseudotime for this seed (full vector, all cells)
            pseudotime_per_seed.append(adata.obs['dpt_pseudotime'].to_numpy())

            # --- Calculate AUC for the current seed ---
            validation_df = adata.obs[['dpt_pseudotime', 'AD']].copy()
            validation_df.dropna(subset=['dpt_pseudotime'], inplace=True)
            y_true = (validation_df['AD'] == 'AD').astype(int)
            y_score = validation_df['dpt_pseudotime']

            auc_score = roc_auc_score(y_true=y_true, y_score=y_score)
            subcluster_aucs.append(auc_score)
            print(f"  --> Seed {seed}: AUC = {auc_score:.4f}")

            # Logistic regression for this seed
            X = sm.add_constant(validation_df['dpt_pseudotime'])
            if len(y_true.unique()) < 2:
                print(f"  Skipping seed {seed}: only one class present after dropping NaNs.")
                continue

            log_reg = sm.Logit(y_true, X).fit(disp=0)
            subcluster_pvals.append(log_reg.pvalues['dpt_pseudotime'])
            subcluster_odds_ratios.append(np.exp(log_reg.params['dpt_pseudotime']))

            # ROC curve for this seed (clinical AD)
            fpr, tpr, _ = roc_curve(y_true=y_true, y_score=y_score)
            all_fprs.append(fpr)
            all_tprs.append(tpr)

            if auc_score > 0.55:
                roc_plot_data.append({'fpr': fpr, 'tpr': tpr, 'auc': auc_score, 'seed': seed})

            # CERAD-based classifier
            # SEA-AD analog of Mathys' ceradsc.isin([1,2]) (definite/probable AD):
            #     'Overall AD neuropathological Change'.isin(['High','Intermediate'])
            try:
                cerad_df = adata.obs[['dpt_pseudotime']].copy()
                # Clear the index name so the next line doesn't create a collision
                # between an index named 'TAG' and a column named 'TAG' (which newer
                # pandas rejects on .merge as ambiguous; latent bug from Mathys notebook).
                cerad_df.index.name = None
                cerad_df['TAG'] = cerad_df.index

                meta_subset = meta_df[['TAG', 'Overall AD neuropathological Change']].dropna()

                cerad_df = cerad_df.merge(meta_subset, on="TAG", how="left")
                cerad_df = cerad_df.dropna(subset=["Overall AD neuropathological Change", "dpt_pseudotime"])
                cerad_df['AD_cerad'] = cerad_df['Overall AD neuropathological Change'].isin(['High', 'Intermediate']).astype(int)

                if cerad_df['AD_cerad'].nunique() < 2:
                    print(f"  CERAD classifier skipped for seed {seed}: only one class present.")
                else:
                    X_cerad = sm.add_constant(cerad_df['dpt_pseudotime'])
                    y_cerad = cerad_df['AD_cerad']
                    cerad_model = sm.Logit(y_cerad, X_cerad).fit(disp=0)

                    cerad_auc = roc_auc_score(y_cerad, cerad_df['dpt_pseudotime'])
                    cerad_or = np.exp(cerad_model.params['dpt_pseudotime'])
                    cerad_p = cerad_model.pvalues['dpt_pseudotime']
                    cerad_conf = cerad_model.conf_int().loc['dpt_pseudotime']
                    cerad_or_lower = np.exp(cerad_conf[0])
                    cerad_or_upper = np.exp(cerad_conf[1])

                    cerad_aucs.append(cerad_auc)
                    cerad_odds_ratios.append(cerad_or)
                    cerad_pvals.append(cerad_p)
                    cerad_or_cis.append((cerad_or_lower, cerad_or_upper))

                    fpr_cerad, tpr_cerad, _ = roc_curve(y_cerad, cerad_df['dpt_pseudotime'])

                    if cerad_auc > 0.55:
                        roc_plot_data_cerad.append({
                            'fpr': fpr_cerad,
                            'tpr': tpr_cerad,
                            'auc': cerad_auc,
                            'seed': seed
                        })

                    print(f"  CERAD AUC (seed {seed}): {cerad_auc:.3f}")
            except Exception as e:
                print(f"  CERAD model failed for seed {seed}: {e}")

        # --- After all seeds are run, aggregate results and create combined outputs ---
        if not subcluster_aucs:
            print(f"No valid results generated for {current_subcluster}. Skipping final aggregation.")
            continue

        # === 1) Aggregate pseudotime across seeds (Option A: rank-normalized) ===
        if pseudotime_per_seed:
            pseudo_array = np.column_stack(pseudotime_per_seed)  # cells x n_seeds_used
            n_cells, n_seeds_used = pseudo_array.shape
            ranked = np.zeros_like(pseudo_array, dtype=float)

            for j in range(n_seeds_used):
                col = pseudo_array[:, j]
                valid = ~np.isnan(col)
                ranks = np.empty_like(col, dtype=float)
                ranks[~valid] = np.nan

                # rank 0..(n_valid-1)
                valid_vals = col[valid]
                order = np.argsort(valid_vals)
                rank_values = np.empty_like(order, dtype=float)
                rank_values[order] = np.arange(len(order), dtype=float)

                ranks[valid] = rank_values
                denom = max(len(order) - 1, 1)
                ranked[:, j] = ranks / denom

            mean_pseudotime = np.nanmean(ranked, axis=1)
            adata_template.obs['mean_pseudotime'] = mean_pseudotime
            adata_template.obs['TAG'] = adata_template.obs.index

            # --- Save GAM input using mean pseudotime ---
            # Mathys 'projid' -> SEA-AD 'Donor ID' (subject identifier)
            gam_df = adata_template.obs[['Donor ID', 'AD', 'TAG']].copy()
            gam_df['AD_binary'] = (gam_df['AD'] == 'AD').astype(int)
            gam_df['dpt_pseudotime'] = mean_pseudotime
            gam_df = gam_df.dropna(subset=['Donor ID', 'AD_binary', 'TAG', 'dpt_pseudotime'])
            gam_df.to_csv(os.path.join(output_path, f"pseudotime_gam_input_{safe_subcluster_name}.csv"), index=False)

            # === 2) PC selection by correlation with mean pseudotime ===
            pc_scores = adata_template.obsm['X_pca'][:, :n_pcs]  # cells x PCs
            pc_correlations = []
            for k in range(pc_scores.shape[1]):
                rho, _ = spearmanr(pc_scores[:, k], mean_pseudotime)
                pc_correlations.append(rho)
            best_pc = int(np.nanargmax(np.abs(pc_correlations)))

            # "Average" gene loadings across seeds (PCA is shared; average = original)
            loadings_matrix = adata_template.varm['PCs']  # genes x PCs
            avg_loadings = loadings_matrix[:, best_pc]

            top_idx = np.argsort(np.abs(avg_loadings))[::-1][:10]
            top_trajectory_genes_from_gt = adata_template.var_names[top_idx].tolist()
            print(f"  Top trajectory-informing genes (PC{best_pc + 1}): {top_trajectory_genes_from_gt}")

            top_50 = np.argsort(np.abs(avg_loadings))[::-1][:50]
            top_50_trajectory_genes_from_gt = adata_template.var_names[top_50].tolist()

            # --- Save top 50 gene loadings in rank order ---
            top50_df = pd.DataFrame({
                'rank': np.arange(1, 51),
                'gene': adata_template.var_names[top_50],
                'loading': avg_loadings[top_50],
                'abs_loading': np.abs(avg_loadings[top_50])
            })

            top50_path = os.path.join(output_path, f"top50_PC_loadings_{safe_subcluster_name}.csv")
            top50_df.to_csv(top50_path, index=False)
            print(f"Saved top-50 PC loadings -> {top50_path}")

            # === 3a) UMAP plot colored by mean pseudotime ===
            plot_filename_umap = os.path.join(output_path, f"umap_pseudotime_{safe_subcluster_name}.png")
            sc.pl.umap(
                adata_template,
                color='mean_pseudotime',
                title=f'{current_subcluster}: Mean Pseudotime Trajectory',
                show=False,
                save=f"_{safe_subcluster_name}_umap.png",
                cmap="viridis"
            )
            shutil.move(f'figures/umap_{safe_subcluster_name}_umap.png', plot_filename_umap)

            # === 3b) Smoothed expression of top genes vs mean pseudotime ===
            if top_trajectory_genes_from_gt:
                fig, ax = plt.subplots(figsize=(10, 6))
                genes_to_plot = top_trajectory_genes_from_gt[:4]
                colors = sns.color_palette("colorblind", n_colors=len(genes_to_plot))

                pseudo_vals = mean_pseudotime
                valid_idx = ~np.isnan(pseudo_vals)
                sorted_idx = np.argsort(pseudo_vals[valid_idx])

                for gene, color in zip(genes_to_plot, colors):
                    expr_vals = adata_template[:, gene].X
                    # Handle sparse
                    if not isinstance(expr_vals, np.ndarray):
                        expr_vals = expr_vals.toarray()
                    expr_vals = np.asarray(expr_vals).flatten()

                    expr_valid = expr_vals[valid_idx][sorted_idx]
                    pseudo_sorted = pseudo_vals[valid_idx][sorted_idx]

                    smoothed = lowess(expr_valid, pseudo_sorted, frac=0.3)
                    smoothed_z = (smoothed[:, 1] - smoothed[:, 1].mean()) / smoothed[:, 1].std()
                    ax.plot(smoothed[:, 0], smoothed_z, label=gene, color=color)

                ax.set_xlabel('Mean Pseudotime (rank-normalized)')
                ax.set_ylabel('Z-scored Expression')
                ax.set_title(f'{current_subcluster}: Smoothed Expression of Top Genes')
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                fig.savefig(os.path.join(output_path, f"lineplot_top_genes_{safe_subcluster_name}.png"),
                            bbox_inches='tight')
                plt.close(fig)

        # === Aggregate CERAD stats across seeds ===
        if cerad_aucs:
            mean_auc_cerad = np.mean(cerad_aucs)
            se_auc_cerad = np.std(cerad_aucs, ddof=1) / np.sqrt(len(cerad_aucs))
            auc_cerad_ci = (mean_auc_cerad - 1.96 * se_auc_cerad, mean_auc_cerad + 1.96 * se_auc_cerad)

            log_cerad_ors = np.log(cerad_odds_ratios)
            mean_log_or_cerad = np.mean(log_cerad_ors)
            se_log_or_cerad = np.std(log_cerad_ors, ddof=1) / np.sqrt(len(log_cerad_ors))
            or_cerad_ci = (
                np.exp(mean_log_or_cerad - 1.96 * se_log_or_cerad),
                np.exp(mean_log_or_cerad + 1.96 * se_log_or_cerad)
            )
            mean_or_cerad = np.exp(mean_log_or_cerad)

            _, combined_p_cerad = combine_pvalues(cerad_pvals, method='fisher')

            print(f"  CERAD mean AUC: {mean_auc_cerad:.3f} (95% CI: {auc_cerad_ci[0]:.3f}-{auc_cerad_ci[1]:.3f})")
            print(f"  CERAD mean OR: {mean_or_cerad:.2f} (95% CI: {or_cerad_ci[0]:.2f}-{or_cerad_ci[1]:.2f}), p={combined_p_cerad:.3g}")

            if mean_auc_cerad > 0.55:
                mean_fpr_cerad = np.linspace(0, 1, 100)
                interp_tprs_cerad = []

                for data in roc_plot_data_cerad:
                    interp_tprs_cerad.append(np.interp(mean_fpr_cerad, data['fpr'], data['tpr']))

                mean_tpr_cerad = np.mean(interp_tprs_cerad, axis=0)
                mean_tpr_cerad[0], mean_tpr_cerad[-1] = 0, 1

                summary_roc_data_cerad.append({
                    'fpr': mean_fpr_cerad,
                    'tpr': mean_tpr_cerad,
                    'Subcluster': current_subcluster,
                    'mean_auc': mean_auc_cerad,
                    'se_auc': se_auc_cerad
                })
        else:
            mean_auc_cerad = np.nan
            auc_cerad_ci = (np.nan, np.nan)
            mean_or_cerad = np.nan
            or_cerad_ci = (np.nan, np.nan)
            combined_p_cerad = np.nan

        # === Aggregate clinical AD stats across seeds ===
        mean_auc = np.mean(subcluster_aucs)
        se_auc = np.std(subcluster_aucs, ddof=1) / np.sqrt(len(subcluster_aucs))
        auc_ci_95 = (mean_auc - 1.96 * se_auc, mean_auc + 1.96 * se_auc)

        log_ors = np.log(subcluster_odds_ratios)
        mean_log_or = np.mean(log_ors)
        se_log_or = np.std(log_ors, ddof=1) / np.sqrt(len(log_ors))
        log_or_ci_95 = (mean_log_or - 1.96 * se_log_or, mean_log_or + 1.96 * se_log_or)
        mean_or = np.exp(mean_log_or)
        or_ci_95 = (np.exp(log_or_ci_95[0]), np.exp(log_or_ci_95[1]))

        _, combined_p_value = combine_pvalues(subcluster_pvals, method='fisher')

        print(f"\n--- Aggregated Results for {current_subcluster} ---")
        print(f"Mean AUC: {mean_auc:.4f} (95% CI: {auc_ci_95[0]:.4f} - {auc_ci_95[1]:.4f})")
        print(f"Mean Odds Ratio: {mean_or:.4f} (95% CI: {or_ci_95[0]:.4f} - {or_ci_95[1]:.4f})")
        print(f"Combined P-value (Fisher's): {combined_p_value:.4g}")

        # === Mean ROC curve across seeds for this subcluster ===
        if mean_auc > 0.55:
            mean_fpr = np.linspace(0, 1, 100)
            interp_tprs = []
            for fpr, tpr in zip(all_fprs, all_tprs):
                interp_tprs.append(np.interp(mean_fpr, fpr, tpr))

            if interp_tprs:
                mean_tpr = np.mean(interp_tprs, axis=0)
                mean_tpr[0], mean_tpr[-1] = 0, 1

                summary_roc_data.append({
                    'fpr': mean_fpr,
                    'tpr': mean_tpr,
                    'Subcluster': current_subcluster,
                    'mean_auc': mean_auc,
                    'se_auc': se_auc
                })

        # Store per-subcluster summary row
        validation_results.append({
            'Subcluster': current_subcluster,
            'mean_auc': mean_auc,
            'auc_ci_lower': auc_ci_95[0],
            'auc_ci_upper': auc_ci_95[1],
            'mean_odds_ratio': mean_or,
            'or_ci_lower': or_ci_95[0],
            'or_ci_upper': or_ci_95[1],
            'combined_p_value': combined_p_value,
            'mean_auc_cerad': mean_auc_cerad,
            'auc_cerad_lower': auc_cerad_ci[0],
            'auc_cerad_upper': auc_cerad_ci[1],
            'odds_ratio_cerad': mean_or_cerad,
            'or_cerad_lower': or_cerad_ci[0],
            'or_cerad_upper': or_cerad_ci[1],
            'p_value_cerad': combined_p_cerad,
            'n_cells': adata_template.n_obs,
            'top_trajectory_genes': ', '.join(top_trajectory_genes_from_gt) if top_trajectory_genes_from_gt else "N/A",
            'top_pc_loadings': ', '.join(top_50_trajectory_genes_from_gt) if top_50_trajectory_genes_from_gt else "N/A",
        })

    # --- Save per-subcluster outputs (this job's piece of the final aggregation) ---
    # The Mathys notebook accumulated across all subclusters in one process and
    # wrote a single publication_summary_statistics.csv + combined ROC plots at
    # the end. With per-subcluster SLURM parallelization, each job writes its
    # own row/ROC arrays; a separate aggregator combines them after all jobs run.
    if not validation_results:
        print(f"\nNo valid results for {args.subcluster}. Nothing to save.")
    else:
        safe = re.sub(r'[^a-zA-Z0-9_.-]', '_', args.subcluster)

        row_df = pd.DataFrame(validation_results)
        row_path = os.path.join(output_path, f"summary_row_{safe}.csv")
        row_df.to_csv(row_path, index=False)
        print(f"Saved per-subcluster summary row -> {row_path}")

        if summary_roc_data:
            d = summary_roc_data[0]
            roc_clin_path = os.path.join(output_path, f"roc_clinical_{safe}.npz")
            np.savez(roc_clin_path,
                     fpr=d['fpr'], tpr=d['tpr'],
                     mean_auc=d['mean_auc'], se_auc=d['se_auc'],
                     subcluster=current_subcluster)
            print(f"Saved per-subcluster clinical ROC -> {roc_clin_path}")

        if summary_roc_data_cerad:
            d = summary_roc_data_cerad[0]
            roc_cerad_path = os.path.join(output_path, f"roc_cerad_{safe}.npz")
            np.savez(roc_cerad_path,
                     fpr=d['fpr'], tpr=d['tpr'],
                     mean_auc=d['mean_auc'], se_auc=d['se_auc'],
                     subcluster=current_subcluster)
            print(f"Saved per-subcluster CERAD ROC -> {roc_cerad_path}")

    print("\n--- Analysis complete. ---")
    print("\n--- Pseudotime robustness analysis complete for all subclusters. ---")


if __name__ == '__main__':
    main()
