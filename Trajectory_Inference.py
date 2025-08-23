import scanpy as sc
import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, roc_curve
from scipy.stats import mannwhitneyu # NEW: Library for statistical testing on violin plots
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests
import scanpy.external as sce
import seaborn as sns
# from gprofiler import GProfiler
# from gseapy import enrichr
# gp = GProfiler(return_dataframe=True)

# --- Set up Paths ---
# Define paths to your data files
expr_path = "/Users/adithyamadduri/Desktop/Projects/Patel_Lab/CellMatrix_with_genenames.parquet"
meta_path = "/Users/adithyamadduri/Desktop/Projects/Patel_Lab/New_CellMetadataSyn1848517.parquet"
deg_path = "/Users/adithyamadduri/Desktop/Projects/Patel_Lab/DE_Outputs_Fixed/"
# deg_path = "/Users/adithyamadduri/Desktop/Projects/Patel_Lab/Clincal_batch_DE_Outputs/"
output_path = "/Users/adithyamadduri/Desktop/Projects/Patel_Lab/Pseudotime_Outputs_Final/" # MODIFIED: Changed output path for new results

# Create output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)

# --- Load and Prepare Data ---
# Read the metadata file
meta_df = pd.read_parquet(meta_path)

# --- Define AD status based on dcfdx ---
meta_df['AD'] = meta_df['dcfdx_lv'].isin([4.0, 5.0]).map({True: 'AD', False: 'Control'})
meta_df['AD'] = meta_df['AD'].astype('category')
print("AD status defined based on 'dcfdx':")
print(meta_df['AD'].value_counts())

# Read the expression matrix (genes x cells)
expr_df = pd.read_parquet(expr_path)
# Set the gene names as the index
expr_df.set_index('index', inplace=True)
# Transpose the matrix to get cells x genes format, which scanpy expects
expr_df_transposed = expr_df.T 

# --- Main Analysis Loop for Each Subcluster ---
# Using capitalized 'Subcluster' as requested
all_subclusters = meta_df['Subcluster'].unique()
# To test on a single subcluster, you can uncomment the line below:
# all_subclusters = ['Ex6']

# NEW: Define the list of random seeds for robustness check
random_seeds = [41, 42, 52, 57, 62, 72, 73]
ground_truth_seed = 42 # The seed for generating all standard plots

# Create a list to store the new validation results
validation_results = []

summary_roc_data = []

summary_roc_data_cerad = []

roc_plot_data_cerad = []

from gseapy import enrichr



# def run_pseudotime_enrichment(adata, current_subcluster, safe_subcluster_name, output_path, top_n_terms=10, celltype_degs=None):
#     adata.obs['quantile'] = pd.qcut(
#         adata.obs['dpt_pseudotime'],
#         q=5,
#         labels=[f'quantile_{i}' for i in range(5)]
#     ).astype('category')

#     sc.tl.rank_genes_groups(
#         adata,
#         groupby='quantile',
#         groups=['quantile_4'],
#         reference='quantile_0',
#         method='wilcoxon'
#     )

    # de_df = sc.get.rank_genes_groups_df(adata, group='quantile_4')
    # deg_genes = de_df[
    #     (de_df['pvals_adj'] < 0.05) & 
    #     (de_df['logfoldchanges'].abs() > 0.25)
    # ]['names'].dropna().astype(str).str.upper().tolist()

    # if celltype_degs is not None: # keep only those also in your parent‐DEG list
    #     deg_genes = [g for g in deg_genes if g in celltype_degs]

    # print(f"[{current_subcluster}] Late>Early DEGs passing FDR<0.05 & logFC>0: {len(deg_genes)}")
    # if len(deg_genes) == 0:
    #     print(f"[{current_subcluster}] No significant late-up genes. Skipping enrichment.")
    #     return

    # enrichr_output = os.path.join(output_path, f"gprofiler_results_{safe_subcluster_name}")
    # os.makedirs(enrichr_output, exist_ok=True)

    # # Run g:Profiler with GO:BP and DisGeNET
    # result_df = gp.profile(
    #     organism="hsapiens",
    #     query=deg_genes,
    #     user_threshold=0.05,
    #     sources=['GO:BP', 'KEGG', 'REAC', 'WP', 'HP'],
    #     no_evidences=False,
    #     background=celltype_degs  # Optional: only if background list is desired
    # )

    # if result_df.empty:
    #     print(f"[{current_subcluster}] No enrichment results returned by g:Profiler.")
    #     return

    # # Save full results
    # result_df.to_csv(os.path.join(enrichr_output, "gprofiler_results.csv"), index=False)

    # # keep only GO:BP and DisGeNET
    # keep = ['HP']

    # df = result_df[result_df["source"].isin(keep)].copy()

    # df["Term_label"] = df["name"]
    

    # # pick the top N across both sets by Adjusted P
    # top_df = df.sort_values("p_value").head(top_n_terms)
    # top_df["Term_label"].astype(str)

    # # single bar‐plot
    # plt.figure(figsize=(8,6))
    # plt.barh(
    # top_df["Term_label"][::-1],
    # -np.log10(top_df["p_value"][::-1]),
    # color="steelblue"
    # )
    # plt.xlabel("-log10 Adjusted P-value")
    # plt.title(f"{current_subcluster}: GO & DisGeNET enrichment")
    # plt.tight_layout()
    # plt.savefig(os.path.join(enrichr_output, f"GO_DisGeNET_barplot_top{top_n_terms}.png"), dpi=300)
    # plt.show()

# 2. Loop through each unique subcluster
for current_subcluster in all_subclusters:
    print(f"\n--- Processing Subcluster: {current_subcluster} ---")

    # --- Subset Data for the Current Subcluster ---
    meta_sub = meta_df[meta_df['Subcluster'] == current_subcluster]

    if meta_sub.shape[0] < 50 or len(meta_sub['AD'].unique()) < 2:
        print(f"Skipping {current_subcluster} due to having fewer than 50 cells or only one group (AD/Control).")
        continue
        
    safe_subcluster_name = re.sub(r'[^a-zA-Z0-9_.-]', '_', current_subcluster)
    subcluster_cells = meta_sub['TAG'].tolist()
    expr_sub = expr_df_transposed.loc[subcluster_cells]

    # --- Create AnnData Object and perform initial processing ---
    adata_template = sc.AnnData(X=expr_sub)
    adata_template.obs = meta_sub.set_index('TAG').loc[adata_template.obs.index]
    
    print("Loading DEGs with p_adj < 0.05...")
    parent_cell_type = adata_template.obs['broad.cell.type'].iloc[0]
    deg_file_path = os.path.join(deg_path, f"poisson_DE_results_{parent_cell_type}.csv")

    if not os.path.exists(deg_file_path):
        print(f"Warning: DEG file not found for parent type {parent_cell_type}. Skipping.")
        continue

    degs_df = pd.read_csv(deg_file_path)
    genes_for_trajectory = degs_df[degs_df['p_adj'] < 0.05]['gene'].tolist()
    available_genes = [gene for gene in genes_for_trajectory if gene in adata_template.var_names]
    
    if len(available_genes) < 10:
        print(f"Skipping {current_subcluster} as there are fewer than 10 DEGs available.")
        continue
        
    adata_template = adata_template[:, available_genes].copy()
    print(f"Using {len(available_genes)} genes to build trajectories.")

    sc.pp.normalize_total(adata_template, target_sum=1e4)
    sc.pp.log1p(adata_template)
    # sc.pp.pca(adata_template, n_comps=min(50, len(available_genes) - 1))
    # sc.pp.neighbors(adata_template, n_neighbors=10, n_pcs=min(50, len(available_genes) - 1))


    # Run PCA first (required before Harmony)
    sc.pp.pca(adata_template, n_comps=min(50, len(available_genes) - 1))


    adata_template.obs['projid'] = adata_template.obs['projid'].astype('category')
    sce.pp.harmony_integrate(adata_template, key='projid')
    adata_template.obsm['X_pca'] = adata_template.obsm['X_pca_harmony']  # overwrite PCA with corrected version

    # Run Harmony on PCA (batch = projid)
    # sce.pp.harmony_integrate(adata_template, key='projid')

    # Then build neighbors using Harmony-corrected PCA
    sc.pp.neighbors(adata_template, use_rep='X_pca')


    sc.tl.leiden(adata_template)
    sc.tl.umap(adata_template)
    sc.tl.paga(adata_template, groups='leiden')
    sc.tl.diffmap(adata_template)

    # --- NEW: Loop over random seeds to test trajectory robustness ---
    subcluster_aucs = []
    roc_plot_data = []
    top_trajectory_genes_from_gt = [] # To store genes from the ground truth run

    subcluster_pvals, subcluster_odds_ratios = [], []
    all_fprs, all_tprs = [], [] # To calculate mean ROC

    cerad_aucs = []
    cerad_odds_ratios = []
    cerad_pvals = []
    cerad_or_cis = []

    for seed in random_seeds:
        print(f"\n  -- Running analysis for seed: {seed} --")
        
        # Create a fresh copy for each run to avoid overwriting results
        adata = adata_template.copy()

        # --- Select Root Cell Based on ceradsc and cogdx (with current seed) ---
        subcluster_meta_dpt = adata.obs.copy()
        subcluster_meta_dpt['ceradsc'] = pd.to_numeric(subcluster_meta_dpt['ceradsc'], errors='coerce')
        subcluster_meta_dpt['cogdx'] = pd.to_numeric(subcluster_meta_dpt['cogdx'], errors='coerce')
        subcluster_meta_dpt.dropna(subset=['ceradsc', 'cogdx'], inplace=True)
        
        if subcluster_meta_dpt.empty:
            print("  Skipping trajectory ordering: No cells with valid ceradsc/cogdx values.")
            # Break the inner loop since this will be true for all seeds
            break

        max_ceradsc = subcluster_meta_dpt['ceradsc'].max()
        root_candidates = subcluster_meta_dpt[subcluster_meta_dpt['ceradsc'] == max_ceradsc]
        min_cogdx = root_candidates['cogdx'].min()
        root_candidates = root_candidates[root_candidates['cogdx'] == min_cogdx]
        
        # MODIFIED: Use the current seed for sampling
        root_cell_id = root_candidates.sample(n=1, random_state=seed).index[0]
        
        adata.uns['iroot'] = np.where(adata.obs_names == root_cell_id)[0][0]

        # --- Calculate Pseudotime (DPT) ---
        sc.tl.dpt(adata)

        

        # --- Calculate AUC for the current seed ---
        validation_df = adata.obs[['dpt_pseudotime', 'AD']].copy()
        validation_df.dropna(subset=['dpt_pseudotime'], inplace=True) # Ensure no NaNs in pseudotime
        y_true = (validation_df['AD'] == 'AD').astype(int)
        y_score = validation_df['dpt_pseudotime']
        
        auc_score = roc_auc_score(y_true=y_true, y_score=y_score)
        subcluster_aucs.append(auc_score)
        print(f"  --> Seed {seed}: AUC = {auc_score:.4f}")

        # Run logistic regression for this seed
        X = sm.add_constant(validation_df['dpt_pseudotime'])
        if len(y_true.unique()) < 2:
            print(f"  Skipping seed {seed}: only one class present after dropping NaNs.")
            continue # Skip this seed if only one outcome is present
            
        log_reg = sm.Logit(y_true, X).fit(disp=0)
        
        # Store the stats from this run
        subcluster_pvals.append(log_reg.pvalues['dpt_pseudotime'])
        subcluster_odds_ratios.append(np.exp(log_reg.params['dpt_pseudotime']))

        # --- Store data for combined ROC plot if AUC > 0.6 ---
        # Calculate and store the ROC curve for this seed
        fpr, tpr, _ = roc_curve(y_true=y_true, y_score=y_score)
        all_fprs.append(fpr)
        all_tprs.append(tpr)

        # Store data for the per-subcluster plot if AUC > 0.6
        if auc_score > 0.6:
            roc_plot_data.append({'fpr': fpr, 'tpr': tpr, 'auc': auc_score, 'seed': seed})

        
        # CERAD-based classifier (moved outside ground_truth_seed)
        try:
            cerad_df = adata.obs[['dpt_pseudotime']].copy()
            cerad_df['TAG'] = cerad_df.index

            meta_subset = meta_df[['TAG', 'ceradsc']].dropna()
            meta_subset['ceradsc'] = pd.to_numeric(meta_subset['ceradsc'], errors='coerce')

            cerad_df = cerad_df.merge(meta_subset, on="TAG", how="left")
            cerad_df = cerad_df.dropna(subset=["ceradsc", "dpt_pseudotime"])
            cerad_df['AD_cerad'] = cerad_df['ceradsc'].isin([1, 2]).astype(int)

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

                if cerad_auc > 0.6:
                    roc_plot_data_cerad.append({
                        'fpr': fpr_cerad,
                        'tpr': tpr_cerad,
                        'auc': cerad_auc,
                        'seed': seed
                    })

                print(f"  CERAD AUC (seed {seed}): {cerad_auc:.3f}")
        except Exception as e:
            print(f"  CERAD model failed for seed {seed}: {e}")

        # --- MODIFIED: Generate plots ONLY for the "ground truth" seed (42) ---
        if seed == ground_truth_seed:
            print(f"  --- Generating plots for ground truth seed: {ground_truth_seed} ---")

            # run_pseudotime_enrichment(
            #     adata=adata,
            #     current_subcluster=current_subcluster,
            #     safe_subcluster_name=safe_subcluster_name,
            #     output_path=output_path,
            #     top_n_terms=10,
            #     celltype_degs=available_genes
            # )
            adata.obs['TAG'] = adata.obs.index
            subset_df = adata.obs[['projid', 'AD', 'dpt_pseudotime', 'TAG']].copy()
            subset_df['AD_binary'] = (subset_df['AD'] == 'AD').astype(int)
            
            subset_df = subset_df.dropna(subset=['projid', 'AD_binary', 'TAG', 'dpt_pseudotime'])

            # Save to CSV for import into R
            subset_df.to_csv(os.path.join(output_path, f"pseudotime_gam_input_{safe_subcluster_name}.csv"), index=False)
            
            # --- Identify Top Genes from PCA Loadings (from ground_truth run) ---
            pc4_loadings = adata.varm['PCs'][:, 3]
            gene_loading_df = pd.DataFrame({
                'gene': adata.var_names,
                'loading_abs': np.abs(pc4_loadings)
            }).sort_values(by='loading_abs', ascending=False)
            top_trajectory_genes_from_gt = gene_loading_df.head(8)['gene'].tolist()
            print(f"  Top trajectory-informing genes: {top_trajectory_genes_from_gt}")

            # Plot 1: UMAP colored by pseudotime
            plot_filename_umap = os.path.join(output_path, f"umap_pseudotime_{safe_subcluster_name}.png")
            sc.pl.umap(
                adata,
                color='dpt_pseudotime',
                title=f'{current_subcluster}: Pseudotime Trajectory',
                show=False,
                save=f"_{safe_subcluster_name}_umap.png",
                cmap="viridis"  # or "plasma", "turbo", etc.
            )
            
            #sc.pl.umap(adata, color='dpt_pseudotime', title=f'{current_subcluster}: Pseudotime Trajectory (Seed {seed})', show=False, save=f"_{safe_subcluster_name}_umap.png")
            os.rename(f'figures/umap_{safe_subcluster_name}_umap.png', plot_filename_umap)

            # Plot 2: ENHANCED Violin plot of pseudotime vs AD status
            plot_filename_violin = os.path.join(output_path, f"violin_pseudotime_{safe_subcluster_name}.png")
            fig, ax = plt.subplots(figsize=(7, 6))
            sc.pl.violin(adata, keys='dpt_pseudotime', groupby='AD', ax=ax, show=False)
            ax.set_ylabel("Pseudotime")
            ax.set_xlabel("Group")
            
            # Add median and statistical test
            ad_pseudo = adata.obs[adata.obs['AD'] == 'AD']['dpt_pseudotime'].dropna()
            ctrl_pseudo = adata.obs[adata.obs['AD'] == 'Control']['dpt_pseudotime'].dropna()
            median_ad = ad_pseudo.median()
            median_ctrl = ctrl_pseudo.median()
            stat, p_val = mannwhitneyu(ad_pseudo, ctrl_pseudo, alternative='two-sided')
            
            ax.set_title(f'Pseudotime Distribution in {current_subcluster}')
            ax.text(0.5, 0.95, f"Median Diff: {median_ad - median_ctrl:.2f}\nMann-Whitney U p-val: {p_val:.2e}", 
                    transform=ax.transAxes, ha='center', va='top', bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))
            
            fig.savefig(plot_filename_violin, bbox_inches='tight')
            plt.close(fig)

            # Plot 3: Smoothed expression of TOP TRAJECTORY-INFORMING genes
            if top_trajectory_genes_from_gt:
                fig, ax = plt.subplots(figsize=(10, 6))

                genes_to_plot = top_trajectory_genes_from_gt[:4]

                # Color-blind friendly palette
                colors = sns.color_palette("colorblind", n_colors=len(genes_to_plot))

                # for gene in top_trajectory_genes_from_gt:
                #     expr_vals = adata[:, gene].X.toarray().flatten()
                #     pseudo_vals = adata.obs['dpt_pseudotime'].values
                #     # Filter out NaNs in pseudotime for LOWESS
                #     valid_idx = ~np.isnan(pseudo_vals)
                #     sorted_idx = np.argsort(pseudo_vals[valid_idx])
                #     smoothed = lowess(expr_vals[valid_idx][sorted_idx], pseudo_vals[valid_idx][sorted_idx], frac=0.3)
                #     smoothed_z = (smoothed[:, 1] - smoothed[:, 1].mean()) / smoothed[:, 1].std()
                #     ax.plot(smoothed[:, 0], smoothed_z, label=gene)

                # ax.set_xlabel('Pseudotime')
                # ax.set_ylabel('Z-scored Expression')
                # ax.set_title(f'{current_subcluster}: Smoothed Expression of Top Genes')
                # ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                # fig.savefig(os.path.join(output_path, f"lineplot_top_genes_{safe_subcluster_name}.png"), bbox_inches='tight')
                # plt.close(fig)

                for gene, color in zip(genes_to_plot, colors):
                    expr_vals = adata[:, gene].X.toarray().flatten()
                    pseudo_vals = adata.obs['dpt_pseudotime'].values
                    valid_idx = ~np.isnan(pseudo_vals)
                    sorted_idx = np.argsort(pseudo_vals[valid_idx])
                    smoothed = lowess(expr_vals[valid_idx][sorted_idx],
                                    pseudo_vals[valid_idx][sorted_idx],
                                    frac=0.3)
                    smoothed_z = (smoothed[:, 1] - smoothed[:, 1].mean()) / smoothed[:, 1].std()
                    ax.plot(smoothed[:, 0], smoothed_z, label=gene, color=color)

                ax.set_xlabel('Pseudotime')
                ax.set_ylabel('Z-scored Expression')
                ax.set_title(f'{current_subcluster}: Smoothed Expression of Top Genes')
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                fig.savefig(os.path.join(output_path, f"lineplot_top_genes_{safe_subcluster_name}.png"),
                            bbox_inches='tight')
                plt.close(fig)

    # --- After all seeds are run, aggregate results and create combined ROC plot ---
    if not subcluster_aucs:
        print(f"No valid results generated for {current_subcluster}. Skipping final aggregation.")
        continue


    if cerad_aucs:
        mean_auc_cerad = np.mean(cerad_aucs)
        se_auc_cerad = np.std(cerad_aucs, ddof=1) / np.sqrt(len(cerad_aucs))
        auc_cerad_ci = (mean_auc_cerad - 1.96 * se_auc_cerad, mean_auc_cerad + 1.96 * se_auc_cerad)

        log_cerad_ors = np.log(cerad_odds_ratios)
        mean_log_or_cerad = np.mean(log_cerad_ors)
        se_log_or_cerad = np.std(log_cerad_ors, ddof=1) / np.sqrt(len(log_cerad_ors))
        or_cerad_ci = (np.exp(mean_log_or_cerad - 1.96 * se_log_or_cerad), np.exp(mean_log_or_cerad + 1.96 * se_log_or_cerad))
        mean_or_cerad = np.exp(mean_log_or_cerad)

        _, combined_p_cerad = combine_pvalues(cerad_pvals, method='fisher')

        print(f"  CERAD mean AUC: {mean_auc_cerad:.3f} (95% CI: {auc_cerad_ci[0]:.3f}–{auc_cerad_ci[1]:.3f})")
        print(f"  CERAD mean OR: {mean_or_cerad:.2f} (95% CI: {or_cerad_ci[0]:.2f}–{or_cerad_ci[1]:.2f}), p={combined_p_cerad:.3g}")

        if mean_auc_cerad > 0.6:
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

    # NEW: Calculate mean and standard error for AUC
    # NEW: Calculate aggregated statistics with confidence intervals
    # Calculate AUC stats
    mean_auc = np.mean(subcluster_aucs)
    se_auc = np.std(subcluster_aucs, ddof=1) / np.sqrt(len(subcluster_aucs))
    auc_ci_95 = (mean_auc - 1.96 * se_auc, mean_auc + 1.96 * se_auc)
    
    # Calculate Odds Ratio stats (using log transform for stability)
    log_ors = np.log(subcluster_odds_ratios)
    mean_log_or = np.mean(log_ors)
    se_log_or = np.std(log_ors, ddof=1) / np.sqrt(len(log_ors))
    log_or_ci_95 = (mean_log_or - 1.96 * se_log_or, mean_log_or + 1.96 * se_log_or)
    mean_or = np.exp(mean_log_or)
    or_ci_95 = (np.exp(log_or_ci_95[0]), np.exp(log_or_ci_95[1]))

    # Combine p-values using Fisher's method
    _, combined_p_value = combine_pvalues(subcluster_pvals, method='fisher')
    
    print(f"\n--- Aggregated Results for {current_subcluster} ---")
    print(f"Mean AUC: {mean_auc:.4f} (95% CI: {auc_ci_95[0]:.4f} - {auc_ci_95[1]:.4f})")
    print(f"Mean Odds Ratio: {mean_or:.4f} (95% CI: {or_ci_95[0]:.4f} - {or_ci_95[1]:.4f})")
    print(f"Combined P-value (Fisher's): {combined_p_value:.4g}")

    # NEW: Calculate and store data for the final summary ROC plot
    if mean_auc > 0.6:
        # We need to interpolate to get a "mean" curve
        mean_fpr = np.linspace(0, 1, 100)
        interp_tprs = []
        # We need the original fpr/tpr from each seed run for this
        for fpr, tpr in zip(all_fprs, all_tprs):
             interp_tprs.append(np.interp(mean_fpr, fpr, tpr))
        
        if interp_tprs:
            mean_tpr = np.mean(interp_tprs, axis=0)
            mean_tpr[0], mean_tpr[-1] = 0, 1 # Ensure curve starts/ends correctly
            
            summary_roc_data.append({
                'fpr': mean_fpr,
                'tpr': mean_tpr,
                'Subcluster': current_subcluster,
                'mean_auc': mean_auc,
                'se_auc': se_auc
            })

        # === CERAD Combined ROC Plot ===
    if summary_roc_data_cerad:
        plt.style.use('seaborn-v0_8-whitegrid')
        fig, ax = plt.subplots(figsize=(10, 8))

        summary_roc_data_cerad.sort(key=lambda x: x['mean_auc'], reverse=True)

        for data in summary_roc_data_cerad:
            ax.plot(data['fpr'], data['tpr'],
                    label=f"{data['Subcluster']} (AUC = {data['mean_auc']:.3f} ± {data['se_auc']:.3f})")

        ax.plot([0, 1], [0, 1], 'k--', label='Chance (AUC = 0.50)')
        ax.set_xlabel('False Positive Rate', fontsize=14)
        ax.set_ylabel('True Positive Rate', fontsize=14)
        ax.set_title('Pathological AD ROC Curves using Disease Trajectory', fontsize=16)
        ax.legend(loc='lower right', frameon=True, fontsize=10)
        plt.savefig(os.path.join(output_path, "summary_roc_plot_CERAD_all_subclusters.png"),
                    bbox_inches='tight', dpi=300)
        plt.close(fig)
        print("Saved CERAD summary ROC plot to: summary_roc_plot_CERAD_all_subclusters.png") 

    # NEW: Generate and save the combined ROC curve plot
    if roc_plot_data:
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.plot([0, 1], [0, 1], 'k--', label='Chance (AUC = 0.50)') # Chance line
        
        for data in roc_plot_data:
            ax.plot(data['fpr'], data['tpr'], label=f'Seed {data["seed"]} (AUC = {data["auc"]:.3f})')
        
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title(f'ROC Curve Robustness for {current_subcluster}')
        ax.legend(loc='lower right')
        fig.savefig(os.path.join(output_path, f"roc_curves_combined_{safe_subcluster_name}.png"), bbox_inches='tight')
        plt.close(fig)
        print(f"Saved combined ROC plot for {current_subcluster}")

    # MODIFIED: Append aggregated results to the summary list
    # validation_results.append({
    #     'Subcluster': current_subcluster,
    #     'mean_auc': mean_auc,
    #     'auc_ci_lower': auc_ci_95[0],
    #     'auc_ci_upper': auc_ci_95[1],
    #     'mean_odds_ratio': mean_or,
    #     'or_ci_lower': or_ci_95[0],
    #     'or_ci_upper': or_ci_95[1],
    #     'combined_p_value': combined_p_value,
    #     'n_cells': adata.n_obs,
    #     'top_trajectory_genes': ', '.join(top_trajectory_genes_from_gt) if top_trajectory_genes_from_gt else "N/A"
    # })

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
        'n_cells': adata.n_obs,
        'top_trajectory_genes': ', '.join(top_trajectory_genes_from_gt) if top_trajectory_genes_from_gt else "N/A"
    })

# --- Save Final Aggregated Validation Results ---
# --- Final Statistical Correction and Summary ---
if not validation_results:
    print("\nNo subclusters were successfully processed. Exiting.")
else:
    # Create the summary DataFrame and sort by best result
    validation_summary_df = pd.DataFrame(validation_results).sort_values(by='mean_auc', ascending=False)
    
    # NEW: Apply Benjamini-Hochberg FDR correction across all subclusters
    if 'combined_p_value' in validation_summary_df.columns:
        # Get p-values, handling potential NaNs
        pvals = validation_summary_df['combined_p_value'].dropna()
        if not pvals.empty:
            reject, q_values, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
            # Add q-values back to the original DataFrame
            validation_summary_df['fdr_q_value'] = np.nan
            validation_summary_df.loc[pvals.index, 'fdr_q_value'] = q_values
    
    if 'p_value_cerad' in validation_summary_df.columns:
        cerad_pvals = validation_summary_df['p_value_cerad'].dropna()
        if not cerad_pvals.empty:
            _, cerad_qvals, _, _ = multipletests(cerad_pvals, alpha=0.05, method='fdr_bh')
            validation_summary_df['fdr_q_value_cerad'] = np.nan
            validation_summary_df.loc[cerad_pvals.index, 'fdr_q_value_cerad'] = cerad_qvals

    # Save the final, statistically rich summary table
    validation_summary_path = os.path.join(output_path, "publication_summary_statistics.csv")
    validation_summary_df.to_csv(validation_summary_path, index=False)
    print(f"\nFinal publication-ready summary saved to: {validation_summary_path}")

    # NEW: Generate the single, combined ROC plot for all significant subclusters
    if summary_roc_data:
        plt.style.use('seaborn-v0_8-whitegrid') # A nice style for publication
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Sort curves by AUC to make the legend prettier
        summary_roc_data.sort(key=lambda x: x['mean_auc'], reverse=True)
        
        for data in summary_roc_data:
            ax.plot(data['fpr'], data['tpr'], label=f"{data['Subcluster']} (AUC = {data['mean_auc']:.3f} \u00B1 {data['se_auc']:.3f})")

        ax.plot([0, 1], [0, 1], 'k--', label='Chance (AUC = 0.50)')
        ax.set_xlabel('False Positive Rate', fontsize=14)
        ax.set_ylabel('True Positive Rate', fontsize=14)
        ax.set_title('Clinical AD ROC Curves using Disease Trajectory', fontsize=16)
        ax.legend(loc='lower right', frameon=True, fontsize=10)
        plt.savefig(os.path.join(output_path, "summary_roc_plot_all_subclusters.png"), bbox_inches='tight', dpi=300)
        plt.close(fig)
        print(f"Saved summary ROC plot to: {os.path.join(output_path, 'summary_roc_plot_all_subclusters.png')}")

print("\n--- Analysis complete. ---")

print("\n--- Pseudotime robustness analysis complete for all subclusters. ---")
