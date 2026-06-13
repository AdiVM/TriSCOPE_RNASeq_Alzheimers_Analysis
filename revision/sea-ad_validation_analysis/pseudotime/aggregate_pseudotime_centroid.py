"""
aggregate_pseudotime_centroid.py

Aggregator for the per-subcluster pseudotime outputs produced by
pseudotime_seaad_centroid.py.

Reads all summary_row_*.csv (one row per subcluster) and roc_*.npz files,
applies BH-FDR correction to the per-subcluster p-values, writes the
publication summary CSV, and renders the two combined ROC plots (one for
clinical AD, one for CERAD/neuropathology) -- replicating exactly the
aggregation block at the bottom of the Mathys pseudotime_hvg notebook.

Run:
    python aggregate_pseudotime_centroid.py
"""

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output_dir',
                   default='/n/groups/patel/adithya/Alz_Outputs/Final_Outputs/seaad_pseudotime/pseudotime_centroid')
    args = p.parse_args()

    out = args.output_dir

    # -------- 1) Concatenate per-subcluster summary rows --------
    summary_files = sorted(glob.glob(os.path.join(out, 'summary_row_*.csv')))
    if not summary_files:
        raise SystemExit(f'No summary_row_*.csv found under {out}')
    print(f'Concatenating {len(summary_files)} per-subcluster summary rows')
    df = pd.concat([pd.read_csv(f) for f in summary_files], ignore_index=True)
    df = df.sort_values('mean_auc', ascending=False).reset_index(drop=True)

    # -------- 2) BH-FDR correction (clinical AD and CERAD), matches Mathys block --------
    if 'combined_p_value' in df.columns:
        pvals = df['combined_p_value'].dropna()
        if not pvals.empty:
            _, q_values, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
            df['fdr_q_value'] = np.nan
            df.loc[pvals.index, 'fdr_q_value'] = q_values

    if 'p_value_cerad' in df.columns:
        cpvals = df['p_value_cerad'].dropna()
        if not cpvals.empty:
            _, cqvals, _, _ = multipletests(cpvals, alpha=0.05, method='fdr_bh')
            df['fdr_q_value_cerad'] = np.nan
            df.loc[cpvals.index, 'fdr_q_value_cerad'] = cqvals

    summary_path = os.path.join(out, 'publication_summary_statistics.csv')
    df.to_csv(summary_path, index=False)
    print(f'Wrote {summary_path}  ({df.shape[0]} subclusters x {df.shape[1]} cols)')

    # -------- 3) Combined ROC plots --------
    plt.style.use('seaborn-v0_8-whitegrid')

    def plot_combined_roc(pattern, title, out_png):
        files = sorted(glob.glob(os.path.join(out, pattern)))
        if not files:
            print(f'[skip] no files match {pattern} -- no ROC plot rendered.')
            return
        entries = []
        for fp in files:
            d = np.load(fp, allow_pickle=True)
            entries.append({
                'fpr': d['fpr'],
                'tpr': d['tpr'],
                'mean_auc': float(d['mean_auc']),
                'se_auc': float(d['se_auc']),
                'Subcluster': re.sub(r'^roc_(clinical|cerad)_|\.npz$', '', os.path.basename(fp)),
            })
        entries.sort(key=lambda x: x['mean_auc'], reverse=True)

        fig, ax = plt.subplots(figsize=(10, 8))
        for e in entries:
            ax.plot(
                e['fpr'], e['tpr'],
                label=f"{e['Subcluster']} (AUC = {e['mean_auc']:.3f} ± {e['se_auc']:.3f})",
            )
        ax.plot([0, 1], [0, 1], 'k--', label='Chance (AUC = 0.50)')
        ax.set_xlabel('False Positive Rate', fontsize=14)
        ax.set_ylabel('True Positive Rate', fontsize=14)
        ax.set_title(title, fontsize=16)
        ax.legend(loc='lower right', frameon=True, fontsize=10)
        plt.savefig(out_png, bbox_inches='tight', dpi=300)
        plt.close(fig)
        print(f'Wrote {out_png}  ({len(entries)} subclusters)')

    plot_combined_roc(
        pattern='roc_clinical_*.npz',
        title='Clinical AD ROC Curves using Disease Trajectory',
        out_png=os.path.join(out, 'summary_roc_plot_all_subclusters.png'),
    )
    plot_combined_roc(
        pattern='roc_cerad_*.npz',
        title='Pathological AD ROC Curves using Disease Trajectory (CERAD)',
        out_png=os.path.join(out, 'summary_roc_plot_CERAD_all_subclusters.png'),
    )

    print('\n--- Aggregation complete. ---')


if __name__ == '__main__':
    main()
