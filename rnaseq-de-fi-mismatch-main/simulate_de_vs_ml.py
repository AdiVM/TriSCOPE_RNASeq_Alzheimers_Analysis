#!/usr/bin/env python3
# ruff: noqa: E402
from __future__ import annotations

import argparse
import os
from pathlib import Path

RUNTIME_CACHE_ROOT = Path.cwd() / ".runtime_cache"
for cache_dir in (
    RUNTIME_CACHE_ROOT,
    RUNTIME_CACHE_ROOT / "mplconfig",
    RUNTIME_CACHE_ROOT / "xdg_cache",
):
    cache_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(RUNTIME_CACHE_ROOT / "mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(RUNTIME_CACHE_ROOT / "xdg_cache"))
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier, early_stopping
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split


COUNTS_PER_MILLION = 1_000_000.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Simulate count-based RNA-seq data from a negative-binomial model and "
            "show why adjusted DE and multivariate LightGBM feature importance can disagree."
        )
    )
    parser.add_argument("--seed", type=int, default=7, help="Random seed.")
    parser.add_argument(
        "--n-per-stratum",
        type=int,
        default=80,
        help=(
            "Samples in each of the four sex x disease strata: "
            "male-control, female-control, male-disease, female-disease."
        ),
    )
    parser.add_argument(
        "--n-noise-genes",
        type=int,
        default=250,
        help="Number of pure noise genes to add to the simulation.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=20,
        help="Number of top LightGBM genes to summarize.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("simulation_output"),
        help="Directory where CSVs, plots, and the text summary will be written.",
    )
    return parser.parse_args()


def prepare_runtime(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / ".mplconfig").mkdir(parents=True, exist_ok=True)


def sample_negative_binomial(
    mean_counts: np.ndarray, dispersion_alpha: float, rng: np.random.Generator
) -> np.ndarray:
    gamma_shape = 1.0 / dispersion_alpha
    gamma_scale = mean_counts * dispersion_alpha
    latent_rate = rng.gamma(shape=gamma_shape, scale=gamma_scale)
    return rng.poisson(latent_rate)


def counts_to_log_cpm(counts_df: pd.DataFrame, metadata_df: pd.DataFrame) -> pd.DataFrame:
    library_size = metadata_df["library_size"].to_numpy(dtype=float)
    cpm = counts_df.to_numpy(dtype=float) / library_size[:, None] * COUNTS_PER_MILLION
    log_cpm = np.log2(cpm + 1.0)
    return pd.DataFrame(log_cpm, columns=counts_df.columns, index=counts_df.index)


def simulate_dataset(
    seed: int,
    n_per_stratum: int,
    n_noise_genes: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)

    n = 4 * n_per_stratum
    sex = np.array(
        [1] * n_per_stratum
        + [0] * n_per_stratum
        + [1] * n_per_stratum
        + [0] * n_per_stratum
    )
    disease = np.array([0] * (2 * n_per_stratum) + [1] * (2 * n_per_stratum))
    age = rng.normal(loc=60.0, scale=9.0, size=n)
    age_z = (age - age.mean()) / age.std(ddof=0)
    sex_pm = 2 * sex - 1

    hidden_subtype = np.zeros(n, dtype=int)
    disease_idx = np.flatnonzero(disease == 1)
    hidden_subtype[disease_idx] = rng.choice([-1, 1], size=len(disease_idx))

    log_size_factor = rng.normal(loc=0.0, scale=0.25, size=n)
    size_factor = np.exp(log_size_factor)
    library_size = np.round(15_000_000 * size_factor).astype(int)

    n_redundant_de = 25
    n_sex_interaction = 5
    n_heterogeneous = 5
    n_threshold_effect = 5
    n_distribution_shift = 5
    n_age_only = 10
    total_genes = (
        n_redundant_de
        + n_sex_interaction
        + n_heterogeneous
        + n_threshold_effect
        + n_distribution_shift
        + n_age_only
        + n_noise_genes
    )

    counts = np.zeros((n, total_genes), dtype=int)
    categories: list[str] = []
    base_mean_counts: list[float] = []
    dispersion_alpha: list[float] = []
    shared_de_module = rng.normal(loc=0.0, scale=0.35, size=n)
    col = 0

    for _ in range(n_redundant_de):
        base = float(rng.lognormal(mean=np.log(80.0), sigma=0.35))
        alpha = float(0.05 + 6.0 / base)
        loading = float(rng.normal(loc=0.95, scale=0.07))
        mean_counts = base * size_factor * np.exp(
            0.42 * disease + 0.15 * age_z + 0.05 * sex_pm + loading * shared_de_module
        )
        counts[:, col] = sample_negative_binomial(mean_counts, alpha, rng)
        categories.append("redundant_de")
        base_mean_counts.append(base)
        dispersion_alpha.append(alpha)
        col += 1

    for _ in range(n_sex_interaction):
        base = float(rng.lognormal(mean=np.log(60.0), sigma=0.30))
        alpha = float(0.06 + 6.5 / base)
        mean_factor = 1.0 + 0.25 * disease * sex_pm
        mean_counts = base * size_factor * mean_factor * np.exp(0.06 * age_z)
        counts[:, col] = sample_negative_binomial(mean_counts, alpha, rng)
        categories.append("sex_interaction")
        base_mean_counts.append(base)
        dispersion_alpha.append(alpha)
        col += 1

    for _ in range(n_heterogeneous):
        base = float(rng.lognormal(mean=np.log(60.0), sigma=0.30))
        alpha = float(0.06 + 6.5 / base)
        mean_factor = 1.0 + 0.70 * disease * hidden_subtype
        mean_counts = base * size_factor * mean_factor
        counts[:, col] = sample_negative_binomial(mean_counts, alpha, rng)
        categories.append("heterogeneous_disease")
        base_mean_counts.append(base)
        dispersion_alpha.append(alpha)
        col += 1

    for _ in range(n_threshold_effect):
        base = float(rng.lognormal(mean=np.log(55.0), sigma=0.30))
        alpha = float(0.08 + 7.5 / base)
        high_tail_gate = rng.binomial(1, 0.20, size=n)
        disease_factor = np.where(high_tail_gate == 1, 3.0, 0.5)
        mean_factor = np.where(disease == 1, disease_factor, 1.0)
        mean_counts = base * size_factor * mean_factor
        counts[:, col] = sample_negative_binomial(mean_counts, alpha, rng)
        categories.append("threshold_effect")
        base_mean_counts.append(base)
        dispersion_alpha.append(alpha)
        col += 1

    for _ in range(n_distribution_shift):
        base = float(rng.lognormal(mean=np.log(55.0), sigma=0.30))
        alpha = float(0.08 + 7.5 / base)
        signs = rng.choice([-1.0, 1.0], size=n)
        control_factor = 1.0 + 0.25 * signs
        disease_factor = 1.0 + 0.75 * signs
        mean_factor = np.where(disease == 1, disease_factor, control_factor)
        mean_counts = base * size_factor * mean_factor
        counts[:, col] = sample_negative_binomial(mean_counts, alpha, rng)
        categories.append("distribution_shift")
        base_mean_counts.append(base)
        dispersion_alpha.append(alpha)
        col += 1

    for _ in range(n_age_only):
        base = float(rng.lognormal(mean=np.log(55.0), sigma=0.35))
        alpha = float(0.08 + 8.0 / base)
        mean_counts = base * size_factor * np.exp(0.55 * age_z)
        counts[:, col] = sample_negative_binomial(mean_counts, alpha, rng)
        categories.append("age_only")
        base_mean_counts.append(base)
        dispersion_alpha.append(alpha)
        col += 1

    for _ in range(n_noise_genes):
        base = float(rng.lognormal(mean=np.log(50.0), sigma=0.55))
        alpha = float(0.10 + 9.0 / base)
        mean_counts = base * size_factor
        counts[:, col] = sample_negative_binomial(mean_counts, alpha, rng)
        categories.append("noise")
        base_mean_counts.append(base)
        dispersion_alpha.append(alpha)
        col += 1

    genes = [f"G{idx:03d}" for idx in range(total_genes)]
    sample_ids = [f"S{idx:03d}" for idx in range(n)]
    counts_df = pd.DataFrame(counts, index=sample_ids, columns=genes)
    metadata_df = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "disease": disease,
            "disease_label": np.where(disease == 1, "disease", "control"),
            "sex": sex,
            "sex_label": np.where(sex == 1, "male", "female"),
            "age": age,
            "age_z": age_z,
            "hidden_subtype": hidden_subtype,
            "subtype_label": np.where(
                disease == 0,
                "control",
                np.where(hidden_subtype == 1, "disease_subtype_A", "disease_subtype_B"),
            ),
            "size_factor": size_factor,
            "log_size_factor": log_size_factor,
            "library_size": library_size,
        },
        index=sample_ids,
    )
    log_cpm_df = counts_to_log_cpm(counts_df, metadata_df)
    gene_info_df = pd.DataFrame(
        {
            "gene": genes,
            "category": categories,
            "base_mean_count": base_mean_counts,
            "dispersion_alpha": dispersion_alpha,
        }
    )
    return counts_df, log_cpm_df, metadata_df, gene_info_df


def make_pydeseq_metadata(metadata_df: pd.DataFrame) -> pd.DataFrame:
    return metadata_df[["age", "sex_label", "disease_label"]].copy()


def build_interaction_contrast(interaction_dds: DeseqDataSet) -> np.ndarray:
    lfc_columns = interaction_dds.varm["LFC"].columns.tolist()
    interaction_column = next(
        column
        for column in lfc_columns
        if "sex_label" in column and "disease_label" in column
    )
    contrast = np.zeros(len(lfc_columns), dtype=float)
    contrast[lfc_columns.index(interaction_column)] = 1.0
    return contrast


def run_differential_expression(counts_df: pd.DataFrame, metadata_df: pd.DataFrame) -> pd.DataFrame:
    pydeseq_metadata = make_pydeseq_metadata(metadata_df)

    additive_dds = DeseqDataSet(
        counts=counts_df,
        metadata=pydeseq_metadata,
        design="~ age + sex_label + disease_label",
        refit_cooks=False,
        quiet=True,
        n_cpus=1,
    )
    additive_dds.deseq2()
    additive_stats = DeseqStats(
        additive_dds,
        contrast=["disease_label", "disease", "control"],
        quiet=True,
        n_cpus=1,
    )
    additive_stats.summary()
    additive_results = additive_stats.results_df.reset_index().rename(
        columns={
            "index": "gene",
            "baseMean": "deseq2_base_mean",
            "log2FoldChange": "de_log2_fc",
            "pvalue": "de_pvalue",
            "padj": "fdr",
        }
    )

    interaction_dds = DeseqDataSet(
        counts=counts_df,
        metadata=pydeseq_metadata,
        design="~ age + sex_label + disease_label + sex_label:disease_label",
        refit_cooks=False,
        quiet=True,
        n_cpus=1,
    )
    interaction_dds.deseq2()
    interaction_stats = DeseqStats(
        interaction_dds,
        contrast=build_interaction_contrast(interaction_dds),
        quiet=True,
        n_cpus=1,
    )
    interaction_stats.summary()
    interaction_results = interaction_stats.results_df.reset_index().rename(
        columns={
            "index": "gene",
            "log2FoldChange": "disease_sex_log2_fc",
            "pvalue": "disease_sex_pvalue",
            "padj": "disease_sex_fdr",
        }
    )

    de_df = additive_results[
        ["gene", "deseq2_base_mean", "de_log2_fc", "de_pvalue", "fdr"]
    ].merge(
        interaction_results[
            ["gene", "disease_sex_log2_fc", "disease_sex_pvalue", "disease_sex_fdr"]
        ],
        on="gene",
        how="left",
    )
    mean_counts = counts_df.mean(axis=0)
    missing_genes = pd.Index(de_df["gene"]).difference(mean_counts.index)
    if not missing_genes.empty:
        missing_preview = ", ".join(map(str, missing_genes[:10]))
        raise ValueError(
            "Differential expression results contain genes that are not present in "
            f"counts_df columns: {missing_preview}"
            + ("..." if len(missing_genes) > 10 else "")
        )
    de_df["mean_count"] = mean_counts.reindex(de_df["gene"]).to_numpy()
    de_df["de_significant"] = de_df["fdr"] < 0.05
    de_df["neg_log10_fdr"] = -np.log10(np.clip(de_df["fdr"], 1e-300, None))
    return de_df


def fit_lightgbm(
    log_cpm_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, float, int]:
    features = log_cpm_df.copy()
    features["age"] = metadata_df["age"]
    features["sex"] = metadata_df["sex"]
    target = metadata_df["disease"]

    x_train_full, x_test, y_train_full, y_test = train_test_split(
        features,
        target,
        test_size=0.35,
        random_state=seed,
        stratify=target,
    )
    x_train, x_valid, y_train, y_valid = train_test_split(
        x_train_full,
        y_train_full,
        test_size=0.25,
        random_state=seed,
        stratify=y_train_full,
    )

    model = LGBMClassifier(
        n_estimators=500,
        learning_rate=0.04,
        num_leaves=31,
        max_depth=5,
        min_child_samples=12,
        subsample=0.80,
        colsample_bytree=0.70,
        reg_alpha=0.20,
        reg_lambda=1.20,
        min_gain_to_split=0.02,
        importance_type="gain",
        n_jobs=1,
        random_state=seed,
        verbosity=-1,
    )
    model.fit(
        x_train,
        y_train,
        eval_set=[(x_valid, y_valid)],
        eval_metric="auc",
        callbacks=[early_stopping(stopping_rounds=30, verbose=False)],
    )

    test_prob = model.predict_proba(x_test)[:, 1]
    auc = roc_auc_score(y_test, test_prob)
    feature_importance_df = pd.DataFrame(
        {
            "feature": model.feature_name_,
            "lightgbm_gain": model.feature_importances_.astype(float),
        }
    ).sort_values("lightgbm_gain", ascending=False, ignore_index=True)

    gene_gain_df = feature_importance_df[
        feature_importance_df["feature"].isin(log_cpm_df.columns)
    ].rename(columns={"feature": "gene"})
    gene_gain_df = gene_gain_df.reset_index(drop=True)
    gene_gain_df["gain_rank"] = np.arange(1, len(gene_gain_df) + 1)
    gene_gain_df["nonzero_gain"] = gene_gain_df["lightgbm_gain"] > 0
    return gene_gain_df, feature_importance_df, auc, int(
        model.best_iteration_ if model.best_iteration_ is not None else model.n_estimators
    )


def make_category_summary(results_df: pd.DataFrame, top_k: int) -> pd.DataFrame:
    top_genes = set(results_df.sort_values("lightgbm_gain", ascending=False).head(top_k)["gene"])
    summary_df = (
        results_df.groupby("category", as_index=False)
        .agg(
            n_genes=("gene", "size"),
            de_significant=("de_significant", "sum"),
            nonzero_gain=("nonzero_gain", "sum"),
            median_fdr=("fdr", "median"),
            median_gain=("lightgbm_gain", "median"),
            median_abs_log2_fc=("de_log2_fc", lambda s: float(np.median(np.abs(s)))),
            median_disease_sex_pvalue=("disease_sex_pvalue", "median"),
            median_mean_count=("mean_count", "median"),
            median_dispersion_alpha=("dispersion_alpha", "median"),
        )
        .sort_values("category", ignore_index=True)
    )
    summary_df[f"top_{top_k}_gain"] = summary_df["category"].map(
        results_df.groupby("category")["gene"]
        .apply(lambda genes: sum(gene in top_genes for gene in genes))
        .to_dict()
    )
    summary_df["de_significant_fraction"] = summary_df["de_significant"] / summary_df["n_genes"]
    summary_df["nonzero_gain_fraction"] = summary_df["nonzero_gain"] / summary_df["n_genes"]
    summary_df[f"top_{top_k}_gain_fraction"] = summary_df[f"top_{top_k}_gain"] / summary_df["n_genes"]
    return summary_df


def format_pvalue(value: float) -> str:
    if value < 1e-3:
        return f"{value:.2e}"
    return f"{value:.3f}"


def build_text_summary(
    results_df: pd.DataFrame,
    category_summary_df: pd.DataFrame,
    auc: float,
    best_iteration: int,
    top_k: int,
) -> str:
    summary_lookup = category_summary_df.set_index("category")
    redundant = summary_lookup.loc["redundant_de"]
    sex_interaction = summary_lookup.loc["sex_interaction"]
    heterogeneous = summary_lookup.loc["heterogeneous_disease"]
    threshold = summary_lookup.loc["threshold_effect"]
    distribution = summary_lookup.loc["distribution_shift"]

    top_gain = results_df.sort_values("lightgbm_gain", ascending=False).head(top_k)
    top_de = results_df.sort_values("fdr").head(top_k)

    lines = [
        "Simulation summary",
        "==================",
        "Count model: negative-binomial counts with sample-specific library sizes",
        (
            "DE model: PyDESeq2 with additive design "
            "~ age + sex_label + disease_label, plus a second PyDESeq2 fit for the "
            "sex-by-disease interaction contrast"
        ),
        "ML model: LightGBM on log2-CPM for all genes + age + sex",
        f"LightGBM test ROC AUC: {auc:.3f}",
        f"Best boosting iteration after early stopping: {best_iteration}",
        "",
        "Why DE genes can have low tree importance",
        "-----------------------------------------",
        (
            f"{int(redundant['de_significant'])}/{int(redundant['n_genes'])} correlated "
            f"'redundant_de' genes are FDR-significant in the PyDESeq2 additive model, but only "
            f"{int(redundant['nonzero_gain'])}/{int(redundant['n_genes'])} receive non-zero "
            f"LightGBM gain and only {int(redundant[f'top_{top_k}_gain'])}/"
            f"{int(redundant['n_genes'])} land in the top {top_k} genes by gain."
        ),
        (
            "These genes all carry overlapping disease information, so the tree ensemble can "
            "split on a few representatives and leave the rest with little or zero marginal gain."
        ),
        "",
        "Why high-importance genes can fail standard DE",
        "----------------------------------------------",
        (
            f"Sex-interaction genes: {int(sex_interaction['de_significant'])}/"
            f"{int(sex_interaction['n_genes'])} are FDR-significant in the additive PyDESeq2 model, "
            f"but {int(sex_interaction['nonzero_gain'])}/"
            f"{int(sex_interaction['n_genes'])} have non-zero gain and "
            f"{int(sex_interaction[f'top_{top_k}_gain'])}/"
            f"{int(sex_interaction['n_genes'])} appear in the top {top_k} by gain. "
            f"The median disease:sex interaction p-value in the interaction PyDESeq2 model is "
            f"{format_pvalue(float(sex_interaction['median_disease_sex_pvalue']))}."
        ),
        (
            f"Heterogeneous disease genes: {int(heterogeneous['de_significant'])}/"
            f"{int(heterogeneous['n_genes'])} are FDR-significant, but "
            f"{int(heterogeneous['nonzero_gain'])}/{int(heterogeneous['n_genes'])} have non-zero "
            f"gain and {int(heterogeneous[f'top_{top_k}_gain'])}/"
            f"{int(heterogeneous['n_genes'])} appear in the top {top_k}. "
            "These genes were simulated so the disease group contains high and low subtypes "
            "whose mean counts cancel on the raw-count scale, even though the pattern remains predictive."
        ),
        (
            f"Threshold-effect genes: {int(threshold['de_significant'])}/"
            f"{int(threshold['n_genes'])} are FDR-significant, but "
            f"{int(threshold['nonzero_gain'])}/{int(threshold['n_genes'])} have non-zero gain "
            f"and {int(threshold[f'top_{top_k}_gain'])}/{int(threshold['n_genes'])} appear in "
            f"the top {top_k}. Disease samples have a rare very-high expression tail balanced "
            "by more common mildly low values, so a tree can learn a high-value rule even when "
            "the average effect is weak."
        ),
        (
            f"Distribution-shift genes: {int(distribution['de_significant'])}/"
            f"{int(distribution['n_genes'])} are FDR-significant, but "
            f"{int(distribution['nonzero_gain'])}/{int(distribution['n_genes'])} have non-zero gain "
            f"and {int(distribution[f'top_{top_k}_gain'])}/{int(distribution['n_genes'])} appear in "
            f"the top {top_k}. Control samples are relatively concentrated near the middle, while "
            "disease samples are more dispersed into the tails with the same mean, which helps "
            "tree splits more than a mean-shift test."
        ),
        "",
        "Interpretation",
        "--------------",
        (
            "The disagreement is expected because adjusted DE asks whether one gene has a disease "
            "main effect after accounting for covariates, one gene at a time, while LightGBM asks "
            "whether a transformed predictor is useful inside a joint nonlinear classifier that can "
            "exploit redundancy, thresholds, and interactions."
        ),
        "",
        f"Top genes by smallest FDR (top {top_k})",
        "---------------------------------------",
    ]
    for row in top_de.itertuples(index=False):
        lines.append(
            f"{row.gene}: {row.category}, FDR={row.fdr:.3e}, gain={row.lightgbm_gain:.2f}"
        )

    lines.extend(["", f"Top genes by LightGBM gain (top {top_k})", "------------------------------------"])
    for row in top_gain.itertuples(index=False):
        lines.append(
            f"{row.gene}: {row.category}, gain={row.lightgbm_gain:.2f}, FDR={row.fdr:.3e}"
        )

    return "\n".join(lines)


def plot_de_vs_gain(
    results_df: pd.DataFrame, example_genes: dict[str, str], output_path: Path
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    palette = {
        "redundant_de": "#1b9e77",
        "sex_interaction": "#d95f02",
        "heterogeneous_disease": "#7570b3",
        "threshold_effect": "#e6ab02",
        "distribution_shift": "#1f78b4",
        "age_only": "#66a61e",
        "noise": "#bdbdbd",
    }

    plot_df = results_df.copy()
    plot_df["log1p_gain"] = np.log1p(plot_df["lightgbm_gain"])

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.scatterplot(
        data=plot_df,
        x="log1p_gain",
        y="neg_log10_fdr",
        hue="category",
        palette=palette,
        s=60,
        alpha=0.85,
        linewidth=0,
        ax=ax,
    )
    ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("log(1 + LightGBM gain importance)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("Count-based DE significance and LightGBM importance can diverge")

    for gene in example_genes.values():
        row = plot_df.loc[plot_df["gene"] == gene].iloc[0]
        ax.text(row["log1p_gain"] + 0.03, row["neg_log10_fdr"] + 0.05, gene, fontsize=10)

    ax.legend(title="Simulated gene type", frameon=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_example_patterns(
    log_cpm_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    example_genes: dict[str, str],
    output_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(2, 3, figsize=(18, 9.5))
    axes = axes.ravel()

    redundant_gene = example_genes["redundant_de"]
    redundant_df = metadata_df.copy()
    redundant_df["expression"] = log_cpm_df.loc[redundant_df.index, redundant_gene]
    sns.boxplot(
        data=redundant_df,
        x="disease_label",
        y="expression",
        color="#9ad8c2",
        ax=axes[0],
    )
    sns.stripplot(
        data=redundant_df,
        x="disease_label",
        y="expression",
        color="black",
        alpha=0.22,
        size=2.3,
        ax=axes[0],
    )
    axes[0].set_title(f"{redundant_gene}: correlated DE block")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Simulated log2-CPM")

    sex_gene = example_genes["sex_interaction"]
    sex_df = metadata_df.copy()
    sex_df["expression"] = log_cpm_df.loc[sex_df.index, sex_gene]
    sns.boxplot(
        data=sex_df,
        x="sex_label",
        y="expression",
        hue="disease_label",
        palette={"control": "#9ecae1", "disease": "#fc9272"},
        ax=axes[1],
    )
    axes[1].set_title(f"{sex_gene}: disease effect differs by sex")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    axes[1].legend(title="", loc="upper left", frameon=True)

    hetero_gene = example_genes["heterogeneous_disease"]
    hetero_df = metadata_df.copy()
    hetero_df["expression"] = log_cpm_df.loc[hetero_df.index, hetero_gene]
    sns.boxplot(
        data=hetero_df,
        x="subtype_label",
        y="expression",
        hue="subtype_label",
        dodge=False,
        palette={
            "control": "#d9d9d9",
            "disease_subtype_A": "#bcbddc",
            "disease_subtype_B": "#9e9ac8",
        },
        ax=axes[2],
    )
    sns.stripplot(
        data=hetero_df,
        x="subtype_label",
        y="expression",
        color="black",
        alpha=0.18,
        size=1.8,
        ax=axes[2],
    )
    axes[2].set_title(f"{hetero_gene}: hidden disease subtypes cancel in mean count")
    axes[2].set_xlabel("")
    axes[2].set_ylabel("")
    axes[2].tick_params(axis="x", rotation=12)
    if axes[2].legend_ is not None:
        axes[2].legend_.remove()

    threshold_gene = example_genes["threshold_effect"]
    threshold_df = metadata_df.copy()
    threshold_df["expression"] = log_cpm_df.loc[threshold_df.index, threshold_gene]
    sns.violinplot(
        data=threshold_df,
        x="disease_label",
        y="expression",
        color="#f0d37a",
        inner="quartile",
        cut=0,
        ax=axes[3],
    )
    sns.stripplot(
        data=threshold_df,
        x="disease_label",
        y="expression",
        color="black",
        alpha=0.16,
        size=1.8,
        ax=axes[3],
    )
    axes[3].set_title(f"{threshold_gene}: rare high disease tail")
    axes[3].set_xlabel("")
    axes[3].set_ylabel("Simulated log2-CPM")

    distribution_gene = example_genes["distribution_shift"]
    distribution_df = metadata_df.copy()
    distribution_df["expression"] = log_cpm_df.loc[distribution_df.index, distribution_gene]
    sns.violinplot(
        data=distribution_df,
        x="disease_label",
        y="expression",
        color="#9ecae1",
        inner="quartile",
        cut=0,
        ax=axes[4],
    )
    sns.stripplot(
        data=distribution_df,
        x="disease_label",
        y="expression",
        color="black",
        alpha=0.16,
        size=1.8,
        ax=axes[4],
    )
    axes[4].set_title(f"{distribution_gene}: wider disease distribution, same center")
    axes[4].set_xlabel("")
    axes[4].set_ylabel("")

    axes[5].axis("off")
    axes[5].text(
        0.02,
        0.92,
        "Mechanisms shown",
        fontsize=13,
        fontweight="bold",
        ha="left",
        va="top",
        transform=axes[5].transAxes,
    )
    axes[5].text(
        0.02,
        0.78,
        "1. Correlated redundant DE genes",
        fontsize=11,
        ha="left",
        va="top",
        transform=axes[5].transAxes,
    )
    axes[5].text(
        0.02,
        0.64,
        "2. Disease-by-sex interaction",
        fontsize=11,
        ha="left",
        va="top",
        transform=axes[5].transAxes,
    )
    axes[5].text(
        0.02,
        0.50,
        "3. Hidden disease subtypes",
        fontsize=11,
        ha="left",
        va="top",
        transform=axes[5].transAxes,
    )
    axes[5].text(
        0.02,
        0.36,
        "4. Nonlinear threshold effect",
        fontsize=11,
        ha="left",
        va="top",
        transform=axes[5].transAxes,
    )
    axes[5].text(
        0.02,
        0.22,
        "5. Distributional change without a mean shift",
        fontsize=11,
        ha="left",
        va="top",
        transform=axes[5].transAxes,
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_two_tails_logic(
    log_cpm_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    heterogeneous_gene: str,
    output_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    plot_df = metadata_df.copy()
    plot_df["expression"] = log_cpm_df.loc[plot_df.index, heterogeneous_gene].to_numpy()

    control_expression = plot_df.loc[plot_df["disease_label"] == "control", "expression"].to_numpy()
    low_threshold, high_threshold = np.quantile(control_expression, [0.20, 0.80])
    x_min = plot_df["expression"].min() - 0.35
    x_max = plot_df["expression"].max() + 0.35

    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw={"width_ratios": [1.15, 1]})

    for ax in axes:
        ax.axvspan(x_min, low_threshold, color="#c6dbef", alpha=0.35, linewidth=0)
        ax.axvspan(low_threshold, high_threshold, color="#d9d9d9", alpha=0.40, linewidth=0)
        ax.axvspan(high_threshold, x_max, color="#dadaeb", alpha=0.40, linewidth=0)
        ax.axvline(low_threshold, color="#4d4d4d", linestyle="--", linewidth=1)
        ax.axvline(high_threshold, color="#4d4d4d", linestyle="--", linewidth=1)
        ax.set_xlim(x_min, x_max)

    sns.kdeplot(
        data=plot_df.loc[plot_df["disease_label"] == "control"],
        x="expression",
        fill=True,
        alpha=0.25,
        linewidth=2,
        color="#7f7f7f",
        label="control",
        ax=axes[0],
    )
    sns.kdeplot(
        data=plot_df.loc[plot_df["disease_label"] == "disease"],
        x="expression",
        fill=True,
        alpha=0.25,
        linewidth=2,
        color="#756bb1",
        label="disease",
        ax=axes[0],
    )
    axes[0].set_title(f"{heterogeneous_gene}: disease is bimodal, control is central")
    axes[0].set_xlabel("Simulated log2-CPM")
    axes[0].set_ylabel("Density")
    axes[0].legend(title="", frameon=True, loc="upper left")
    ymax = axes[0].get_ylim()[1]
    axes[0].text((x_min + low_threshold) / 2, ymax * 0.95, "Disease tail", ha="center", va="top", fontsize=10)
    axes[0].text((low_threshold + high_threshold) / 2, ymax * 0.95, "Control-rich middle", ha="center", va="top", fontsize=10)
    axes[0].text((high_threshold + x_max) / 2, ymax * 0.95, "Disease tail", ha="center", va="top", fontsize=10)

    subtype_order = ["disease_subtype_B", "control", "disease_subtype_A"]
    subtype_palette = {
        "control": "#6b6b6b",
        "disease_subtype_A": "#756bb1",
        "disease_subtype_B": "#9e9ac8",
    }
    subtype_rows = {label: idx for idx, label in enumerate(subtype_order)}
    rng = np.random.default_rng(123)
    scatter_df = plot_df.copy()
    scatter_df["y"] = scatter_df["subtype_label"].map(subtype_rows).astype(float)
    scatter_df["y"] += rng.uniform(-0.12, 0.12, size=len(scatter_df))

    sns.scatterplot(
        data=scatter_df,
        x="expression",
        y="y",
        hue="subtype_label",
        hue_order=subtype_order,
        palette=subtype_palette,
        s=28,
        alpha=0.80,
        edgecolor="none",
        ax=axes[1],
    )
    axes[1].set_title("One-gene logic: middle = control, both tails = disease")
    axes[1].set_xlabel("Simulated log2-CPM")
    axes[1].set_ylabel("")
    axes[1].set_yticks([0, 1, 2])
    axes[1].set_yticklabels(["Disease subtype B", "Control", "Disease subtype A"])
    axes[1].legend(title="", frameon=True, loc="lower right")
    axes[1].text((x_min + low_threshold) / 2, 2.38, "Low tail", ha="center", va="bottom", fontsize=10)
    axes[1].text((low_threshold + high_threshold) / 2, 2.38, "Middle", ha="center", va="bottom", fontsize=10)
    axes[1].text((high_threshold + x_max) / 2, 2.38, "High tail", ha="center", va="bottom", fontsize=10)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    prepare_runtime(args.output_dir)

    counts_df, log_cpm_df, metadata_df, gene_info_df = simulate_dataset(
        seed=args.seed,
        n_per_stratum=args.n_per_stratum,
        n_noise_genes=args.n_noise_genes,
    )
    de_df = run_differential_expression(counts_df, metadata_df)
    gain_df, all_feature_importance_df, auc, best_iteration = fit_lightgbm(
        log_cpm_df,
        metadata_df,
        seed=args.seed,
    )

    results_df = (
        gene_info_df.merge(de_df, on="gene", how="left")
        .merge(gain_df, on="gene", how="left")
        .sort_values(["category", "gene"], ignore_index=True)
    )
    category_summary_df = make_category_summary(results_df, top_k=args.top_k)

    example_genes = {
        "redundant_de": results_df.loc[
            results_df["category"] == "redundant_de"
        ].sort_values("fdr", ascending=True)["gene"].iloc[0],
        "sex_interaction": results_df.loc[
            results_df["category"] == "sex_interaction"
        ].sort_values("lightgbm_gain", ascending=False)["gene"].iloc[0],
        "heterogeneous_disease": results_df.loc[
            results_df["category"] == "heterogeneous_disease"
        ].sort_values("lightgbm_gain", ascending=False)["gene"].iloc[0],
        "threshold_effect": results_df.loc[
            results_df["category"] == "threshold_effect"
        ].sort_values("lightgbm_gain", ascending=False)["gene"].iloc[0],
        "distribution_shift": results_df.loc[
            results_df["category"] == "distribution_shift"
        ].sort_values("lightgbm_gain", ascending=False)["gene"].iloc[0],
    }

    summary_text = build_text_summary(
        results_df=results_df,
        category_summary_df=category_summary_df,
        auc=auc,
        best_iteration=best_iteration,
        top_k=args.top_k,
    )

    metadata_df.to_csv(args.output_dir / "sample_metadata.csv", index=False)
    counts_df.to_csv(args.output_dir / "counts_matrix.csv", index_label="sample_id")
    log_cpm_df.to_csv(args.output_dir / "log_cpm_matrix.csv", index_label="sample_id")
    results_df.to_csv(args.output_dir / "gene_level_results.csv", index=False)
    category_summary_df.to_csv(args.output_dir / "category_summary.csv", index=False)
    all_feature_importance_df.to_csv(args.output_dir / "all_predictor_importance.csv", index=False)
    (args.output_dir / "summary.txt").write_text(summary_text + "\n", encoding="utf-8")

    plot_de_vs_gain(
        results_df=results_df,
        example_genes=example_genes,
        output_path=args.output_dir / "de_vs_gain_scatter.png",
    )
    plot_example_patterns(
        log_cpm_df=log_cpm_df,
        metadata_df=metadata_df,
        example_genes=example_genes,
        output_path=args.output_dir / "mechanism_examples.png",
    )
    plot_two_tails_logic(
        log_cpm_df=log_cpm_df,
        metadata_df=metadata_df,
        heterogeneous_gene=example_genes["heterogeneous_disease"],
        output_path=args.output_dir / "two_tails_logic.png",
    )

    print(summary_text)
    print("")
    print(f"Wrote outputs to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
