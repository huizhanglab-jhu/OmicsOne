from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_feature_correlation(
    feature_corr: pd.DataFrame,
    *,
    omics1: str,
    omics2: str,
    output_path: str | Path,
    dpi: int = 180,
) -> None:
    values = feature_corr["Gene Correlation"].dropna().astype(float)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    sns.histplot(values, bins=40, kde=True, color="#3b7ddd", ax=ax)
    ax.axvline(values.median(), color="#c44536", label=f"median={values.median():.3f}")
    ax.set_xlim(-1, 1)
    ax.set_title(f"Matched-feature correlation\nomics1={omics1}; omics2={omics2}")
    ax.set_xlabel(f"Spearman rho: {omics1} vs {omics2}")
    ax.legend(frameon=False)
    _save(fig, output_path, dpi)


def plot_sample_correlation(
    sample_corr: pd.DataFrame,
    *,
    omics1: str,
    omics2: str,
    output_path: str | Path,
    dpi: int = 180,
) -> None:
    values = sample_corr["Corr"].dropna().astype(float).sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(9, 4))
    colors = np.where(values.ge(0), "#3b7ddd", "#c44536")
    ax.bar(np.arange(len(values)), values, color=colors, width=0.85)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_title(f"Sample-wise correlation\nomics1={omics1}; omics2={omics2}")
    ax.set_xlabel("Samples sorted by rho")
    ax.set_ylabel(f"Spearman rho: {omics1} vs {omics2}")
    ax.set_xticks([])
    _save(fig, output_path, dpi)


def plot_gene_sample_correlation(
    gene_sample_corr: pd.DataFrame,
    *,
    omics1: str,
    omics2: str,
    output_path: str | Path,
    dpi: int = 180,
) -> None:
    data = gene_sample_corr[["Corr_omics1", "Corr_omics2"]].dropna().astype(float)
    fig, ax = plt.subplots(figsize=(5.5, 5))
    if len(data) > 1000:
        density = ax.hexbin(data["Corr_omics1"], data["Corr_omics2"], gridsize=45, cmap="viridis", mincnt=1)
        fig.colorbar(density, ax=ax, label="features")
    else:
        ax.scatter(data["Corr_omics1"], data["Corr_omics2"], s=14, alpha=0.55, color="#3b7ddd", linewidths=0)
    ax.axhline(0, color="0.4", linewidth=0.8, linestyle="--")
    ax.axvline(0, color="0.4", linewidth=0.8, linestyle="--")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_title(f"Gene-sample correlation\nomics1={omics1}; omics2={omics2}")
    ax.set_xlabel(f"Corr(sample rho, {omics1} feature) [omics1]")
    ax.set_ylabel(f"Corr(sample rho, {omics2} feature) [omics2]")
    _save(fig, output_path, dpi)


def plot_top_feature_correlations(
    feature_corr: pd.DataFrame,
    *,
    omics1: str,
    omics2: str,
    output_path: str | Path,
    count: int = 12,
    dpi: int = 180,
) -> None:
    values = feature_corr[["Gene Correlation"]].dropna().astype(float)
    selected = pd.concat(
        [values.nlargest(count, "Gene Correlation"), values.nsmallest(count, "Gene Correlation")]
    )
    selected = selected.loc[~selected.index.duplicated()].sort_values("Gene Correlation")
    labels = [str(value)[:32] for value in selected.index]
    fig, ax = plt.subplots(figsize=(7.5, max(4.5, 0.24 * len(selected) + 1.5)))
    colors = np.where(selected["Gene Correlation"].ge(0), "#3b7ddd", "#c44536")
    ax.barh(labels, selected["Gene Correlation"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlim(-1, 1)
    ax.set_title(f"Strongest matched-feature correlations\nomics1={omics1}; omics2={omics2}")
    ax.set_xlabel(f"Spearman rho: {omics1} vs {omics2}")
    _save(fig, output_path, dpi)


def plot_cluster_ari(
    ari: pd.DataFrame,
    *,
    omics1: str,
    omics2: str,
    output_path: str | Path,
    dpi: int = 180,
) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(ari, vmin=-1, vmax=1, cmap="RdBu_r", annot=True, fmt=".2f", ax=ax)
    ax.set_title(f"Cluster agreement (ARI)\nomics1={omics1}; omics2={omics2}")
    _save(fig, output_path, dpi)


def plot_pair_overview(
    feature_corr: pd.DataFrame,
    sample_corr: pd.DataFrame,
    gene_sample_corr: pd.DataFrame,
    ari: pd.DataFrame | None,
    *,
    omics1: str,
    omics2: str,
    output_path: str | Path,
    dpi: int = 180,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    fig.suptitle(f"OmicsX overview\nomics1={omics1}; omics2={omics2}", fontsize=16, y=0.995)

    feature_values = feature_corr["Gene Correlation"].dropna().astype(float)
    sns.histplot(feature_values, bins=40, kde=True, color="#3b7ddd", ax=axes[0, 0])
    axes[0, 0].axvline(feature_values.median(), color="#c44536", label=f"median={feature_values.median():.3f}")
    axes[0, 0].set_xlim(-1, 1)
    axes[0, 0].set_title("Matched-feature correlation")
    axes[0, 0].set_xlabel(f"Spearman rho: {omics1} vs {omics2}")
    axes[0, 0].legend(frameon=False)

    sample_values = sample_corr["Corr"].dropna().astype(float).sort_values(ascending=False)
    colors = np.where(sample_values.ge(0), "#3b7ddd", "#c44536")
    axes[0, 1].bar(np.arange(len(sample_values)), sample_values, color=colors, width=0.85)
    axes[0, 1].axhline(0, color="black", linewidth=0.8)
    axes[0, 1].set_title("Sample-wise correlation")
    axes[0, 1].set_xlabel("Samples sorted by rho")
    axes[0, 1].set_ylabel(f"Spearman rho: {omics1} vs {omics2}")
    axes[0, 1].set_xticks([])

    gene_sample = gene_sample_corr[["Corr_omics1", "Corr_omics2"]].dropna().astype(float)
    if len(gene_sample) > 1000:
        axes[1, 0].hexbin(
            gene_sample["Corr_omics1"], gene_sample["Corr_omics2"], gridsize=35, cmap="viridis", mincnt=1
        )
    else:
        axes[1, 0].scatter(
            gene_sample["Corr_omics1"], gene_sample["Corr_omics2"], s=12, alpha=0.55, color="#3b7ddd"
        )
    axes[1, 0].axhline(0, color="0.4", linewidth=0.8, linestyle="--")
    axes[1, 0].axvline(0, color="0.4", linewidth=0.8, linestyle="--")
    axes[1, 0].set_xlim(-1, 1)
    axes[1, 0].set_ylim(-1, 1)
    axes[1, 0].set_title("Gene-sample correlation")
    axes[1, 0].set_xlabel(f"{omics1} feature correlation [omics1]")
    axes[1, 0].set_ylabel(f"{omics2} feature correlation [omics2]")

    if ari is None or ari.empty:
        axes[1, 1].text(0.5, 0.5, "Cluster ARI unavailable", ha="center", va="center")
        axes[1, 1].axis("off")
    else:
        sns.heatmap(ari, vmin=-1, vmax=1, cmap="RdBu_r", annot=True, fmt=".2f", cbar=False, ax=axes[1, 1])
        axes[1, 1].set_title("Cluster agreement (ARI)")

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, output_path, dpi)


def plot_all_pair_medians(metrics: pd.DataFrame, *, output_dir: str | Path, dpi: int = 180) -> None:
    output = Path(output_dir)
    for column, title, filename, color in (
        ("gene_corr_median", "Median matched-feature Spearman rho", "all_pairs_gene_corr_median.png", "#3b7ddd"),
        ("sample_corr_median", "Median sample-wise Spearman rho", "all_pairs_sample_corr_median.png", "#2f9e44"),
    ):
        data = metrics.dropna(subset=[column]).sort_values(column)
        labels = [f"{row.pair}\nomics1={row.omics1}; omics2={row.omics2}" for row in data.itertuples()]
        fig, ax = plt.subplots(figsize=(10, max(5, 0.5 * len(data))))
        ax.barh(labels, data[column], color=color)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlabel(title)
        ax.set_title(f"OmicsX: {title} by pair")
        _save(fig, output / filename, dpi)


def _save(fig: plt.Figure, path: str | Path, dpi: int) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
