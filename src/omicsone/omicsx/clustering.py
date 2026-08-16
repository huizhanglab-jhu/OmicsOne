from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage


def cluster_samples(
    matrix: pd.DataFrame,
    *,
    name: str,
    output_path: str | Path,
    n_clusters: int = 2,
    dpi: int = 180,
) -> pd.DataFrame:
    """Cluster samples with Ward linkage on row-standardized feature vectors."""
    if n_clusters < 2:
        raise ValueError("n_clusters must be at least 2")
    standardized = _standardize_for_clustering(matrix)
    if standardized.shape[0] < 2:
        raise ValueError("At least two variable features are required for sample clustering")
    if standardized.shape[1] < 2:
        raise ValueError("At least two samples are required for sample clustering")

    sample_linkage = linkage(standardized.T.to_numpy(), method="ward", metric="euclidean")
    raw_labels = fcluster(sample_linkage, t=n_clusters, criterion="maxclust")
    labels = pd.DataFrame(
        {"group": [f"C{value}" for value in raw_labels]},
        index=standardized.columns,
    )
    labels.index.name = "Sample"

    palette = dict(zip(sorted(labels["group"].unique()), sns.color_palette("Set2", labels["group"].nunique())))
    col_colors = labels["group"].map(palette)
    grid = sns.clustermap(
        standardized,
        col_linkage=sample_linkage,
        col_colors=col_colors,
        cmap="RdBu_r",
        vmin=-2,
        vmax=2,
        figsize=(6, 5),
        xticklabels=False,
        yticklabels=False,
    )
    grid.ax_heatmap.set_xlabel("Samples")
    grid.ax_heatmap.set_ylabel("Features")
    grid.fig.suptitle(name, y=1.02)
    grid.fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(grid.fig)
    return labels.sort_index()


def adjusted_rand_matrix(cluster_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Compare sample cluster assignments across OmicsX cluster views."""
    names = list(cluster_tables)
    result = pd.DataFrame(index=names, columns=names, dtype=float)
    for left in names:
        for right in names:
            joined = cluster_tables[left][["group"]].join(
                cluster_tables[right][["group"]],
                how="inner",
                lsuffix="_left",
                rsuffix="_right",
            )
            if joined.empty:
                result.loc[left, right] = np.nan
            else:
                result.loc[left, right] = _adjusted_rand_score(joined["group_left"], joined["group_right"])
    return result


def _standardize_for_clustering(matrix: pd.DataFrame) -> pd.DataFrame:
    numeric = matrix.apply(pd.to_numeric, errors="coerce").copy()
    row_medians = numeric.median(axis=1, skipna=True)
    numeric = numeric.T.fillna(row_medians).T.dropna(axis=0, how="any")
    means = numeric.mean(axis=1)
    standard_deviations = numeric.std(axis=1, ddof=1)
    keep = standard_deviations.gt(0) & np.isfinite(standard_deviations)
    return numeric.loc[keep].sub(means.loc[keep], axis=0).div(standard_deviations.loc[keep], axis=0)


def _adjusted_rand_score(labels1: pd.Series, labels2: pd.Series) -> float:
    contingency = pd.crosstab(labels1.to_numpy(), labels2.to_numpy()).to_numpy(dtype=np.int64)
    n_samples = int(contingency.sum())
    if n_samples < 2:
        return 1.0

    sum_cells = float(_comb2(contingency).sum())
    sum_rows = float(_comb2(contingency.sum(axis=1)).sum())
    sum_columns = float(_comb2(contingency.sum(axis=0)).sum())
    total_pairs = float(_comb2(np.array([n_samples], dtype=np.int64))[0])
    expected = sum_rows * sum_columns / total_pairs
    maximum = 0.5 * (sum_rows + sum_columns)
    denominator = maximum - expected
    if denominator == 0:
        return 1.0
    return (sum_cells - expected) / denominator


def _comb2(values: np.ndarray) -> np.ndarray:
    return values * (values - 1) // 2
