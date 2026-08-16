from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


def load_omics_matrix(path: str | Path) -> pd.DataFrame:
    """Read a feature-by-sample TSV without discarding partially observed rows."""
    matrix = pd.read_csv(path, sep="\t", header=0, index_col=0)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.loc[
        ~matrix.index.duplicated(keep="first"),
        ~matrix.columns.duplicated(keep="first"),
    ]
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    return matrix.dropna(axis=0, how="all").dropna(axis=1, how="all")


def align_matrices(omics1: pd.DataFrame, omics2: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return deterministically ordered common features and samples."""
    features = sorted(set(omics1.index.astype(str)) & set(omics2.index.astype(str)))
    samples = sorted(set(omics1.columns.astype(str)) & set(omics2.columns.astype(str)))
    return omics1.loc[features, samples].copy(), omics2.loc[features, samples].copy()


def feature_correlations(
    omics1: pd.DataFrame,
    omics2: pd.DataFrame,
    *,
    min_valid_pairs: int = 10,
) -> pd.DataFrame:
    """Calculate one Spearman correlation per matched feature."""
    _validate_aligned(omics1, omics2)
    rows: list[list[object]] = []
    for feature in omics1.index:
        x = omics1.loc[feature].to_numpy(dtype=float)
        y = omics2.loc[feature].to_numpy(dtype=float)
        result = _complete_spearman(x, y, min_valid_pairs)
        if result is None:
            continue
        rho, p_value, n_valid = result
        rows.append(
            [
                feature,
                rho,
                p_value,
                _coefficient_of_variation(x[np.isfinite(x) & np.isfinite(y)]),
                _coefficient_of_variation(y[np.isfinite(x) & np.isfinite(y)]),
                n_valid,
            ]
        )

    columns = ["Feature", "Gene Correlation", "P", "CV1", "CV2", "N"]
    result = pd.DataFrame(rows, columns=columns).set_index("Feature")
    if result.empty:
        result["BH adjusted P"] = pd.Series(dtype=float)
        return result
    result["BH adjusted P"] = multipletests(result["P"].to_numpy(dtype=float), method="fdr_bh")[1]
    return result.sort_values(["Gene Correlation", "P"], ascending=[False, True])


def sample_correlations(
    omics1: pd.DataFrame,
    omics2: pd.DataFrame,
    *,
    min_valid_pairs: int = 10,
) -> pd.DataFrame:
    """Calculate one Spearman correlation per matched sample across features."""
    _validate_aligned(omics1, omics2)
    rows: list[list[object]] = []
    for sample in omics1.columns:
        result = _complete_spearman(
            omics1[sample].to_numpy(dtype=float),
            omics2[sample].to_numpy(dtype=float),
            min_valid_pairs,
        )
        if result is None:
            continue
        rho, p_value, n_valid = result
        rows.append([sample, rho, p_value, n_valid])
    result = pd.DataFrame(rows, columns=["Sample", "Corr", "p-value", "N"]).set_index("Sample")
    if result.empty:
        return result
    return result.sort_values(["Corr", "p-value"], ascending=[False, True])


def gene_sample_correlations(
    sample_corr: pd.DataFrame,
    omics1: pd.DataFrame,
    omics2: pd.DataFrame,
    *,
    min_valid_pairs: int = 10,
) -> pd.DataFrame:
    """Relate feature abundance to sample-level cross-omics concordance."""
    _validate_aligned(omics1, omics2)
    if "Corr" not in sample_corr.columns:
        raise ValueError("sample_corr must contain a 'Corr' column")
    sample_rho = sample_corr["Corr"].reindex(omics1.columns).to_numpy(dtype=float)
    rows: list[list[object]] = []
    for feature in omics1.index:
        values1 = omics1.loc[feature].to_numpy(dtype=float)
        values2 = omics2.loc[feature].to_numpy(dtype=float)
        n1 = int((np.isfinite(sample_rho) & np.isfinite(values1)).sum())
        n2 = int((np.isfinite(sample_rho) & np.isfinite(values2)).sum())
        result1 = _complete_spearman(sample_rho, values1, min_valid_pairs)
        result2 = _complete_spearman(sample_rho, values2, min_valid_pairs)
        if result1 is None and result2 is None:
            continue
        rho1, p1, _ = result1 or (np.nan, np.nan, n1)
        rho2, p2, _ = result2 or (np.nan, np.nan, n2)
        rows.append([feature, rho1, p1, rho2, p2, n1, n2])
    result = pd.DataFrame(
        rows,
        columns=["Feature", "Corr_omics1", "P_omics1", "Corr_omics2", "P_omics2", "N_omics1", "N_omics2"],
    ).set_index("Feature")
    if result.empty:
        return result
    return result.sort_values(["Corr_omics1", "Corr_omics2"], ascending=[False, False], na_position="last")


def select_high_low_features(
    feature_corr: pd.DataFrame,
    *,
    fraction: float = 0.05,
    max_variable_features: int = 1000,
) -> tuple[list[str], list[str]]:
    """Select correlation tails after a symmetric variability filter."""
    required = {"Gene Correlation", "CV1", "CV2"}
    missing = required - set(feature_corr.columns)
    if missing:
        raise ValueError(f"feature_corr is missing columns: {', '.join(sorted(missing))}")
    if not 0 < fraction <= 0.5:
        raise ValueError("fraction must be greater than 0 and no more than 0.5")
    if max_variable_features < 1:
        raise ValueError("max_variable_features must be at least 1")
    if feature_corr.empty:
        return [], []

    working = feature_corr.copy()
    working["Variability"] = working[["CV1", "CV2"]].max(axis=1, skipna=True)
    working = working.dropna(subset=["Gene Correlation", "Variability"])
    if working.empty:
        return [], []
    variable = working.nlargest(min(max_variable_features, len(working)), "Variability")
    tail_count = max(1, int(round(len(working) * fraction)))
    high = set(working.nlargest(tail_count, "Gene Correlation").index) & set(variable.index)
    low = set(working.nsmallest(tail_count, "Gene Correlation").index) & set(variable.index)
    return sorted(map(str, high)), sorted(map(str, low))


def _validate_aligned(omics1: pd.DataFrame, omics2: pd.DataFrame) -> None:
    if not omics1.index.equals(omics2.index) or not omics1.columns.equals(omics2.columns):
        raise ValueError("omics1 and omics2 must have identical feature and sample order; call align_matrices first")


def _complete_spearman(
    x: Iterable[float],
    y: Iterable[float],
    min_valid_pairs: int,
) -> tuple[float, float, int] | None:
    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)
    valid = np.isfinite(x_array) & np.isfinite(y_array)
    n_valid = int(valid.sum())
    if n_valid < min_valid_pairs:
        return None
    x_valid = x_array[valid]
    y_valid = y_array[valid]
    if np.unique(x_valid).size < 2 or np.unique(y_valid).size < 2:
        return None
    rho, p_value = spearmanr(x_valid, y_valid)
    if not np.isfinite(rho) or not np.isfinite(p_value):
        return None
    return float(rho), float(p_value), n_valid


def _coefficient_of_variation(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        return np.nan
    mean = float(np.mean(values))
    if mean == 0 or not np.isfinite(mean):
        return np.nan
    return float(np.std(values, ddof=1) / abs(mean))
