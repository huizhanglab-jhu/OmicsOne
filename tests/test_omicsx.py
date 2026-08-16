from __future__ import annotations

import configparser
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from omicsone.omicsx.analysis import (
    align_matrices,
    feature_correlations,
    gene_sample_correlations,
    load_omics_matrix,
    sample_correlations,
    select_high_low_features,
)
from omicsone.omicsx.clustering import adjusted_rand_matrix, cluster_samples
from omicsone.omicsx.pipeline import run_omicsx_config


def test_load_and_align_keep_partial_rows_and_sort_overlap(tmp_path: Path):
    input_path = tmp_path / "matrix.tsv"
    pd.DataFrame(
        {"s2": [1.0, np.nan, np.nan], "s1": [2.0, 3.0, np.nan]},
        index=["g2", "g1", "all_missing"],
    ).to_csv(input_path, sep="\t")

    loaded = load_omics_matrix(input_path)

    assert list(loaded.index) == ["g2", "g1"]
    assert loaded.loc["g1", "s2"] != loaded.loc["g1", "s2"]

    other = pd.DataFrame({"s1": [4.0, 5.0], "s2": [6.0, 7.0]}, index=["g1", "g2"])
    aligned1, aligned2 = align_matrices(loaded, other)
    assert list(aligned1.index) == ["g1", "g2"]
    assert list(aligned1.columns) == ["s1", "s2"]
    pd.testing.assert_index_equal(aligned1.index, aligned2.index)
    pd.testing.assert_index_equal(aligned1.columns, aligned2.columns)


def test_feature_correlations_use_pairwise_complete_values_and_fdr():
    omics1 = pd.DataFrame(
        {
            "s1": [1.0, 1.0, 4.0],
            "s2": [2.0, 2.0, 4.0],
            "s3": [np.nan, 3.0, 4.0],
            "s4": [4.0, 4.0, 4.0],
        },
        index=["positive", "negative", "constant"],
    )
    omics2 = pd.DataFrame(
        {
            "s1": [2.0, 4.0, 5.0],
            "s2": [4.0, 3.0, 5.0],
            "s3": [6.0, 2.0, 5.0],
            "s4": [8.0, 1.0, 5.0],
        },
        index=omics1.index,
    )

    result = feature_correlations(omics1, omics2, min_valid_pairs=3)

    assert list(result.index) == ["positive", "negative"]
    assert result.loc["positive", "Gene Correlation"] == pytest.approx(1.0)
    assert result.loc["positive", "N"] == 3
    assert result.loc["negative", "Gene Correlation"] == pytest.approx(-1.0)
    assert result["BH adjusted P"].between(0, 1).all()


def test_sample_and_gene_sample_correlation_axes_are_explicit():
    columns = [f"s{i}" for i in range(6)]
    sample_rho = pd.DataFrame({"Corr": [-0.9, -0.5, -0.1, 0.1, 0.5, 0.9]}, index=columns)
    omics1 = pd.DataFrame(
        [sample_rho["Corr"].to_numpy(), np.arange(6, dtype=float), np.ones(6)],
        index=["driver", "other", "constant_omics1"],
        columns=columns,
    )
    omics2 = pd.DataFrame(
        [-sample_rho["Corr"].to_numpy(), np.arange(6, dtype=float), sample_rho["Corr"].to_numpy()],
        index=omics1.index,
        columns=columns,
    )

    result = gene_sample_correlations(sample_rho, omics1, omics2, min_valid_pairs=4)

    assert result.loc["driver", "Corr_omics1"] == pytest.approx(1.0)
    assert result.loc["driver", "Corr_omics2"] == pytest.approx(-1.0)
    assert result.loc["driver", "N_omics1"] == 6
    assert result.loc["driver", "N_omics2"] == 6
    assert np.isnan(result.loc["constant_omics1", "Corr_omics1"])
    assert result.loc["constant_omics1", "N_omics1"] == 6


def test_sample_correlations_skip_constant_samples():
    omics1 = pd.DataFrame({"positive": [1, 2, 3, 4], "negative": [1, 2, 3, 4], "constant": [1, 1, 1, 1]})
    omics2 = pd.DataFrame({"positive": [2, 4, 6, 8], "negative": [8, 6, 4, 2], "constant": [2, 3, 4, 5]})

    result = sample_correlations(omics1, omics2, min_valid_pairs=4)

    assert list(result.index) == ["positive", "negative"]
    assert result.loc["positive", "Corr"] == pytest.approx(1.0)
    assert result.loc["negative", "Corr"] == pytest.approx(-1.0)


def test_feature_selection_uses_variability_from_both_omics():
    frame = pd.DataFrame(
        {
            "Gene Correlation": [0.9, 0.8, -0.8, -0.9],
            "CV1": [0.1, 0.2, 0.1, 0.2],
            "CV2": [10.0, 0.1, 9.0, 0.1],
        },
        index=["high", "other_high", "low", "other_low"],
    )

    high, low = select_high_low_features(frame, fraction=0.5, max_variable_features=2)

    assert high == ["high"]
    assert low == ["low"]


def test_clustering_is_reproducible_and_ari_is_symmetric(tmp_path: Path):
    matrix = pd.DataFrame(
        {
            "a1": [4.0, 5.0, 1.0],
            "a2": [4.2, 4.8, 1.1],
            "b1": [-4.0, -5.0, -1.0],
            "b2": [-4.1, -4.8, -1.2],
        },
        index=["f1", "f2", "f3"],
    )
    first = cluster_samples(matrix, name="test", output_path=tmp_path / "first.png", n_clusters=2, dpi=72)
    second = cluster_samples(matrix, name="test", output_path=tmp_path / "second.png", n_clusters=2, dpi=72)

    pd.testing.assert_frame_equal(first, second)
    assert first["group"].nunique() == 2
    assert first.loc["a1", "group"] == first.loc["a2", "group"]
    assert first.loc["b1", "group"] == first.loc["b2", "group"]
    ari = adjusted_rand_matrix({"first": first, "second": second})
    assert ari.loc["first", "second"] == pytest.approx(1.0)
    assert ari.equals(ari.T)


def test_config_replay_completes_valid_pair_and_skips_low_overlap(tmp_path: Path):
    rng = np.random.default_rng(42)
    features = [f"f{i}" for i in range(12)]
    samples = [f"s{i}" for i in range(12)]
    omics1 = pd.DataFrame(rng.normal(size=(12, 12)), index=features, columns=samples)
    omics2 = omics1 * 0.7 + pd.DataFrame(rng.normal(scale=0.2, size=(12, 12)), index=features, columns=samples)
    for index in range(12):
        omics1.iloc[index, index] = np.nan
    fusion = pd.DataFrame([rng.normal(size=12)], index=["f0"], columns=samples)

    paths = {}
    for name, matrix in (("rna", omics1), ("protein", omics2), ("fusion", fusion)):
        path = tmp_path / f"{name}.tsv"
        matrix.to_csv(path, sep="\t")
        paths[name] = path

    output_dir = tmp_path / "out"
    config = configparser.ConfigParser()
    config["input"] = {f"{name}_path": str(path) for name, path in paths.items()}
    config["output"] = {"output_dir": str(output_dir)}
    config["settings"] = {
        "min_overlap_features": "10",
        "min_overlap_samples": "10",
        "min_feature_pairs": "10",
        "min_sample_pairs": "10",
        "generate_clustering": "false",
        "save_aligned_matrices": "false",
        "dpi": "72",
    }
    config_path = tmp_path / "config.ini"
    with config_path.open("w", encoding="utf-8") as handle:
        config.write(handle)

    result = run_omicsx_config(config_path, print_json=False)

    assert result["status"] == "completed"
    assert result["completed"] == 1
    assert result["skipped"] == 2
    summary = pd.read_csv(output_dir / "omicsX_batch_summary.tsv", sep="\t")
    assert summary.set_index("pair").loc["rna__protein", "status"] == "completed"
    feature_table = pd.read_csv(
        output_dir / "rna__protein" / "gene_correlation" / "gene_wise_corr.tsv",
        sep="\t",
        index_col=0,
    )
    assert len(feature_table) == 12
    assert (output_dir / "rna__protein" / "visualizations" / "omicsX_pair_overview.png").stat().st_size > 1000
    assert (output_dir / "config.ini").exists()
