from __future__ import annotations

import configparser
import itertools
import json
import shutil
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from omicsone.omicsx.analysis import (
    align_matrices,
    feature_correlations,
    gene_sample_correlations,
    load_omics_matrix,
    sample_correlations,
    select_high_low_features,
)
from omicsone.omicsx.clustering import adjusted_rand_matrix, cluster_samples
from omicsone.omicsx.plotting import (
    plot_all_pair_medians,
    plot_cluster_ari,
    plot_feature_correlation,
    plot_gene_sample_correlation,
    plot_pair_overview,
    plot_sample_correlation,
    plot_top_feature_correlations,
)


@dataclass(frozen=True)
class OmicsXSettings:
    min_overlap_features: int = 10
    min_overlap_samples: int = 10
    min_feature_pairs: int = 10
    min_sample_pairs: int = 10
    high_low_fraction: float = 0.05
    max_variable_features: int = 1000
    n_clusters: int = 2
    dpi: int = 180
    generate_clustering: bool = True
    save_aligned_matrices: bool = True
    fail_on_pair_error: bool = True

    def validate(self) -> None:
        integer_minimums = {
            "min_overlap_features": (self.min_overlap_features, 1),
            "min_overlap_samples": (self.min_overlap_samples, 1),
            "min_feature_pairs": (self.min_feature_pairs, 2),
            "min_sample_pairs": (self.min_sample_pairs, 2),
            "max_variable_features": (self.max_variable_features, 1),
            "n_clusters": (self.n_clusters, 2),
            "dpi": (self.dpi, 72),
        }
        for name, (value, minimum) in integer_minimums.items():
            if value < minimum:
                raise ValueError(f"{name} must be at least {minimum}")
        if not 0 < self.high_low_fraction <= 0.5:
            raise ValueError("high_low_fraction must be greater than 0 and no more than 0.5")


def run_omicsx_config(
    config_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    print_json: bool = True,
) -> dict[str, Any]:
    """Run all configured OmicsX pairs and return a replay manifest summary."""
    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"OmicsX config does not exist: {config_file}")
    parser = configparser.ConfigParser()
    parser.read(config_file, encoding="utf-8-sig")
    inputs = _read_inputs(parser)
    settings = _read_settings(parser)
    settings.validate()
    resolved_output = Path(output_dir) if output_dir else _configured_output_dir(parser)
    resolved_output.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(config_file, resolved_output / "input_config.ini")
    _write_effective_config(parser, resolved_output, settings)

    matrices: dict[str, pd.DataFrame] = {}
    input_records: list[dict[str, Any]] = []
    for name, path in inputs.items():
        matrix = load_omics_matrix(path)
        matrices[name] = matrix
        input_records.append(
            {
                "name": name,
                "path": str(path),
                "features": int(matrix.shape[0]),
                "samples": int(matrix.shape[1]),
            }
        )

    pair_records: list[dict[str, Any]] = []
    metrics: list[dict[str, Any]] = []
    requested_pairs = _read_pairs(parser, list(inputs))
    for omics1, omics2 in requested_pairs:
        pair_name = f"{omics1}__{omics2}"
        pair_dir = resolved_output / pair_name
        try:
            record, metric = _run_pair(
                matrices[omics1],
                matrices[omics2],
                omics1=omics1,
                omics2=omics2,
                pair_dir=pair_dir,
                settings=settings,
            )
            pair_records.append(record)
            if metric is not None:
                metrics.append(metric)
        except Exception as exc:
            pair_dir.mkdir(parents=True, exist_ok=True)
            error_text = traceback.format_exc()
            (pair_dir / "ERROR.txt").write_text(error_text, encoding="utf-8")
            pair_records.append(
                {
                    "pair": pair_name,
                    "omics1": omics1,
                    "omics2": omics2,
                    "status": "failed",
                    "error": str(exc),
                    "output_dir": str(pair_dir),
                }
            )

    metrics_frame = pd.DataFrame(metrics)
    summary_frame = pd.DataFrame(pair_records)
    summary_frame.to_csv(resolved_output / "omicsX_batch_summary.tsv", sep="\t", index=False)
    metrics_frame.to_csv(resolved_output / "pair_visualization_metrics.tsv", sep="\t", index=False)
    if not metrics_frame.empty:
        plot_all_pair_medians(metrics_frame, output_dir=resolved_output, dpi=settings.dpi)

    manifest = {
        "config": str(config_file),
        "effective_config": str(resolved_output / "config.ini"),
        "output_dir": str(resolved_output),
        "settings": asdict(settings),
        "inputs": input_records,
        "pairs": pair_records,
        "counts": {
            status: sum(record.get("status") == status for record in pair_records)
            for status in ("completed", "skipped", "failed")
        },
    }
    manifest_path = resolved_output / "omicsX_batch_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    result = {
        "status": "failed" if manifest["counts"]["failed"] else "completed",
        "output_dir": str(resolved_output),
        "manifest": str(manifest_path),
        "summary": str(resolved_output / "omicsX_batch_summary.tsv"),
        **manifest["counts"],
    }
    if print_json:
        print(json.dumps(result, indent=2), flush=True)
    if settings.fail_on_pair_error and manifest["counts"]["failed"]:
        failed_names = [record["pair"] for record in pair_records if record.get("status") == "failed"]
        raise RuntimeError(f"OmicsX failed for pair(s): {', '.join(failed_names)}")
    return result


def _run_pair(
    matrix1: pd.DataFrame,
    matrix2: pd.DataFrame,
    *,
    omics1: str,
    omics2: str,
    pair_dir: Path,
    settings: OmicsXSettings,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    aligned1, aligned2 = align_matrices(matrix1, matrix2)
    pair_name = f"{omics1}__{omics2}"
    record: dict[str, Any] = {
        "pair": pair_name,
        "omics1": omics1,
        "omics2": omics2,
        "overlap_features": int(aligned1.shape[0]),
        "overlap_samples": int(aligned1.shape[1]),
        "output_dir": str(pair_dir),
    }
    if aligned1.shape[0] < settings.min_overlap_features or aligned1.shape[1] < settings.min_overlap_samples:
        record.update(
            {
                "status": "skipped",
                "reason": (
                    f"Requires at least {settings.min_overlap_features} overlapping features and "
                    f"{settings.min_overlap_samples} overlapping samples."
                ),
            }
        )
        return record, None

    gene_dir = pair_dir / "gene_correlation"
    sample_dir = pair_dir / "sample_correlation"
    cluster_dir = pair_dir / "sample_clustering"
    visualization_dir = pair_dir / "visualizations"
    for directory in (gene_dir, sample_dir, cluster_dir, visualization_dir):
        directory.mkdir(parents=True, exist_ok=True)

    pair_info = {
        "pair": pair_name,
        "omics1": omics1,
        "omics2": omics2,
        "omics1_source": "left side of pair name and x-axis in gene-sample plots",
        "omics2_source": "right side of pair name and y-axis in gene-sample plots",
    }
    (pair_dir / "pair_info.json").write_text(json.dumps(pair_info, indent=2), encoding="utf-8")
    if settings.save_aligned_matrices:
        aligned1.to_csv(pair_dir / f"{omics1}_aligned.tsv", sep="\t")
        aligned2.to_csv(pair_dir / f"{omics2}_aligned.tsv", sep="\t")

    feature_corr = feature_correlations(aligned1, aligned2, min_valid_pairs=settings.min_feature_pairs)
    if feature_corr.empty:
        record.update({"status": "skipped", "reason": "No features passed pairwise-complete correlation filters."})
        return record, None
    feature_corr.to_csv(gene_dir / "gene_wise_corr.tsv", sep="\t")
    plot_feature_correlation(
        feature_corr,
        omics1=omics1,
        omics2=omics2,
        output_path=gene_dir / "gene_correlation_hist.png",
        dpi=settings.dpi,
    )

    sample_corr = sample_correlations(aligned1, aligned2, min_valid_pairs=settings.min_sample_pairs)
    if sample_corr.empty:
        record.update({"status": "skipped", "reason": "No samples passed pairwise-complete correlation filters."})
        return record, None
    sample_corr.to_csv(sample_dir / "sample_wise_correlation.tsv", sep="\t")
    gene_sample_corr = gene_sample_correlations(
        sample_corr,
        aligned1,
        aligned2,
        min_valid_pairs=settings.min_feature_pairs,
    )
    gene_sample_corr.to_csv(sample_dir / "gene_sample_correlation.tsv", sep="\t")
    plot_sample_correlation(
        sample_corr,
        omics1=omics1,
        omics2=omics2,
        output_path=sample_dir / "sample_correlation_bar.png",
        dpi=settings.dpi,
    )
    plot_gene_sample_correlation(
        gene_sample_corr,
        omics1=omics1,
        omics2=omics2,
        output_path=sample_dir / "gene_sample_correlation_scatter.png",
        dpi=settings.dpi,
    )
    plot_top_feature_correlations(
        feature_corr,
        omics1=omics1,
        omics2=omics2,
        output_path=visualization_dir / "top_feature_correlations.png",
        dpi=settings.dpi,
    )

    high_features, low_features = select_high_low_features(
        feature_corr,
        fraction=settings.high_low_fraction,
        max_variable_features=settings.max_variable_features,
    )
    pd.Series(high_features, name="feature_high").to_csv(cluster_dir / "feature_high.tsv", sep="\t", index=False)
    pd.Series(low_features, name="feature_low").to_csv(cluster_dir / "feature_low.tsv", sep="\t", index=False)

    cluster_tables: dict[str, pd.DataFrame] = {}
    if settings.generate_clustering:
        jobs = (
            ("high_omics1", aligned1, high_features),
            ("high_omics2", aligned2, high_features),
            ("low_omics1", aligned1, low_features),
            ("low_omics2", aligned2, low_features),
        )
        for label, matrix, features in jobs:
            if len(features) < 2:
                continue
            table = cluster_samples(
                matrix.loc[features],
                name=f"{label}: omics1={omics1}; omics2={omics2}",
                output_path=cluster_dir / f"{label}.png",
                n_clusters=settings.n_clusters,
                dpi=settings.dpi,
            )
            table.to_csv(cluster_dir / f"{label}_clusters.tsv", sep="\t")
            cluster_tables[label] = table

    ari: pd.DataFrame | None = None
    if len(cluster_tables) >= 2:
        ari = adjusted_rand_matrix(cluster_tables)
        ari.to_csv(cluster_dir / "cluster_ARI.tsv", sep="\t")
        plot_cluster_ari(
            ari,
            omics1=omics1,
            omics2=omics2,
            output_path=cluster_dir / "ARI.png",
            dpi=settings.dpi,
        )
    plot_pair_overview(
        feature_corr,
        sample_corr,
        gene_sample_corr,
        ari,
        omics1=omics1,
        omics2=omics2,
        output_path=visualization_dir / "omicsX_pair_overview.png",
        dpi=settings.dpi,
    )

    record.update(
        {
            "status": "completed",
            "feature_corr_features": int(len(feature_corr)),
            "sample_corr_samples": int(len(sample_corr)),
            "gene_sample_features": int(len(gene_sample_corr)),
            "high_feature_count": int(len(high_features)),
            "low_feature_count": int(len(low_features)),
            "cluster_view_count": int(len(cluster_tables)),
        }
    )
    metric = {
        "pair": pair_name,
        "omics1": omics1,
        "omics2": omics2,
        "gene_corr_median": float(feature_corr["Gene Correlation"].median()),
        "gene_corr_mean": float(feature_corr["Gene Correlation"].mean()),
        "sample_corr_median": float(sample_corr["Corr"].median()),
        "sample_corr_mean": float(sample_corr["Corr"].mean()),
    }
    return record, metric


def _read_inputs(parser: configparser.ConfigParser) -> dict[str, Path]:
    if not parser.has_section("input"):
        raise ValueError("OmicsX config must contain an [input] section")
    inputs: dict[str, Path] = {}
    for key, raw_value in parser.items("input"):
        if not key.endswith("_path") or key == "metadata_path":
            continue
        value = _clean_value(raw_value)
        if value:
            name = key[: -len("_path")]
            path = Path(value)
            if not path.exists():
                raise FileNotFoundError(f"Input matrix for {name} does not exist: {path}")
            inputs[name] = path
    if len(inputs) < 2:
        raise ValueError("OmicsX config must define at least two matrix *_path entries under [input]")
    return inputs


def _read_settings(parser: configparser.ConfigParser) -> OmicsXSettings:
    defaults = OmicsXSettings()
    if not parser.has_section("settings"):
        return defaults
    section = parser["settings"]
    return OmicsXSettings(
        min_overlap_features=section.getint("min_overlap_features", fallback=defaults.min_overlap_features),
        min_overlap_samples=section.getint("min_overlap_samples", fallback=defaults.min_overlap_samples),
        min_feature_pairs=section.getint("min_feature_pairs", fallback=defaults.min_feature_pairs),
        min_sample_pairs=section.getint("min_sample_pairs", fallback=defaults.min_sample_pairs),
        high_low_fraction=section.getfloat("high_low_fraction", fallback=defaults.high_low_fraction),
        max_variable_features=section.getint("max_variable_features", fallback=defaults.max_variable_features),
        n_clusters=section.getint("n_clusters", fallback=defaults.n_clusters),
        dpi=section.getint("dpi", fallback=defaults.dpi),
        generate_clustering=section.getboolean("generate_clustering", fallback=defaults.generate_clustering),
        save_aligned_matrices=section.getboolean("save_aligned_matrices", fallback=defaults.save_aligned_matrices),
        fail_on_pair_error=section.getboolean("fail_on_pair_error", fallback=defaults.fail_on_pair_error),
    )


def _read_pairs(parser: configparser.ConfigParser, input_names: list[str]) -> list[tuple[str, str]]:
    raw_pairs = parser.get("settings", "pairs", fallback="").strip()
    if not raw_pairs:
        return list(itertools.combinations(input_names, 2))
    pairs: list[tuple[str, str]] = []
    for raw_pair in raw_pairs.split(","):
        separator = ":" if ":" in raw_pair else "__"
        parts = [part.strip() for part in raw_pair.split(separator)]
        if len(parts) != 2 or not all(parts):
            raise ValueError(f"Invalid OmicsX pair '{raw_pair}'. Use omics1:omics2.")
        if any(part not in input_names for part in parts):
            raise ValueError(f"Unknown OmicsX input in pair '{raw_pair}'")
        pairs.append((parts[0], parts[1]))
    return pairs


def _configured_output_dir(parser: configparser.ConfigParser) -> Path:
    if not parser.has_option("output", "output_dir"):
        raise ValueError("OmicsX config must define output_dir under [output] or pass --output-dir")
    return Path(_clean_value(parser.get("output", "output_dir")))


def _write_effective_config(
    source: configparser.ConfigParser,
    output_dir: Path,
    settings: OmicsXSettings,
) -> None:
    effective = configparser.ConfigParser()
    for section in source.sections():
        effective[section] = dict(source.items(section))
    if not effective.has_section("output"):
        effective.add_section("output")
    effective.set("output", "output_dir", str(output_dir))
    if not effective.has_section("settings"):
        effective.add_section("settings")
    for key, value in asdict(settings).items():
        effective.set("settings", key, str(value).lower() if isinstance(value, bool) else str(value))
    with (output_dir / "config.ini").open("w", encoding="utf-8") as handle:
        effective.write(handle)


def _clean_value(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}:
        return value[1:-1].strip()
    return value
