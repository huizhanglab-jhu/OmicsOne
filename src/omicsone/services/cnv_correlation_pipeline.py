from __future__ import annotations

import html
import configparser
import json
import os
import tempfile
import time
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Literal

from omicsone.services.spearman_omics import OmicsInput, compute_paired_omics_spearman
from omicsone.plots.cnv_correlation import (
    CnvCorrelationFigureResult,
    generate_cnv_correlation_figures,
)
from omicsone.replay.cnv_correlation_clean_figures import (
    run_cnv_correlation_clean_figures,
)


PipelineStatus = Literal["queued", "running", "completed", "failed"]


def create_pipeline_job(params: dict[str, Any]) -> dict[str, Any]:
    normalized = _normalize_params(params)
    _validate_pipeline_inputs(normalized)

    output_dir = Path(normalized["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    job_id = _new_job_id(normalized.get("cohort"))
    now = _now()
    record = {
        "job_id": job_id,
        "status": "queued",
        "created_at": now,
        "started_at": None,
        "ended_at": None,
        "current_step": "queued",
        "output_dir": str(output_dir),
        "html_report": None,
        "json_report": None,
        "error": None,
        "request": normalized,
        "jobs": [],
    }
    _write_job_record(job_id, record)
    return record


def run_pipeline_job(job_id: str) -> None:
    record = read_pipeline_job(job_id)
    params = record["request"]
    output_dir = Path(params["output_dir"])
    jobs = _planned_jobs(params)
    results: list[dict[str, Any]] = []
    start_all = time.perf_counter()

    try:
        _update_job(
            job_id,
            status="running",
            started_at=_now(),
            current_step="validating input files",
        )
        _validate_pipeline_inputs(params)
        universe_inputs = _universe_inputs(params)

        for job in jobs:
            job_start = time.perf_counter()
            pair_dir = output_dir / job["name"]
            _update_job(
                job_id,
                current_step=f"computing {job['name']} Spearman correlation",
            )

            pair_result = compute_paired_omics_spearman(
                inputs=[
                    OmicsInput(path=Path(params["cnv_path"]), data_type="cnv"),
                    OmicsInput(
                        path=Path(job["target_path"]),
                        data_type=job["target_type"],
                    ),
                ],
                output_dir=pair_dir,
                min_valid_pairs=params["min_valid_pairs"],
                output_prefix=params["cohort"],
                universe_inputs=universe_inputs,
            )
            pair = pair_result.pairs[0]

            _update_job(
                job_id,
                current_step=f"generating {job['name']} figures",
            )
            figure_result = generate_cnv_correlation_figures(
                correlation_file=pair.correlation_file,
                cnv_for_corr_file=pair_result.matrix_files["cnv"],
                target_for_corr_file=pair_result.matrix_files[job["target_type"]],
                fasta_file=Path(params["fasta_file"]),
                chromosomes_file=Path(params["chromosomes_file"]),
                cytoband_file=Path(params["cytoband_file"]),
                gistic_file=(
                    Path(params["gistic_path"])
                    if params.get("gistic_path") is not None
                    else None
                ),
                output_dir=pair_dir / "figures",
                target_type=job["target_type"],
                correlation_threshold=params["correlation_threshold"],
                output_prefix=f"{params['cohort']}_{job['name']}",
                chunksize=params["chunksize"],
                dpi=params["dpi"],
            )
            clean_figure_result = None
            if params.get("generate_clean_figures") and figure_result.gistic_counts_file:
                _update_job(
                    job_id,
                    current_step=f"generating {job['name']} clean figures",
                )
                clean_config = _write_clean_figure_config(
                    params=params,
                    job=job,
                    figure_result=figure_result,
                )
                clean_figure_result = run_cnv_correlation_clean_figures(
                    clean_config,
                    print_json=False,
                )

            results.append(
                _build_pair_report(
                    job=job,
                    elapsed_seconds=time.perf_counter() - job_start,
                    pair_result=pair_result,
                    correlation_file=pair.correlation_file,
                    result_rows=pair.result_rows,
                    figure_result=figure_result,
                    clean_figure_result=clean_figure_result,
                )
            )
            _update_job(job_id, jobs=results)

        total_elapsed = time.perf_counter() - start_all
        _update_job(job_id, current_step="writing reports")
        report = _build_pipeline_report(
            job_id=job_id,
            params=params,
            jobs=results,
            total_elapsed_seconds=total_elapsed,
        )
        json_report = output_dir / "pipeline_report.json"
        html_report = output_dir / "pipeline_report.html"
        json_report.write_text(json.dumps(report, indent=2), encoding="utf-8")
        html_report.write_text(_render_html_report(report), encoding="utf-8")

        _update_job(
            job_id,
            status="completed",
            ended_at=_now(),
            current_step="completed",
            html_report=str(html_report),
            json_report=str(json_report),
            jobs=results,
        )
    except Exception as exc:  # pragma: no cover - exercised by API jobs
        _update_job(
            job_id,
            status="failed",
            ended_at=_now(),
            current_step="failed",
            error=f"{type(exc).__name__}: {exc}",
        )


def read_pipeline_job(job_id: str) -> dict[str, Any]:
    path = _job_record_path(job_id)
    if not path.is_file():
        raise FileNotFoundError(f"Pipeline job does not exist: {job_id}")
    return json.loads(path.read_text(encoding="utf-8"))


def _normalize_params(params: dict[str, Any]) -> dict[str, Any]:
    normalized = dict(params)
    cohort = normalized.get("cohort") or "omicsone"
    normalized["cohort"] = str(cohort)
    normalized["use_three_way_common_genes"] = bool(
        normalized.get("use_three_way_common_genes", False)
    )
    normalized["generate_clean_figures"] = bool(
        normalized.get("generate_clean_figures", True)
    )
    normalized["output_dir"] = str(Path(normalized["output_dir"]).expanduser())
    for key in (
        "cnv_path",
        "rna_path",
        "protein_path",
        "gistic_path",
        "fasta_file",
        "chromosomes_file",
        "cytoband_file",
    ):
        if normalized.get(key) is not None:
            normalized[key] = str(Path(normalized[key]).expanduser())
    return normalized


def _validate_pipeline_inputs(params: dict[str, Any]) -> None:
    if params.get("rna_path") is None and params.get("protein_path") is None:
        raise ValueError("At least one of rna_path or protein_path must be provided")

    required = [
        "cnv_path",
        "fasta_file",
        "chromosomes_file",
        "cytoband_file",
    ]
    if params.get("rna_path") is not None:
        required.append("rna_path")
    if params.get("protein_path") is not None:
        required.append("protein_path")
    if params.get("gistic_path") is not None:
        required.append("gistic_path")

    missing = [
        f"{key}: {params.get(key)}"
        for key in required
        if not Path(params[key]).is_file()
    ]
    if missing:
        raise FileNotFoundError("Missing input files: " + "; ".join(missing))
    if params.get("use_three_way_common_genes") and (
        params.get("rna_path") is None or params.get("protein_path") is None
    ):
        raise ValueError(
            "use_three_way_common_genes requires cnv_path, rna_path, and protein_path"
        )


def _universe_inputs(params: dict[str, Any]) -> list[OmicsInput] | None:
    if not params.get("use_three_way_common_genes"):
        return None
    return [
        OmicsInput(path=Path(params["cnv_path"]), data_type="cnv"),
        OmicsInput(path=Path(params["protein_path"]), data_type="protein"),
        OmicsInput(path=Path(params["rna_path"]), data_type="rna"),
    ]


def _planned_jobs(params: dict[str, Any]) -> list[dict[str, str]]:
    jobs = []
    if params.get("protein_path") is not None:
        jobs.append(
            {
                "name": "cnv_vs_protein",
                "target_type": "protein",
                "target_path": params["protein_path"],
            }
        )
    if params.get("rna_path") is not None:
        jobs.append(
            {
                "name": "cnv_vs_rna",
                "target_type": "rna",
                "target_path": params["rna_path"],
            }
        )
    return jobs


def _write_clean_figure_config(
    *,
    params: dict[str, Any],
    job: dict[str, str],
    figure_result: CnvCorrelationFigureResult,
) -> Path:
    figures_dir = figure_result.output_dir
    config = configparser.ConfigParser()
    config["inputs"] = {
        "figures_dir": str(figures_dir),
        "annotated_correlations": str(figure_result.annotated_correlations_file),
        "cnv_distribution_counts": str(figure_result.cnv_distribution_counts_file),
        "gistic_counts": str(figure_result.gistic_counts_file or ""),
        "chromosomes_file": params["chromosomes_file"],
    }
    config["output"] = {
        "clean_dir": str(figures_dir),
        "prefix": f"{params['cohort']}_{job['name']}",
        "target_type": job["target_type"],
        "use_simple_names": "true",
    }
    config["plot"] = {
        "font_family": "Liberation Sans",
        "figure_dpi": str(params["dpi"]),
        "tiff_dpi": "600",
        "heatmap_vmin": "-0.5",
        "heatmap_vmax": "0.5",
        "self_cis_mode": "positive_only",
        "local_definition": "same chromosome arm, excluding self-cis",
        "distal_definition": "different chromosome arm",
    }
    config_path = figures_dir / "replay_clean_config.ini"
    with config_path.open("w", encoding="utf-8") as handle:
        config.write(handle)
    return config_path


def _build_pair_report(
    *,
    job: dict[str, str],
    elapsed_seconds: float,
    pair_result: Any,
    correlation_file: Path,
    result_rows: int,
    figure_result: CnvCorrelationFigureResult,
    clean_figure_result: dict[str, Any] | None = None,
) -> dict[str, Any]:
    png_files = {
        "corr_heatmap": figure_result.corr_heatmap_file,
        "cnv_distribution": figure_result.cnv_distribution_file,
        "chromosome": figure_result.chromosome_file,
        "gistic": figure_result.gistic_file,
        "combined": figure_result.combined_file,
    }
    return {
        "name": job["name"],
        "target_type": job["target_type"],
        "elapsed_seconds": elapsed_seconds,
        "backend": pair_result.backend,
        "common_gene_count": pair_result.common_gene_count,
        "common_sample_count": pair_result.common_sample_count,
        "matrix_shapes": {
            data_type: list(shape)
            for data_type, shape in pair_result.matrix_shapes.items()
        },
        "matrix_files": {
            data_type: str(path)
            for data_type, path in pair_result.matrix_files.items()
        },
        "correlation_file": str(correlation_file),
        "result_rows": result_rows,
        "filtered_correlation_count": figure_result.filtered_correlation_count,
        "annotated_correlation_count": figure_result.annotated_correlation_count,
        "annotated_correlations_file": str(
            figure_result.annotated_correlations_file
        ),
        "cnv_distribution_counts_file": str(
            figure_result.cnv_distribution_counts_file
        ),
        "gistic_counts_file": (
            str(figure_result.gistic_counts_file)
            if figure_result.gistic_counts_file is not None
            else None
        ),
        "png_files": _stringify_paths(png_files),
        "pdf_files": _stringify_paths(
            {
                name: path.with_suffix(".pdf") if path is not None else None
                for name, path in png_files.items()
            }
        ),
        "tiff_files": _stringify_paths(
            {
                name: path.with_suffix(".tiff") if path is not None else None
                for name, path in png_files.items()
            }
        ),
        "clean_figures": clean_figure_result,
    }


def _build_pipeline_report(
    *,
    job_id: str,
    params: dict[str, Any],
    jobs: list[dict[str, Any]],
    total_elapsed_seconds: float,
) -> dict[str, Any]:
    used_files = {
        "cnv": params["cnv_path"],
        "rna": params.get("rna_path"),
        "protein": params.get("protein_path"),
        "gistic": params.get("gistic_path"),
        "fasta": params["fasta_file"],
        "chromosomes": params["chromosomes_file"],
        "cytoband": params["cytoband_file"],
    }
    return {
        "job_id": job_id,
        "cohort": params["cohort"],
        "output_dir": params["output_dir"],
        "min_valid_pairs": params["min_valid_pairs"],
        "correlation_threshold": params["correlation_threshold"],
        "dpi": params["dpi"],
        "gene_universe_mode": (
            "cnv_rna_protein_common"
            if params.get("use_three_way_common_genes")
            else "pairwise_common"
        ),
        "total_elapsed_seconds": total_elapsed_seconds,
        "used_files": used_files,
        "steps": [
            "Validate input files.",
            (
                "Use the shared CNV/RNA/protein gene and sample intersection."
                if params.get("use_three_way_common_genes")
                else "Use the pairwise CNV/target gene and sample intersection."
            ),
            "Preprocess CNV and target omics matrices.",
            "Save header/index *_for_corr.txt files.",
            "Convert matrices to numeric-only Rust Spearman inputs.",
            "Compute Spearman correlations with the Rust backend.",
            "Generate CNV correlation figures as PNG, editable-font PDF, and 600-DPI TIFF.",
            (
                "Generate standardized clean CNV correlation figures with simple filenames."
                if params.get("generate_clean_figures")
                else "Skip standardized clean CNV correlation figures."
            ),
            "Write JSON and HTML reports.",
        ],
        "jobs": jobs,
    }


def _render_html_report(report: dict[str, Any]) -> str:
    output_dir = Path(report["output_dir"])

    def rel(path: str | Path | None) -> str:
        if path is None:
            return ""
        p = Path(path)
        try:
            return p.relative_to(output_dir).as_posix()
        except ValueError:
            return p.as_posix()

    def link(path: str | Path | None) -> str:
        if path is None:
            return ""
        value = html.escape(rel(path))
        return f'<a href="{value}">{value}</a>'

    parts = [
        "<!doctype html>",
        "<html><head><meta charset='utf-8'>",
        f"<title>{html.escape(report['cohort'])} CNV Correlation Pipeline</title>",
        "<style>body{font-family:'Liberation Sans',Arial,sans-serif;max-width:1200px;margin:32px auto;line-height:1.45}table{border-collapse:collapse;width:100%;margin:12px 0}th,td{border:1px solid #ddd;padding:6px;vertical-align:top}th{background:#f3f3f3}img{max-width:100%;border:1px solid #ddd;margin:8px 0 18px}.meta{color:#555}</style>",
        "</head><body>",
        f"<h1>{html.escape(report['cohort'])} CNV Correlation Pipeline</h1>",
        f"<p class='meta'>Job ID: {html.escape(report['job_id'])}</p>",
        f"<p class='meta'>Output folder: {html.escape(report['output_dir'])}</p>",
        f"<p class='meta'>Total elapsed: {report['total_elapsed_seconds']:.1f} seconds</p>",
        "<h2>Used Files</h2><table><tr><th>Role</th><th>Path</th><th>Exists</th><th>Size bytes</th></tr>",
    ]
    for name, path_value in report["used_files"].items():
        if path_value is None:
            exists = False
            size = ""
            path_text = ""
        else:
            path = Path(path_value)
            exists = path.exists()
            size = path.stat().st_size if exists and path.is_file() else ""
            path_text = str(path)
        parts.append(
            f"<tr><td>{html.escape(name)}</td><td>{html.escape(path_text)}</td><td>{exists}</td><td>{size}</td></tr>"
        )
    parts.append("</table><h2>Process Steps</h2><ol>")
    for step in report["steps"]:
        parts.append(f"<li>{html.escape(step)}</li>")
    parts.append("</ol>")

    for job in report["jobs"]:
        parts.append(f"<h2>{html.escape(job['name'])}</h2><table>")
        rows = [
            ("Backend", job["backend"]),
            ("Elapsed seconds", f"{job['elapsed_seconds']:.1f}"),
            ("Common genes", job["common_gene_count"]),
            ("Common samples", job["common_sample_count"]),
            ("Correlation rows", job["result_rows"]),
            ("Filtered correlations", job["filtered_correlation_count"]),
            ("Annotated correlations", job["annotated_correlation_count"]),
            ("Correlation file", link(job["correlation_file"])),
            ("Annotated correlations", link(job["annotated_correlations_file"])),
            ("CNV distribution counts", link(job["cnv_distribution_counts_file"])),
            ("GISTIC counts", link(job["gistic_counts_file"])),
        ]
        for data_type, path_value in job["matrix_files"].items():
            rows.append((f"Matrix file: {data_type}", link(path_value)))
        for label, value in rows:
            parts.append(f"<tr><th>{html.escape(str(label))}</th><td>{value}</td></tr>")
        parts.append("</table><h3>Figures</h3>")
        for name, png_path in job["png_files"].items():
            if png_path is None:
                continue
            pdf_path = job["pdf_files"][name]
            tiff_path = job["tiff_files"][name]
            parts.append(
                f"<figure><figcaption>{html.escape(name)} | PDF: {link(pdf_path)} | TIFF: {link(tiff_path)}</figcaption>"
                f"<img src='{html.escape(rel(png_path))}' alt='{html.escape(name)}'></figure>"
            )
    parts.append("</body></html>")
    return "\n".join(parts)


def _update_job(job_id: str, **updates: Any) -> None:
    record = read_pipeline_job(job_id)
    record.update(updates)
    _write_job_record(job_id, record)


def _write_job_record(job_id: str, record: dict[str, Any]) -> None:
    path = _job_record_path(job_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(".tmp")
    tmp_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    tmp_path.replace(path)


def _job_record_path(job_id: str) -> Path:
    return _job_registry_dir() / f"{job_id}.json"


def _job_registry_dir() -> Path:
    configured = os.environ.get("OMICSONE_PIPELINE_JOB_DIR")
    if configured:
        return Path(configured)
    return Path(tempfile.gettempdir()) / "omicsone_pipeline_jobs"


def _new_job_id(cohort: str | None) -> str:
    safe_cohort = "".join(
        char if char.isalnum() or char in ("-", "_") else "_"
        for char in (cohort or "omicsone")
    )
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{timestamp}_{safe_cohort}_{uuid.uuid4().hex[:8]}"


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _stringify_paths(paths: dict[str, Path | None]) -> dict[str, str | None]:
    return {
        name: str(path) if path is not None else None
        for name, path in paths.items()
    }
