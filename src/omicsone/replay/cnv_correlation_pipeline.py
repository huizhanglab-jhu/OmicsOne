from __future__ import annotations

import argparse
import configparser
import json
import os
import shutil
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg", force=True)

from omicsone.services.cnv_correlation_pipeline import (
    create_pipeline_job,
    read_pipeline_job,
    run_pipeline_job,
)


DEFAULT_OUTPUT_CONFIG_NAME = "config.ini"
DEFAULT_TEMPLATE_NAME = "default_cnv_correlation_pipeline_config.ini"


_STRING_KEYS = {
    "cohort",
    "name",
    "cnv_path",
    "rna_path",
    "protein_path",
    "gistic_path",
    "fasta_file",
    "fasta_path",
    "chromosomes_file",
    "chromosomes_path",
    "cytoband_file",
    "cytoband_path",
    "output_dir",
    "out_dir",
}

_INT_KEYS = {"min_valid_pairs", "chunksize", "dpi"}
_FLOAT_KEYS = {"correlation_threshold"}
_BOOL_KEYS = {"generate_clean_figures", "use_three_way_common_genes", "use_simple_names"}


def default_cnv_correlation_pipeline_config_path() -> Path:
    """Return the packaged default CNV correlation pipeline config template."""
    return Path(__file__).with_name(DEFAULT_TEMPLATE_NAME)


def ensure_output_cnv_correlation_pipeline_config(
    output_dir: str | Path,
    *,
    config_name: str = DEFAULT_OUTPUT_CONFIG_NAME,
    template_path: str | Path | None = None,
) -> Path:
    """Create output_dir/config.ini from the packaged template if it is missing."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    config_path = output_path / config_name
    if config_path.exists():
        return config_path

    source_path = Path(template_path) if template_path else default_cnv_correlation_pipeline_config_path()
    if not source_path.exists():
        raise FileNotFoundError(f"Default CNV correlation pipeline config template does not exist: {source_path}")
    shutil.copyfile(source_path, config_path)
    return config_path


def load_cnv_correlation_pipeline_config(config_path: str | Path) -> dict[str, Any]:
    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Config file does not exist: {config_file}")

    parser = configparser.ConfigParser()
    parser.read(config_file, encoding="utf-8-sig")

    payload: dict[str, Any] = {}
    for section in ("task", "input", "files", "paths", "analysis", "settings", "output", "plot"):
        if not parser.has_section(section):
            continue

        for key in _STRING_KEYS:
            if parser.has_option(section, key):
                value = _clean_config_value(parser.get(section, key))
                if value:
                    payload[key] = value

        for key in _INT_KEYS:
            if parser.has_option(section, key):
                value = _clean_config_value(parser.get(section, key))
                if value:
                    payload[key] = int(value)

        for key in _FLOAT_KEYS:
            if parser.has_option(section, key):
                value = _clean_config_value(parser.get(section, key))
                if value:
                    payload[key] = float(value)
        for key in _BOOL_KEYS:
            if parser.has_option(section, key):
                payload[key] = parser.getboolean(section, key)

    return _normalize_payload(payload)


def run_cnv_correlation_pipeline(
    input_ini: str | Path,
    *,
    output_dir: str | Path | None = None,
    config_name: str = DEFAULT_OUTPUT_CONFIG_NAME,
    template_path: str | Path | None = None,
    print_json: bool = True,
    **overrides: Any,
) -> dict[str, Any]:
    """
    Run the CNV correlation pipeline from a replay-friendly INI file.

    The adapter runs the same service used by the FastAPI endpoint, but it
    executes synchronously so a replay workflow can call it as a normal Python
    node or command-line step.
    """
    input_path = Path(input_ini)
    input_payload = load_cnv_correlation_pipeline_config(input_path)

    output_value = output_dir or input_payload.get("output_dir") or input_payload.get("out_dir")
    if not output_value:
        raise ValueError("Config must define output_dir under [paths], [output], or pass --output-dir.")
    resolved_output_dir = Path(output_value)

    ensure_output_cnv_correlation_pipeline_config(
        resolved_output_dir,
        config_name=config_name,
        template_path=template_path,
    )

    payload = {
        "min_valid_pairs": 4,
        "correlation_threshold": 0.6,
        "chunksize": 50000,
        "dpi": 600,
        "use_three_way_common_genes": False,
        "generate_clean_figures": True,
        **input_payload,
        **overrides,
        "output_dir": str(resolved_output_dir),
    }
    payload = _normalize_payload(payload)
    _validate_required_settings(payload)

    previous_job_dir = os.environ.get("OMICSONE_PIPELINE_JOB_DIR")
    os.environ["OMICSONE_PIPELINE_JOB_DIR"] = str(Path(payload["output_dir"]) / "_replay_jobs")
    try:
        job = create_pipeline_job(payload)
        run_pipeline_job(job["job_id"])
        completed = read_pipeline_job(job["job_id"])
    finally:
        if previous_job_dir is None:
            os.environ.pop("OMICSONE_PIPELINE_JOB_DIR", None)
        else:
            os.environ["OMICSONE_PIPELINE_JOB_DIR"] = previous_job_dir

    result = {
        "input_ini": str(input_path),
        "config_ini": str(Path(payload["output_dir"]) / config_name),
        "job_record": str(Path(payload["output_dir"]) / "_replay_jobs" / f"{completed['job_id']}.json"),
        **completed,
    }
    if result.get("status") != "completed":
        raise RuntimeError(f"CNV correlation pipeline failed: {result.get('error')}")
    if print_json:
        print(json.dumps(result, indent=2), flush=True)
    return result


def _clean_config_value(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}:
        return value[1:-1].strip()
    return value


def _normalize_payload(payload: dict[str, Any]) -> dict[str, Any]:
    normalized = dict(payload)
    if "name" in normalized and "cohort" not in normalized:
        normalized["cohort"] = normalized.pop("name")
    if "out_dir" in normalized and "output_dir" not in normalized:
        normalized["output_dir"] = normalized.pop("out_dir")
    if "fasta_path" in normalized and "fasta_file" not in normalized:
        normalized["fasta_file"] = normalized.pop("fasta_path")
    if "chromosomes_path" in normalized and "chromosomes_file" not in normalized:
        normalized["chromosomes_file"] = normalized.pop("chromosomes_path")
    if "cytoband_path" in normalized and "cytoband_file" not in normalized:
        normalized["cytoband_file"] = normalized.pop("cytoband_path")

    for key in ("rna_path", "protein_path", "gistic_path"):
        if normalized.get(key) == "":
            normalized[key] = None
    return normalized


def _validate_required_settings(payload: dict[str, Any]) -> None:
    required = ["cohort", "cnv_path", "output_dir", "fasta_file", "chromosomes_file", "cytoband_file"]
    missing = [key for key in required if not payload.get(key)]
    if not payload.get("rna_path") and not payload.get("protein_path"):
        missing.append("rna_path or protein_path")
    if missing:
        raise ValueError(f"Missing required CNV correlation pipeline setting(s): {', '.join(missing)}")


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the OmicsOne CNV correlation pipeline for local replay.")
    parser.add_argument("--config", required=True, help="Path to the replay input config.ini.")
    parser.add_argument("--output-dir", help="Override output_dir from the config.")
    parser.add_argument("--cohort", help="Override cohort/name from the config.")
    parser.add_argument("--min-valid-pairs", type=int, help="Override minimum valid pairs for Spearman.")
    parser.add_argument("--correlation-threshold", type=float, help="Override figure correlation threshold.")
    parser.add_argument("--chunksize", type=int, help="Override correlation file read chunk size.")
    parser.add_argument("--dpi", type=int, help="Override figure DPI.")
    parser.add_argument(
        "--use-three-way-common-genes",
        action="store_true",
        help="Restrict pairwise correlations to the shared CNV/RNA/protein gene universe.",
    )
    parser.add_argument("--quiet", action="store_true", help="Do not print the JSON result to stdout.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    overrides: dict[str, Any] = {}
    for arg_name, key in (
        ("cohort", "cohort"),
        ("min_valid_pairs", "min_valid_pairs"),
        ("correlation_threshold", "correlation_threshold"),
        ("chunksize", "chunksize"),
        ("dpi", "dpi"),
    ):
        value = getattr(args, arg_name)
        if value is not None:
            overrides[key] = value
    if args.use_three_way_common_genes:
        overrides["use_three_way_common_genes"] = True

    run_cnv_correlation_pipeline(
        args.config,
        output_dir=args.output_dir,
        print_json=not args.quiet,
        **overrides,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

