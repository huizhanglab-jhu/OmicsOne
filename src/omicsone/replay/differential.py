from __future__ import annotations

import configparser
import json
import shutil
from pathlib import Path
from typing import Any

from omicsone.services.volcano_enrichment import (
    VolcanoEnrichmentResult,
    generate_volcano_enrichment,
    load_volcano_config,
)


DEFAULT_OUTPUT_CONFIG_NAME = "config.ini"
DEFAULT_TEMPLATE_NAME = "default_differential_config.ini"

_INPUT_STRING_KEYS = {
    "normal_path",
    "tumor_path",
    "fasta_path",
    "out_dir",
    "output_dir",
    "cohort",
    "omics",
    "job_name",
    "output_prefix",
    "api_url",
    "title",
    "enrichment_title",
}


def default_differential_config_path() -> Path:
    """Return the packaged default differential-analysis config template."""
    return Path(__file__).with_name(DEFAULT_TEMPLATE_NAME)


def _read_input_ini(input_ini: str | Path) -> dict[str, Any]:
    input_path = Path(input_ini)
    if not input_path.exists():
        raise FileNotFoundError(f"Input INI does not exist: {input_path}")

    parser = configparser.ConfigParser()
    parser.read(input_path, encoding="utf-8-sig")

    payload: dict[str, Any] = {}
    for section in ("input", "files", "analysis", "output", "plot"):
        if not parser.has_section(section):
            continue
        for key in _INPUT_STRING_KEYS:
            if parser.has_option(section, key):
                value = parser.get(section, key).strip()
                if value:
                    payload[key] = value

    # Let the standard OmicsOne parser handle optional typed settings if users put
    # them in input.ini, while keeping the path aliases above available for local replay.
    payload.update(load_volcano_config(input_path))
    return payload


def ensure_output_differential_config(
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

    source_path = Path(template_path) if template_path else default_differential_config_path()
    if not source_path.exists():
        raise FileNotFoundError(f"Default differential config template does not exist: {source_path}")
    shutil.copyfile(source_path, config_path)
    return config_path


def _resolve_diff_output_dir(input_payload: dict[str, Any], explicit_output_dir: str | Path | None) -> Path:
    if explicit_output_dir is not None:
        return Path(explicit_output_dir)

    base_output_dir = input_payload.get("out_dir") or input_payload.get("output_dir")
    if not base_output_dir:
        raise ValueError("input.ini must define out_dir or output_dir under [input], [files], or [output].")

    output_path = Path(str(base_output_dir))
    if output_path.name.lower() == "diff":
        return output_path
    return output_path / "diff"


def _result_to_dict(result: VolcanoEnrichmentResult, config_path: Path, input_path: Path) -> dict[str, Any]:
    return {
        "cohort": result.cohort,
        "omics": result.omics,
        "output_dir": str(result.output_dir),
        "input_ini": str(input_path),
        "config_ini": str(config_path),
        "diff_tsv": str(result.diff_tsv),
        "combined_matrix_tsv": str(result.combined_matrix_tsv),
        "up_genes_tsv": str(result.up_genes_tsv),
        "down_genes_tsv": str(result.down_genes_tsv),
        "total_genes_tsv": str(result.total_genes_tsv),
        "up_enrichment_tsv": str(result.up_enrichment_tsv),
        "down_enrichment_tsv": str(result.down_enrichment_tsv),
        "enrichment_plot_tsv": str(result.enrichment_plot_tsv),
        "volcano_png": str(result.volcano_png),
        "volcano_pdf": str(result.volcano_pdf),
        "volcano_tiff": str(result.volcano_tiff),
        "enrichment_png": str(result.enrichment_png),
        "enrichment_pdf": str(result.enrichment_pdf),
        "enrichment_tiff": str(result.enrichment_tiff),
        "report_html": str(result.report_html),
        "result_log": str(result.result_log),
        "n8n_js": str(result.n8n_js),
        "feature_count": result.feature_count,
        "diff_feature_count": result.diff_feature_count,
        "up_count": result.up_count,
        "down_count": result.down_count,
        "pure_up_gene_count": result.pure_up_gene_count,
        "pure_down_gene_count": result.pure_down_gene_count,
        "total_gene_count": result.total_gene_count,
        "up_enrichment_count": result.up_enrichment_count,
        "down_enrichment_count": result.down_enrichment_count,
        "method": result.method,
        "fdr_cutoff": result.fdr_cutoff,
        "log2fc_cutoff": result.log2fc_cutoff,
        "gene_sets": result.gene_sets,
    }


def run_differential_analysis(
    input_ini: str | Path,
    *,
    output_dir: str | Path | None = None,
    config_name: str = DEFAULT_OUTPUT_CONFIG_NAME,
    template_path: str | Path | None = None,
    print_json: bool = True,
    **overrides: Any,
) -> dict[str, Any]:
    """
    Run differential volcano/enrichment analysis from replay-friendly INI files.

    input.ini supplies data paths and the output directory. The output directory
    owns an editable config.ini; if it is missing, the packaged
    default_differential_config.ini is copied there before execution.
    """
    input_path = Path(input_ini)
    input_payload = _read_input_ini(input_path)
    diff_output_dir = _resolve_diff_output_dir(input_payload, output_dir)
    input_payload["output_dir"] = str(diff_output_dir)
    input_payload.pop("out_dir", None)

    config_path = ensure_output_differential_config(
        input_payload["output_dir"],
        config_name=config_name,
        template_path=template_path,
    )
    config_payload = load_volcano_config(config_path)
    payload = {**config_payload, **input_payload, **overrides}

    missing = [key for key in ("normal_path", "tumor_path", "output_dir") if not payload.get(key)]
    if missing:
        raise ValueError(f"Missing required differential input setting(s): {', '.join(missing)}")

    result = generate_volcano_enrichment(**payload)
    result_dict = _result_to_dict(result, config_path, input_path)
    if print_json:
        print(json.dumps(result_dict, indent=2), flush=True)
    return result_dict


def run_phospho_differential_analysis(
    input_ini: str | Path,
    *,
    output_dir: str | Path | None = None,
    config_name: str = DEFAULT_OUTPUT_CONFIG_NAME,
    template_path: str | Path | None = None,
    print_json: bool = True,
    **overrides: Any,
) -> dict[str, Any]:
    """
    Run phosphosite-level differential volcano/enrichment analysis.

    This is the replay-safe adapter for phosphoproteomics matrices. It keeps
    full phosphosite feature IDs such as ENSG|ENSP|S123 intact unless the caller
    explicitly overrides strip_feature_version.
    """
    phospho_defaults = {
        "omics": "Phospho",
        "strip_feature_version": False,
    }
    phospho_defaults.update(overrides)
    return run_differential_analysis(
        input_ini,
        output_dir=output_dir,
        config_name=config_name,
        template_path=template_path,
        print_json=print_json,
        **phospho_defaults,
    )

