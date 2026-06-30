from __future__ import annotations

import configparser
import json
import shutil
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

from omicsone.services.mutation_figures import (
    MutationFigureResult,
    generate_hnsc_mutation_figures,
)


DEFAULT_OUTPUT_CONFIG_NAME = "config.ini"
DEFAULT_TEMPLATE_NAME = "default_mutation_config.ini"
DEFAULT_API_URL = "http://127.0.0.1:8001/api/v1/mutations/heatmap/figures"

_INPUT_STRING_KEYS = {
    "mutation_excel_path",
    "maf_path",
    "out_dir",
    "output_dir",
    "omics",
    "cohort",
    "heatmap_filename",
    "mutation_type_filename",
    "output_prefix",
    "api_url",
}


def default_mutation_config_path() -> Path:
    """Return the packaged default mutation-figure config template."""
    return Path(__file__).with_name(DEFAULT_TEMPLATE_NAME)


def load_mutation_config(config_path: str | Path) -> dict[str, Any]:
    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Config file does not exist: {config_file}")

    parser = configparser.ConfigParser()
    parser.read(config_file, encoding="utf-8-sig")

    payload: dict[str, Any] = {}
    string_keys = {
        "input": ["mutation_excel_path", "maf_path", "omics", "cohort"],
        "files": ["mutation_excel_path", "maf_path"],
        "output": ["output_dir", "out_dir", "heatmap_filename", "mutation_type_filename", "output_prefix"],
        "api": ["api_url"],
    }
    float_keys = {
        "settings": ["mutation_threshold"],
        "analysis": ["mutation_threshold"],
    }

    for section, keys in string_keys.items():
        for key in keys:
            if parser.has_option(section, key):
                value = parser.get(section, key).strip()
                if value:
                    payload[key] = value

    for section, keys in float_keys.items():
        for key in keys:
            if parser.has_option(section, key):
                value = parser.get(section, key).strip()
                if value:
                    payload[key] = float(value)

    return payload


def _read_input_ini(input_ini: str | Path) -> dict[str, Any]:
    input_path = Path(input_ini)
    if not input_path.exists():
        raise FileNotFoundError(f"Input INI does not exist: {input_path}")

    parser = configparser.ConfigParser()
    parser.read(input_path, encoding="utf-8-sig")

    payload: dict[str, Any] = {}
    for section in ("input", "files", "settings", "analysis", "output", "api"):
        if not parser.has_section(section):
            continue
        for key in _INPUT_STRING_KEYS:
            if parser.has_option(section, key):
                value = parser.get(section, key).strip()
                if value:
                    payload[key] = value

    payload.update(load_mutation_config(input_path))
    return payload


def ensure_output_mutation_config(
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

    source_path = Path(template_path) if template_path else default_mutation_config_path()
    if not source_path.exists():
        raise FileNotFoundError(f"Default mutation config template does not exist: {source_path}")
    shutil.copyfile(source_path, config_path)
    return config_path


def _resolve_mutation_output_dir(input_payload: dict[str, Any], explicit_output_dir: str | Path | None) -> Path:
    if explicit_output_dir is not None:
        return Path(explicit_output_dir)

    base_output_dir = input_payload.get("out_dir") or input_payload.get("output_dir")
    if not base_output_dir:
        raise ValueError("input.ini must define output_dir under [input], [files], or [output].")

    output_path = Path(str(base_output_dir))
    if "mutation" in output_path.name.lower():
        return output_path
    return output_path / "mutations"


def _request_payload(payload: dict[str, Any]) -> dict[str, Any]:
    keys = {
        "mutation_excel_path",
        "maf_path",
        "output_dir",
        "mutation_threshold",
        "cohort",
        "heatmap_filename",
        "mutation_type_filename",
        "output_prefix",
    }
    return {key: payload[key] for key in keys if key in payload}


def _result_to_dict(result: MutationFigureResult, config_path: Path, input_path: Path) -> dict[str, Any]:
    return {
        "cohort": result.cohort,
        "omics": "SomaticMutation",
        "output_dir": str(result.heatmap_pdf.parent),
        "input_ini": str(input_path),
        "config_ini": str(config_path),
        "heatmap_pdf": str(result.heatmap_pdf),
        "mutation_type_pdf": str(result.mutation_type_pdf),
        "result_log": str(result.result_log),
        "gene_summary_tsv": str(result.gene_summary_tsv),
        "binary_matrix_tsv": str(result.binary_matrix_tsv),
        "mutation_type_matrix_tsv": str(result.mutation_type_matrix_tsv),
        "encoded_matrix_tsv": str(result.encoded_matrix_tsv),
        "sample_annotations_tsv": str(result.sample_annotations_tsv),
        "sample_annotation_colors_tsv": str(result.sample_annotation_colors_tsv),
        "mutation_color_table_tsv": str(result.mutation_color_table_tsv),
        "filtered_maf_tsv": str(result.filtered_maf_tsv),
        "gene_count": result.gene_count,
        "sample_count": result.sample_count,
        "total_maf_rows": result.total_maf_rows,
        "filtered_maf_rows": result.filtered_maf_rows,
        "found_mutations": json.dumps(result.found_mutations),
        "heatmap_width_inch": result.heatmap_width_inch,
        "heatmap_height_inch": result.heatmap_height_inch,
        "heatmap_aspect": result.heatmap_aspect,
        "mutation_type_width_inch": result.mutation_type_width_inch,
        "mutation_type_height_inch": result.mutation_type_height_inch,
        "mutation_type_aspect": result.mutation_type_aspect,
        "sample_gene_ratio": result.sample_gene_ratio,
    }


def _build_payload(
    input_ini: str | Path,
    *,
    output_dir: str | Path | None = None,
    config_name: str = DEFAULT_OUTPUT_CONFIG_NAME,
    template_path: str | Path | None = None,
    **overrides: Any,
) -> tuple[dict[str, Any], Path, Path]:
    input_path = Path(input_ini)
    input_payload = _read_input_ini(input_path)
    mutation_output_dir = _resolve_mutation_output_dir(input_payload, output_dir)
    input_payload["output_dir"] = str(mutation_output_dir)
    input_payload.pop("out_dir", None)

    config_path = ensure_output_mutation_config(
        input_payload["output_dir"],
        config_name=config_name,
        template_path=template_path,
    )
    config_payload = load_mutation_config(config_path)
    payload = {**config_payload, **input_payload, **overrides}

    missing = [key for key in ("mutation_excel_path", "maf_path", "output_dir") if not payload.get(key)]
    if missing:
        raise ValueError(f"Missing required mutation input setting(s): {', '.join(missing)}")

    return payload, config_path, input_path


def run_mutation_figures(
    input_ini: str | Path,
    *,
    output_dir: str | Path | None = None,
    config_name: str = DEFAULT_OUTPUT_CONFIG_NAME,
    template_path: str | Path | None = None,
    print_json: bool = True,
    **overrides: Any,
) -> dict[str, Any]:
    """
    Run mutation heatmap and mutation-type figures from replay-friendly INI files.

    input.ini supplies data paths and the output directory. The output directory
    owns an editable config.ini; if it is missing, the packaged
    default_mutation_config.ini is copied there before execution.
    """
    payload, config_path, input_path = _build_payload(
        input_ini,
        output_dir=output_dir,
        config_name=config_name,
        template_path=template_path,
        **overrides,
    )
    request_payload = _request_payload(payload)

    result = generate_hnsc_mutation_figures(**request_payload)
    result_dict = _result_to_dict(result, config_path, input_path)
    if print_json:
        print(json.dumps(result_dict, indent=2), flush=True)
    return result_dict


def post_mutation_figures_api(
    input_ini: str | Path,
    *,
    output_dir: str | Path | None = None,
    config_name: str = DEFAULT_OUTPUT_CONFIG_NAME,
    template_path: str | Path | None = None,
    api_url: str | None = None,
    timeout_seconds: float = 600,
    print_json: bool = True,
    **overrides: Any,
) -> dict[str, Any]:
    """POST the mutation-figure payload to a running OmicsOne API."""
    payload, config_path, input_path = _build_payload(
        input_ini,
        output_dir=output_dir,
        config_name=config_name,
        template_path=template_path,
        **overrides,
    )
    target_url = api_url or payload.get("api_url") or DEFAULT_API_URL
    request_payload = _request_payload(payload)
    body = json.dumps(request_payload).encode("utf-8")
    request = urllib.request.Request(
        str(target_url),
        data=body,
        headers={"Content-Type": "application/json"},
        method="POST",
    )

    try:
        with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
            response_payload = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"OmicsOne API returned HTTP {exc.code}: {detail}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Could not reach OmicsOne API at {target_url}: {exc.reason}") from exc

    result_dict = {
        **response_payload,
        "output_dir": request_payload["output_dir"],
        "input_ini": str(input_path),
        "config_ini": str(config_path),
        "api_url": str(target_url),
    }
    if print_json:
        print(json.dumps(result_dict, indent=2), flush=True)
    return result_dict

