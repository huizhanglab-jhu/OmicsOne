from __future__ import annotations

import configparser
import json
import shutil
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg", force=True)

from omicsone.services.boxplot_figures import (
    BoxplotFiguresResult,
    generate_boxplot_figures,
)


DEFAULT_OUTPUT_CONFIG_NAME = "config.ini"
DEFAULT_TEMPLATE_NAME = "default_boxplot_config.ini"

_INPUT_STRING_KEYS = {
    "normal_path",
    "tumor_path",
    "fasta_path",
    "out_dir",
    "output_dir",
    "cohort",
    "omics",
    "genes_path",
    "output_prefix",
    "ylabel",
    "font_family",
    "tumor_color",
    "normal_color",
}


def default_boxplot_config_path() -> Path:
    """Return the packaged default boxplot config template."""
    return Path(__file__).with_name(DEFAULT_TEMPLATE_NAME)


def _split_config_list(value: str) -> list[str]:
    return [item.strip() for item in value.replace("\n", ",").split(",") if item.strip()]


def _read_genes_path(genes_path: str | Path) -> list[str]:
    path = Path(genes_path)
    if not path.exists():
        raise FileNotFoundError(f"Gene list file does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"Gene list path is not a file: {path}")

    genes: list[str] = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        genes.extend(_split_config_list(line.replace("\t", ",")))

    if genes and genes[0].strip().lower() in {"gene", "genes", "gene_symbol", "symbol", "feature", "feature_id"}:
        genes = genes[1:]
    if not genes:
        raise ValueError(f"Gene list file does not contain any genes: {path}")
    return genes


def load_boxplot_config(config_path: str | Path) -> dict[str, Any]:
    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Config file does not exist: {config_file}")

    parser = configparser.ConfigParser()
    parser.read(config_file, encoding="utf-8-sig")

    payload: dict[str, Any] = {}
    string_keys = {
        "input": ["normal_path", "tumor_path", "fasta_path", "cohort", "omics", "genes_path"],
        "output": ["output_dir", "output_prefix"],
        "plot": ["ylabel", "font_family", "tumor_color", "normal_color"],
    }
    float_keys = {"plot": ["width", "height", "font_size"]}
    int_keys = {"plot": ["dpi"]}
    bool_keys = {"plot": ["editable_pdf_text"]}
    bool_keys["output"] = ["write_n8n_script"]
    list_keys = {"genes": ["genes"], "analysis": ["genes"], "plot": ["genes"]}

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

    for section, keys in int_keys.items():
        for key in keys:
            if parser.has_option(section, key):
                value = parser.get(section, key).strip()
                if value:
                    payload[key] = int(value)

    for section, keys in bool_keys.items():
        for key in keys:
            if parser.has_option(section, key):
                payload[key] = parser.getboolean(section, key)

    for section, keys in list_keys.items():
        for key in keys:
            if parser.has_option(section, key):
                values = _split_config_list(parser.get(section, key))
                if values:
                    payload[key] = values

    return payload


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

    payload.update(load_boxplot_config(input_path))
    return payload


def ensure_output_boxplot_config(
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

    source_path = Path(template_path) if template_path else default_boxplot_config_path()
    if not source_path.exists():
        raise FileNotFoundError(f"Default boxplot config template does not exist: {source_path}")
    shutil.copyfile(source_path, config_path)
    return config_path


def _resolve_boxplot_output_dir(input_payload: dict[str, Any], explicit_output_dir: str | Path | None) -> Path:
    if explicit_output_dir is not None:
        return Path(explicit_output_dir)

    base_output_dir = input_payload.get("out_dir") or input_payload.get("output_dir")
    if not base_output_dir:
        raise ValueError("input.ini must define out_dir or output_dir under [input], [files], or [output].")

    output_path = Path(str(base_output_dir))
    if "boxplot" in output_path.name.lower():
        return output_path
    return output_path / "boxplots"


def _result_to_dict(result: BoxplotFiguresResult, config_path: Path, input_path: Path) -> dict[str, Any]:
    boxplot_pdfs = [str(path) for path in result.boxplot_pdfs]
    records = [
        {
            "gene": record.gene,
            "gene_id": record.gene_id,
            "pdf": str(record.pdf),
            "pvalue": record.pvalue,
            "significance": record.significance,
            "normal_count": record.normal_count,
            "tumor_count": record.tumor_count,
        }
        for record in result.records
    ]
    return {
        "cohort": result.cohort,
        "omics": result.omics,
        "output_dir": str(result.output_dir),
        "input_ini": str(input_path),
        "config_ini": str(config_path),
        "boxplot_pdfs": json.dumps(boxplot_pdfs),
        "summary_tsv": str(result.summary_tsv),
        "result_log": str(result.result_log),
        "n8n_js": str(result.n8n_js),
        "generated_count": result.generated_count,
        "missing_genes": json.dumps(result.missing_genes),
        "records": json.dumps(records),
    }


def run_boxplot_figures(
    input_ini: str | Path,
    *,
    output_dir: str | Path | None = None,
    config_name: str = DEFAULT_OUTPUT_CONFIG_NAME,
    template_path: str | Path | None = None,
    print_json: bool = True,
    **overrides: Any,
) -> dict[str, Any]:
    """
    Run normal-vs-tumor boxplot figure generation from replay-friendly INI files.

    input.ini supplies data paths and the output directory. The output directory
    owns an editable config.ini; if it is missing, the packaged
    default_boxplot_config.ini is copied there before execution.
    """
    input_path = Path(input_ini)
    input_payload = _read_input_ini(input_path)
    boxplot_output_dir = _resolve_boxplot_output_dir(input_payload, output_dir)
    input_payload["output_dir"] = str(boxplot_output_dir)
    input_payload.pop("out_dir", None)

    config_path = ensure_output_boxplot_config(
        input_payload["output_dir"],
        config_name=config_name,
        template_path=template_path,
    )
    config_payload = load_boxplot_config(config_path)
    payload = {**config_payload, **input_payload, **overrides}

    missing = [key for key in ("normal_path", "tumor_path", "output_dir") if not payload.get(key)]
    if not payload.get("genes") and not payload.get("genes_path"):
        missing.append("genes_path")
    if missing:
        raise ValueError(f"Missing required boxplot input setting(s): {', '.join(missing)}")

    if not payload.get("genes"):
        payload["genes"] = _read_genes_path(payload["genes_path"])
    payload.pop("genes_path", None)

    result = generate_boxplot_figures(**payload)
    result_dict = _result_to_dict(result, config_path, input_path)
    if print_json:
        print(json.dumps(result_dict, indent=2), flush=True)
    return result_dict

