from __future__ import annotations

import configparser
import json
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib

from omicsone.services.pathway_scatter import (
    PathwayScatterResult,
    _matched_highlight_label,
    _missing_highlights,
    _plot_points,
    _read_gene_map,
    _read_gmt_terms,
    _read_highlight_sites,
    _feature_gene_id,
    load_pathway_scatter_config,
    run_pathway_scatter_analysis,
)


DEFAULT_OUTPUT_CONFIG_NAME = "pathway_scatter_config.ini"
DEFAULT_TEMPLATE_NAME = "default_pathway_scatter_config.ini"

_INPUT_STRING_KEYS = {
    "protein_diff_path",
    "phospho_diff_path",
    "phospho_run_dir",
    "pure_up_genes_path",
    "pure_down_genes_path",
    "fasta_path",
    "gmt_path",
    "out_dir",
    "output_dir",
    "cohort",
    "omics",
    "pathway_table_path",
    "highlight_table_path",
    "highlights_table_path",
    "source_mode",
    "protein_log2fc_column",
    "phospho_log2fc_column",
    "gene_selection_mode",
    "title_template",
    "font_family",
}

PATHWAY_SELECTION_MODES = {"filtered", "all_pathway"}


def default_pathway_scatter_config_path() -> Path:
    """Return the packaged default pathway-scatter config template."""
    return Path(__file__).with_name(DEFAULT_TEMPLATE_NAME)


def _split_config_list(value: str) -> list[str]:
    return [item.strip() for item in value.replace("\n", ",").split(",") if item.strip()]


def _read_table(path: str | Path) -> pd.DataFrame:
    table_path = Path(path)
    if not table_path.exists():
        raise FileNotFoundError(f"Pathway/highlight table does not exist: {table_path}")
    sep = "," if table_path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(table_path, sep=sep)


def _pathway_payload_from_table(path: str | Path) -> dict[str, dict[str, str]]:
    table = _read_table(path)
    column_map = {column.strip().lower(): column for column in table.columns}

    folder_column = (
        column_map.get("pathway_folder")
        or column_map.get("folder_name")
        or column_map.get("folder")
        or column_map.get("name")
    )
    term_column = (
        column_map.get("pathway_term")
        or column_map.get("term")
        or column_map.get("pathway")
        or column_map.get("gmt_term")
    )
    highlight_column = (
        column_map.get("highlight_path")
        or column_map.get("highlights_path")
        or column_map.get("highlight_sites_path")
        or column_map.get("highlight_file")
    )

    missing = []
    if not term_column:
        missing.append("pathway_term")
    if not highlight_column:
        missing.append("highlight_path")
    if missing:
        raise ValueError(
            "Pathway/highlight table must contain columns for "
            f"{', '.join(missing)}. Optional folder columns: pathway_folder, folder_name, folder."
        )

    pathways: dict[str, str] = {}
    highlight_paths: dict[str, str] = {}
    for row_index, row in table.iterrows():
        term = str(row[term_column]).strip()
        highlight_path = str(row[highlight_column]).strip()
        folder = str(row[folder_column]).strip() if folder_column else term
        if not term or not highlight_path or term.lower() == "nan" or highlight_path.lower() == "nan":
            continue
        if not folder or folder.lower() == "nan":
            folder = term
        if folder in pathways:
            raise ValueError(f"Duplicate pathway folder '{folder}' in row {row_index + 2} of {path}")
        pathways[folder] = term
        highlight_paths[folder] = highlight_path

    if not pathways:
        raise ValueError(f"Pathway/highlight table does not contain any usable rows: {path}")
    return {"pathways": pathways, "highlight_paths": highlight_paths}


def load_pathway_scatter_replay_config(config_path: str | Path) -> dict[str, Any]:
    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Config file does not exist: {config_file}")

    parser = configparser.ConfigParser()
    parser.optionxform = str
    parser.read(config_file, encoding="utf-8-sig")

    payload: dict[str, Any] = {}
    string_keys = {
        "input": [
            "protein_diff_path",
            "phospho_diff_path",
            "phospho_run_dir",
            "pure_up_genes_path",
            "pure_down_genes_path",
            "fasta_path",
            "gmt_path",
            "pathway_table_path",
            "highlight_table_path",
            "highlights_table_path",
            "source_mode",
            "protein_log2fc_column",
            "phospho_log2fc_column",
        ],
        "output": ["output_dir", "out_dir"],
        "selection": ["gene_selection_mode"],
        "plot": ["title_template", "font_family"],
    }
    float_keys = {"plot": ["width", "height", "font_scale", "point_size", "highlight_point_size"]}
    int_keys = {"plot": ["dpi", "tiff_dpi"]}
    bool_keys = {"plot": ["show_title"]}

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

    if parser.has_section("pathways"):
        pathways = {key: value for key, value in parser.items("pathways") if key.strip() and value.strip()}
        if pathways:
            payload["pathways"] = pathways
    if parser.has_section("highlights"):
        highlight_paths = {key: value for key, value in parser.items("highlights") if key.strip() and value.strip()}
        if highlight_paths:
            payload["highlight_paths"] = highlight_paths
    if parser.has_section("tables"):
        table_paths = {key: value for key, value in parser.items("tables") if key.strip() and value.strip()}
        if table_paths:
            payload["table_paths"] = table_paths
    if parser.has_section("directions"):
        directions = {key: value.strip().lower() for key, value in parser.items("directions") if key.strip() and value.strip()}
        if directions:
            payload["directions"] = directions

    table_path = payload.get("pathway_table_path") or payload.get("highlight_table_path") or payload.get("highlights_table_path")
    if table_path:
        payload.update(_pathway_payload_from_table(table_path))

    return payload


def _read_input_ini(input_ini: str | Path) -> dict[str, Any]:
    input_path = Path(input_ini)
    if not input_path.exists():
        raise FileNotFoundError(f"Input INI does not exist: {input_path}")

    parser = configparser.ConfigParser()
    parser.optionxform = str
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

    payload.update(load_pathway_scatter_replay_config(input_path))
    return payload


def ensure_output_pathway_scatter_config(
    output_dir: str | Path,
    *,
    config_name: str = DEFAULT_OUTPUT_CONFIG_NAME,
    template_path: str | Path | None = None,
) -> Path:
    """Create output_dir/pathway_scatter_config.ini from the packaged template if missing."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    config_path = output_path / config_name
    if config_path.exists():
        return config_path

    source_path = Path(template_path) if template_path else default_pathway_scatter_config_path()
    if not source_path.exists():
        raise FileNotFoundError(f"Default pathway-scatter config template does not exist: {source_path}")
    shutil.copyfile(source_path, config_path)
    return config_path


def _resolve_pathway_output_dir(input_payload: dict[str, Any], explicit_output_dir: str | Path | None) -> Path:
    if explicit_output_dir is not None:
        return Path(explicit_output_dir)

    base_output_dir = input_payload.get("out_dir") or input_payload.get("output_dir")
    if not base_output_dir:
        raise ValueError("input.ini must define out_dir or output_dir under [input], [files], or [output].")

    output_path = Path(str(base_output_dir))
    lower_name = output_path.name.lower()
    if "pathway" in lower_name or "scatter" in lower_name:
        return output_path
    return output_path / "pathway_scatter"


def _merge_payloads(*payloads: dict[str, Any]) -> dict[str, Any]:
    merged: dict[str, Any] = {}
    for payload in payloads:
        for key, value in payload.items():
            if value is None:
                continue
            if isinstance(value, dict) and not value:
                continue
            if isinstance(value, list) and not value:
                continue
            merged[key] = value
    return merged


def _result_to_dict(results: list[PathwayScatterResult], config_path: Path, input_path: Path, output_dir: Path) -> dict[str, Any]:
    records = [
        {
            "pathway": result.pathway,
            "output_dir": str(result.output_dir),
            "points_tsv": str(result.points_tsv),
            "highlight_tsv": str(result.highlight_tsv),
            "missing_highlights_tsv": str(result.missing_highlights_tsv),
            "highlight_txt": str(result.highlight_txt),
            "png": str(result.png),
            "pdf": str(result.pdf),
            "tiff": str(result.tiff),
            "point_count": result.point_count,
            "highlight_count": result.highlight_count,
            "trend_slope": result.trend_slope,
            "trend_intercept": result.trend_intercept,
            "trend_r": result.trend_r,
        }
        for result in results
    ]
    return {
        "output_dir": str(output_dir),
        "input_ini": str(input_path),
        "config_ini": str(config_path),
        "summary_tsv": str(output_dir / "pathway_scatter_summary.tsv"),
        "pathway_count": len(results),
        "records": json.dumps(records),
    }


def run_pathway_scatter_plots(
    input_ini: str | Path,
    *,
    output_dir: str | Path | None = None,
    config_name: str = DEFAULT_OUTPUT_CONFIG_NAME,
    template_path: str | Path | None = None,
    print_json: bool = True,
    **overrides: Any,
) -> dict[str, Any]:
    """
    Run protein-vs-phosphosite pathway scatter plots from replay-friendly INI files.

    The node input can define [pathways] and [highlights] sections directly, or
    provide a pathway/highlight TSV with columns pathway_term and highlight_path.
    Optional folder columns include pathway_folder, folder_name, folder, or name.
    """
    input_path = Path(input_ini)
    input_payload = _read_input_ini(input_path)
    pathway_output_dir = _resolve_pathway_output_dir(input_payload, output_dir)
    input_payload["output_dir"] = str(pathway_output_dir)
    input_payload.pop("out_dir", None)

    config_path = ensure_output_pathway_scatter_config(
        input_payload["output_dir"],
        config_name=config_name,
        template_path=template_path,
    )
    config_payload = load_pathway_scatter_replay_config(config_path)
    payload = _merge_payloads(config_payload, input_payload, overrides)

    missing = [key for key in ("protein_diff_path", "phospho_diff_path", "output_dir") if not payload.get(key)]
    if not payload.get("pathways"):
        missing.append("pathways or pathway_table_path")
    if not payload.get("highlight_paths"):
        missing.append("highlights or highlight_table_path")
    if missing:
        raise ValueError(f"Missing required pathway-scatter input setting(s): {', '.join(missing)}")

    payload.pop("pathway_table_path", None)
    payload.pop("highlight_table_path", None)
    payload.pop("highlights_table_path", None)

    results = run_pathway_scatter_analysis(**payload)
    result_dict = _result_to_dict(results, config_path, input_path, pathway_output_dir)
    if print_json:
        print(json.dumps(result_dict, indent=2), flush=True)
    return result_dict


def run_pathway_scatter_analysis_from_config(config_path: str | Path) -> list[PathwayScatterResult]:
    """Compatibility wrapper for direct service-style config files."""
    return run_pathway_scatter_analysis(**load_pathway_scatter_config(config_path))


def _safe_name(value: str) -> str:
    return value.replace(" ", "_").replace("-", "_").lower()


def _read_gene_list(path: str | Path) -> set[str]:
    gene_path = Path(path)
    if not gene_path.exists():
        raise FileNotFoundError(f"Gene list file does not exist: {gene_path}")
    genes = pd.read_csv(gene_path, sep="\t")
    if "gene" in genes.columns:
        values = genes["gene"]
    else:
        values = genes.iloc[:, 0]
    return set(values.dropna().astype(str).str.strip())


def _diff_with_gene(path: str | Path, gene_map: dict[str, str]) -> pd.DataFrame:
    diff = pd.read_csv(path, sep="\t")
    if "Feature" not in diff.columns:
        diff = diff.rename(columns={diff.columns[0]: "Feature"})
    diff["Gene"] = [gene_map.get(_feature_gene_id(feature)) for feature in diff["Feature"]]
    return diff


def _site_from_feature(feature: str) -> str:
    parts = str(feature).split("|")
    return parts[2] if len(parts) > 2 else ""


def _protein_log2fc_by_gene_from_diff(protein_diff: pd.DataFrame, protein_log2fc_column: str) -> pd.Series:
    if protein_log2fc_column not in protein_diff.columns:
        raise ValueError(f"Missing protein log2FC column: {protein_log2fc_column}")
    protein = protein_diff[pd.notna(protein_diff["Gene"])].copy()
    return protein.groupby("Gene")[protein_log2fc_column].mean()


def _selected_pathway_genes(
    *,
    pathway_term: str,
    direction: str,
    mode: str,
    terms: dict[str, set[str]],
    pure_up_genes: set[str],
    pure_down_genes: set[str],
) -> set[str]:
    if pathway_term not in terms:
        raise ValueError(f"Pathway term not found in GMT: {pathway_term}")
    pathway_genes = terms[pathway_term]
    if mode == "all_pathway":
        return set(pathway_genes)
    if mode != "filtered":
        raise ValueError(
            "gene_selection_mode must be one of: "
            f"{', '.join(sorted(PATHWAY_SELECTION_MODES))}."
        )
    if direction.lower() == "down":
        return pure_down_genes & pathway_genes
    return pure_up_genes & pathway_genes


def _build_phosphosite_protein_table(
    *,
    phospho_diff: pd.DataFrame,
    protein_log2fc: pd.Series,
    pathway_term: str,
    direction: str,
    selected_genes: set[str],
    phospho_log2fc_column: str,
) -> pd.DataFrame:
    if phospho_log2fc_column not in phospho_diff.columns:
        raise ValueError(f"Missing phosphosite log2FC column: {phospho_log2fc_column}")
    selected = phospho_diff[phospho_diff["Gene"].isin(selected_genes)].copy()
    selected["protein Log2FC(mean)"] = selected["Gene"].map(protein_log2fc)
    selected["phosphosite Log2FC(mean)"] = pd.to_numeric(selected[phospho_log2fc_column], errors="coerce")
    selected["phosphosite/protein"] = selected["phosphosite Log2FC(mean)"] / selected["protein Log2FC(mean)"]
    selected["pathway"] = pathway_term
    selected["selection_direction"] = direction
    selected["protein_matched"] = selected["protein Log2FC(mean)"].notna()

    optional_columns = [column for column in ["FDR", "Significance"] if column in selected.columns]
    table = selected.loc[
        :,
        [
            "Gene",
            "protein Log2FC(mean)",
            "phosphosite Log2FC(mean)",
            "Feature",
            "phosphosite/protein",
            "pathway",
            "selection_direction",
            "protein_matched",
            *optional_columns,
        ],
    ].rename(columns={"Feature": "Phosphosite"})
    return table.sort_values(["Gene", "Phosphosite"], kind="stable").reset_index(drop=True)


def _read_phosphosite_protein_table(path: str | Path) -> pd.DataFrame:
    table_path = Path(path)
    if not table_path.exists():
        raise FileNotFoundError(f"Filtered phosphosite-protein table does not exist: {table_path}")
    table = pd.read_csv(table_path, sep="\t")
    required = {"Gene", "protein Log2FC(mean)", "phosphosite Log2FC(mean)", "Phosphosite"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(
            f"Filtered phosphosite-protein table {table_path} is missing column(s): "
            f"{', '.join(sorted(missing))}"
        )
    if "protein_matched" not in table.columns:
        table["protein_matched"] = table["protein Log2FC(mean)"].notna()
    return table


def _points_from_phosphosite_protein_table(table: pd.DataFrame, highlight_sites: set[str]) -> pd.DataFrame:
    rows = []
    for _index, row in table.iterrows():
        gene = str(row["Gene"]).strip()
        feature = str(row["Phosphosite"]).strip()
        site = _site_from_feature(feature)
        gene_site = f"{gene}|{site}" if gene and site else ""
        highlight_label = _matched_highlight_label(gene_site, highlight_sites) if gene_site else None
        protein_value = pd.to_numeric(row["protein Log2FC(mean)"], errors="coerce")
        phospho_value = pd.to_numeric(row["phosphosite Log2FC(mean)"], errors="coerce")
        ratio = (
            float(phospho_value) / float(protein_value)
            if pd.notna(protein_value) and pd.notna(phospho_value) and float(protein_value) != 0
            else np.nan
        )
        rows.append(
            {
                "Feature": feature,
                "Gene": gene,
                "Site": site,
                "GeneSite": gene_site,
                "HighlightLabel": highlight_label if highlight_label else "",
                "protein_log2_tumor_over_nat": protein_value,
                "phosphosite_log2_tumor_over_nat": phospho_value,
                "phosphosite_to_protein_log2fc_ratio": ratio,
                "highlight": highlight_label is not None,
            }
        )
    return pd.DataFrame(rows)


def _finite_plot_points(points: pd.DataFrame) -> pd.DataFrame:
    x = pd.to_numeric(points["protein_log2_tumor_over_nat"], errors="coerce")
    y = pd.to_numeric(points["phosphosite_log2_tumor_over_nat"], errors="coerce")
    return points[np.isfinite(x) & np.isfinite(y)].copy()


def _table_pipeline_result_to_dict(
    rows: list[dict[str, Any]],
    config_path: Path,
    input_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    return {
        "output_dir": str(output_dir),
        "input_ini": str(input_path),
        "config_ini": str(config_path),
        "summary_tsv": str(output_dir / "pathway_scatter_summary.tsv"),
        "table_summary_tsv": str(output_dir / "phosphosite_protein_table_summary.tsv"),
        "pathway_count": len(rows),
        "records": json.dumps(rows),
    }


def run_phosphosite_protein_pathway_pipeline(
    input_ini: str | Path,
    *,
    output_dir: str | Path | None = None,
    config_name: str = DEFAULT_OUTPUT_CONFIG_NAME,
    template_path: str | Path | None = None,
    print_json: bool = True,
    **overrides: Any,
) -> dict[str, Any]:
    """
    Generate phosphosite-protein pathway tables and matching scatter plots.

    gene_selection_mode controls the gene universe:
    - filtered: use pure_up_genes for up pathways and pure_down_genes for down pathways.
    - all_pathway: use all genes in the GMT pathway that appear in the phospho diff table.

    Pathway direction defaults to up except pathways with "myogenesis" in the
    pathway term or folder name, which default to down. Override per pathway with
    a [directions] config section keyed by pathway folder name.
    """
    input_path = Path(input_ini)
    input_payload = _read_input_ini(input_path)
    pathway_output_dir = _resolve_pathway_output_dir(input_payload, output_dir)
    input_payload["output_dir"] = str(pathway_output_dir)
    input_payload.pop("out_dir", None)

    config_path = ensure_output_pathway_scatter_config(
        input_payload["output_dir"],
        config_name=config_name,
        template_path=template_path,
    )
    config_payload = load_pathway_scatter_replay_config(config_path)
    payload = _merge_payloads(config_payload, input_payload, overrides)

    table_paths = payload.get("table_paths", {})
    missing = [key for key in ("output_dir",) if not payload.get(key)]
    if not table_paths:
        missing.extend([key for key in ("protein_diff_path", "phospho_diff_path") if not payload.get(key)])
    if not payload.get("pathways"):
        missing.append("pathways or pathway_table_path")
    if not payload.get("highlight_paths"):
        missing.append("highlights or highlight_table_path")
    if table_paths:
        missing_tables = [folder_name for folder_name in payload.get("pathways", {}) if folder_name not in table_paths]
        if missing_tables:
            missing.append(f"tables for: {', '.join(missing_tables)}")
    if missing:
        raise ValueError(f"Missing required pathway pipeline input setting(s): {', '.join(missing)}")

    output_path = Path(payload["output_dir"])
    output_path.mkdir(parents=True, exist_ok=True)
    mode = str(payload.get("gene_selection_mode", "filtered")).strip().lower()
    if mode not in PATHWAY_SELECTION_MODES:
        raise ValueError(
            "gene_selection_mode must be one of: "
            f"{', '.join(sorted(PATHWAY_SELECTION_MODES))}."
        )

    if table_paths:
        terms: dict[str, set[str]] = {}
        protein_log2fc = pd.Series(dtype=float)
        phospho_diff = pd.DataFrame()
        pure_up_genes: set[str] = set()
        pure_down_genes: set[str] = set()
    else:
        gene_map = _read_gene_map(payload.get("fasta_path"))
        terms = _read_gmt_terms(payload.get("gmt_path"))
        protein_diff = _diff_with_gene(payload["protein_diff_path"], gene_map)
        phospho_diff = _diff_with_gene(payload["phospho_diff_path"], gene_map)
        protein_log2fc = _protein_log2fc_by_gene_from_diff(
            protein_diff,
            payload.get("protein_log2fc_column", "Log2FC(mean)"),
        )
        phospho_run_dir = Path(payload["phospho_run_dir"]) if payload.get("phospho_run_dir") else None
        pure_up_path = payload.get("pure_up_genes_path") or (phospho_run_dir / "pure_up_genes.tsv" if phospho_run_dir else None)
        pure_down_path = payload.get("pure_down_genes_path") or (phospho_run_dir / "pure_down_genes.tsv" if phospho_run_dir else None)
        pure_up_genes = _read_gene_list(pure_up_path) if pure_up_path else set()
        pure_down_genes = _read_gene_list(pure_down_path) if pure_down_path else set()

    directions = payload.get("directions", {})
    summary_rows: list[dict[str, Any]] = []
    table_summary_rows: list[dict[str, Any]] = []
    rc_updates = {
        "font.family": "sans-serif",
        "font.sans-serif": [payload.get("font_family", "Liberation Sans"), "Liberation Sans", "Arial", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
    with pd.ExcelWriter(output_path / "phosphosite_protein_pathway_tables.xlsx", engine="openpyxl") as writer:
        with matplotlib.rc_context(rc_updates):
            for folder_name, pathway_term in payload["pathways"].items():
                direction = str(directions.get(folder_name, "")).strip().lower()
                if not direction:
                    direction = "down" if "myogenesis" in f"{folder_name} {pathway_term}".lower() else "up"
                source_table_path = None
                if table_paths:
                    source_table_path = Path(table_paths[folder_name])
                    table = _read_phosphosite_protein_table(source_table_path)
                    selected_genes = set(table["Gene"].dropna().astype(str).str.strip())
                    if "selection_direction" in table.columns and table["selection_direction"].notna().any():
                        direction = str(table["selection_direction"].dropna().iloc[0]).strip().lower()
                else:
                    selected_genes = _selected_pathway_genes(
                        pathway_term=pathway_term,
                        direction=direction,
                        mode=mode,
                        terms=terms,
                        pure_up_genes=pure_up_genes,
                        pure_down_genes=pure_down_genes,
                    )
                    table = _build_phosphosite_protein_table(
                        phospho_diff=phospho_diff,
                        protein_log2fc=protein_log2fc,
                        pathway_term=pathway_term,
                        direction=direction,
                        selected_genes=selected_genes,
                        phospho_log2fc_column=payload.get("phospho_log2fc_column", "Log2FC(mean)"),
                    )

                pathway_dir = output_path / folder_name
                pathway_dir.mkdir(parents=True, exist_ok=True)
                safe = _safe_name(folder_name)
                table_tsv = source_table_path or (output_path / f"{safe}_phosphosite_protein.tsv")
                if source_table_path is None:
                    table.to_csv(table_tsv, sep="\t", index=False)
                table.to_excel(writer, sheet_name=f"{safe[:28]}", index=False)

                highlight_txt = Path(payload["highlight_paths"][folder_name])
                highlight_sites = _read_highlight_sites(highlight_txt)
                copied_highlight_txt = pathway_dir / "highlight_sites.txt"
                if highlight_txt.exists() and highlight_txt.resolve() != copied_highlight_txt.resolve():
                    shutil.copyfile(highlight_txt, copied_highlight_txt)
                else:
                    copied_highlight_txt = highlight_txt

                all_points = _points_from_phosphosite_protein_table(table, highlight_sites)
                plot_points = _finite_plot_points(all_points)
                all_points_tsv = pathway_dir / f"{safe}_filtered_table_points.tsv"
                points_tsv = pathway_dir / f"{safe}_scatter_points.tsv"
                highlight_tsv = pathway_dir / f"{safe}_highlight_points.tsv"
                missing_tsv = pathway_dir / f"{safe}_missing_highlights.tsv"
                png_path = pathway_dir / f"{safe}_protein_phosphosite_scatter.png"
                pdf_path = pathway_dir / f"{safe}_protein_phosphosite_scatter.pdf"
                tiff_path = pathway_dir / f"{safe}_protein_phosphosite_scatter.tiff"

                all_points.to_csv(all_points_tsv, sep="\t", index=False)
                plot_points.to_csv(points_tsv, sep="\t", index=False)
                plot_points[plot_points["highlight"]].to_csv(highlight_tsv, sep="\t", index=False)
                _missing_highlights(plot_points, highlight_sites).to_csv(missing_tsv, sep="\t", index=False)

                slope, intercept, r = _plot_points(
                    plot_points,
                    pathway_term,
                    png_path,
                    pdf_path,
                    tiff_path,
                    width=payload.get("width", 5.0),
                    height=payload.get("height", 4.0),
                    dpi=payload.get("dpi", 300),
                    tiff_dpi=payload.get("tiff_dpi", 600),
                    show_title=payload.get("show_title", True),
                    title_template=payload.get("title_template", "{pathway}"),
                    log2fc_label="mean",
                    font_scale=payload.get("font_scale", 1.0),
                    point_size=payload.get("point_size", 10.0),
                    highlight_point_size=payload.get("highlight_point_size", 18.0),
                )

                table_summary = {
                    "pathway": pathway_term,
                    "folder": folder_name,
                    "gene_selection_mode": mode,
                    "direction": direction,
                    "selected_gene_count": len(selected_genes),
                    "source_row_count": table.shape[0],
                    "protein_matched_gene_count": len(set(table.loc[table["protein_matched"], "Gene"])),
                    "protein_missing_gene_count": len(set(table.loc[~table["protein_matched"], "Gene"])) if table_paths else len(selected_genes - set(protein_log2fc.index)),
                    "protein_missing_genes": ";".join(sorted(set(table.loc[~table["protein_matched"], "Gene"]))) if table_paths else ";".join(sorted(selected_genes - set(protein_log2fc.index))),
                    "table_tsv": str(table_tsv),
                }
                table_summary_rows.append(table_summary)

                summary_rows.append(
                    {
                        **table_summary,
                        "output_dir": str(pathway_dir),
                        "plotted_point_count": plot_points.shape[0],
                        "dropped_nonfinite_count": all_points.shape[0] - plot_points.shape[0],
                        "highlight_count": int(plot_points["highlight"].sum()),
                        "trend_slope": slope,
                        "trend_intercept": intercept,
                        "trend_r": r,
                        "png": str(png_path),
                        "pdf": str(pdf_path),
                        "tiff": str(tiff_path),
                        "missing_highlights_tsv": str(missing_tsv),
                    }
                )
        pd.DataFrame(table_summary_rows).to_excel(writer, sheet_name="summary", index=False)

    pd.DataFrame(table_summary_rows).to_csv(output_path / "phosphosite_protein_table_summary.tsv", sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(output_path / "pathway_scatter_summary.tsv", sep="\t", index=False)

    result_dict = _table_pipeline_result_to_dict(summary_rows, config_path, input_path, output_path)
    if print_json:
        print(json.dumps(result_dict, indent=2), flush=True)
    return result_dict

