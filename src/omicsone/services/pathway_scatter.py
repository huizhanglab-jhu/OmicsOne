from __future__ import annotations

import configparser
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
try:
    from adjustText import adjust_text
except ImportError:  # pragma: no cover - optional plotting enhancement
    adjust_text = None

from omicsone.services.volcano_enrichment import DEFAULT_FASTA_PATH, _feature_gene_id, _read_gene_map


DEFAULT_HALLMARK_GMT = Path.home() / ".cache" / "gseapy" / "Enrichr.MSigDB_Hallmark_2020.gmt"


@dataclass(frozen=True)
class PathwayScatterResult:
    pathway: str
    output_dir: Path
    points_tsv: Path
    highlight_tsv: Path
    missing_highlights_tsv: Path
    highlight_txt: Path
    png: Path
    pdf: Path
    tiff: Path
    point_count: int
    highlight_count: int
    trend_slope: float | None
    trend_intercept: float | None
    trend_r: float | None


def _read_gmt_terms(gmt_path: str | Path) -> dict[str, set[str]]:
    terms: dict[str, set[str]] = {}
    with Path(gmt_path).open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term, _description, *genes = parts
            terms[term] = {gene for gene in genes if gene}
    return terms


def _read_highlight_sites(path: str | Path) -> set[str]:
    sites = set()
    with Path(path).open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            value = line.strip()
            if value and not value.startswith("#"):
                sites.add(value)
    return sites


def _site_key(feature_id: str, gene_map: dict[str, str]) -> tuple[str | None, str | None, str | None]:
    parts = str(feature_id).split("|")
    gene_id = _feature_gene_id(str(feature_id))
    gene = gene_map.get(gene_id)
    site = parts[2] if len(parts) > 2 else None
    label = f"{gene}|{site}" if gene and site else None
    return gene, site, label


def _normalize_site_label(value: str) -> str:
    gene, _, site = value.partition("|")
    return f"{gene.upper()}|{site.upper()}"


def _highlight_match(label: str, highlight_sites: set[str]) -> bool:
    return _matched_highlight_label(label, highlight_sites) is not None


def _matched_highlight_label(label: str, highlight_sites: set[str]) -> str | None:
    normalized = _normalize_site_label(label)
    normalized_highlights = {_normalize_site_label(site): site for site in highlight_sites}
    if normalized in normalized_highlights:
        return normalized_highlights[normalized]

    gene, _, site = normalized.partition("|")
    for highlight, original_highlight in normalized_highlights.items():
        highlight_gene, _, highlight_site = highlight.partition("|")
        if site == highlight_site and gene.startswith(highlight_gene):
            return original_highlight.upper()
    return None


def _protein_log2fc_by_gene(
    protein_diff_path: str | Path,
    gene_map: dict[str, str],
    log2fc_column: str,
) -> pd.Series:
    protein = pd.read_csv(protein_diff_path, sep="\t")
    if "Feature" not in protein.columns:
        protein = protein.rename(columns={protein.columns[0]: "Feature"})
    if log2fc_column not in protein.columns:
        raise ValueError(f"Missing protein log2FC column: {log2fc_column}")
    protein["Gene"] = [gene_map.get(_feature_gene_id(feature)) for feature in protein["Feature"]]
    protein = protein[pd.notna(protein["Gene"])].copy()
    return protein.groupby("Gene")[log2fc_column].mean()


def _phospho_points(
    phospho_diff_path: str | Path,
    protein_log2fc: pd.Series,
    pathway_genes: set[str],
    gene_map: dict[str, str],
    highlight_sites: set[str],
    log2fc_column: str,
) -> pd.DataFrame:
    phospho = pd.read_csv(phospho_diff_path, sep="\t")
    if "Feature" not in phospho.columns:
        phospho = phospho.rename(columns={phospho.columns[0]: "Feature"})
    if log2fc_column not in phospho.columns:
        raise ValueError(f"Missing phosphosite log2FC column: {log2fc_column}")

    rows = []
    for _index, row in phospho.iterrows():
        feature = row["Feature"]
        phospho_log2fc = row[log2fc_column]
        gene, site, label = _site_key(feature, gene_map)
        if not gene or not site or not label:
            continue
        if gene not in pathway_genes or gene not in protein_log2fc.index:
            continue
        protein_value = float(protein_log2fc.loc[gene])
        phospho_value = float(phospho_log2fc)
        ratio = phospho_value / protein_value if protein_value != 0 else np.nan
        highlight_label = _matched_highlight_label(label, highlight_sites)
        is_highlight = highlight_label is not None
        rows.append(
            {
                "Feature": feature,
                "Gene": gene,
                "Site": site,
                "GeneSite": label,
                "HighlightLabel": highlight_label if highlight_label else "",
                "protein_log2_tumor_over_nat": protein_value,
                "phosphosite_log2_tumor_over_nat": phospho_value,
                "phosphosite_to_protein_log2fc_ratio": ratio,
                "highlight": is_highlight,
            }
        )

    return pd.DataFrame(rows)


def _missing_highlights(points: pd.DataFrame, highlight_sites: set[str]) -> pd.DataFrame:
    matched = set(points.loc[points["highlight"], "GeneSite"]) if not points.empty else set()
    rows = []
    for site in sorted(highlight_sites):
        if not any(_highlight_match(match, {site}) for match in matched):
            rows.append({"requested_highlight": site})
    return pd.DataFrame(rows)


def _fit_trend(points: pd.DataFrame) -> tuple[float | None, float | None, float | None]:
    if points.shape[0] < 2:
        return None, None, None
    x = points["protein_log2_tumor_over_nat"].to_numpy(dtype=float)
    y = points["phosphosite_log2_tumor_over_nat"].to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 2:
        return None, None, None
    slope, intercept = np.polyfit(x[valid], y[valid], deg=1)
    r = float(np.corrcoef(x[valid], y[valid])[0, 1])
    return float(slope), float(intercept), r


def _format_site_label(value: str) -> str:
    gene, _, site = str(value).partition("|")
    return f"{gene} | {site.upper()}" if gene and site else str(value)


def _annotation_candidates(prefer_right: bool) -> list[tuple[int, int]]:
    right = [
        (22, 8),
        (28, 0),
        (22, -8),
        (38, 16),
        (38, -16),
        (56, 28),
        (56, -28),
        (76, 42),
        (76, -42),
        (96, 56),
        (96, -56),
        (118, 72),
        (118, -72),
        (140, 90),
        (140, -90),
    ]
    left = [
        (-22, 8),
        (-28, 0),
        (-22, -8),
        (-38, 16),
        (-38, -16),
        (-56, 28),
        (-56, -28),
        (-76, 42),
        (-76, -42),
        (-96, 56),
        (-96, -56),
        (-118, 72),
        (-118, -72),
        (-140, 90),
        (-140, -90),
    ]
    vertical = [
        (0, 28),
        (0, -28),
        (16, 36),
        (-16, 36),
        (16, -36),
        (-16, -36),
        (0, 54),
        (0, -54),
    ]
    return right + left + vertical if prefer_right else left + right + vertical


def _score_annotation_bbox(
    bbox: Any,
    *,
    axes_bbox: Any,
    point_pixels: np.ndarray,
    existing_bboxes: list[Any],
) -> float:
    padded = bbox.expanded(1.06, 1.18)
    outside = (
        padded.x0 < axes_bbox.x0
        or padded.x1 > axes_bbox.x1
        or padded.y0 < axes_bbox.y0
        or padded.y1 > axes_bbox.y1
    )
    score = 100.0 if outside else 0.0
    if point_pixels.size:
        inside_x = (point_pixels[:, 0] >= padded.x0) & (point_pixels[:, 0] <= padded.x1)
        inside_y = (point_pixels[:, 1] >= padded.y0) & (point_pixels[:, 1] <= padded.y1)
        score += float((inside_x & inside_y).sum()) * 20.0
    for existing in existing_bboxes:
        if padded.overlaps(existing.expanded(1.03, 1.12)):
            score += 100000.0
    return score


def _annotate_highlights(
    ax: plt.Axes,
    points: pd.DataFrame,
    highlights: pd.DataFrame,
    *,
    font_scale: float,
) -> None:
    if highlights.empty:
        return

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    x_span = max(float(x_max - x_min), 1.0)
    y_span = max(float(y_max - y_min), 1.0)
    left_x = x_min + 0.08 * x_span
    right_x = x_max - 0.08 * x_span
    lower = y_min + 0.08 * y_span
    upper = y_max - 0.08 * y_span
    min_sep = 0.115 * y_span * min(font_scale, 1.35)
    point_xy = points[["protein_log2_tumor_over_nat", "phosphosite_log2_tumor_over_nat"]].to_numpy(dtype=float)
    point_xy = point_xy[np.isfinite(point_xy).all(axis=1)]

    def side_score(row: pd.Series, side: str) -> float:
        x = float(row["protein_log2_tumor_over_nat"])
        y = float(row["phosphosite_log2_tumor_over_nat"])
        if side == "left":
            x0, x1 = left_x, left_x + 0.30 * x_span
            col_x = left_x
        else:
            x0, x1 = right_x - 0.30 * x_span, right_x
            col_x = right_x
        y0, y1 = y - 0.055 * y_span, y + 0.055 * y_span
        local_points = 0
        if point_xy.size:
            local_points = int(((point_xy[:, 0] >= x0) & (point_xy[:, 0] <= x1) & (point_xy[:, 1] >= y0) & (point_xy[:, 1] <= y1)).sum())
        distance_penalty = abs(col_x - x) / x_span
        return local_points * 8.0 + distance_penalty

    placed: dict[int, tuple[str, float]] = {}
    groups: dict[str, list[tuple[int, float]]] = {"left": [], "right": []}
    for row_index, row in highlights.iterrows():
        side = "left" if side_score(row, "left") <= side_score(row, "right") else "right"
        y = float(row["phosphosite_log2_tumor_over_nat"])
        groups[side].append((row_index, min(max(y, lower), upper)))

    for side, entries in groups.items():
        entries = sorted(entries, key=lambda item: item[1])
        if not entries:
            continue
        spread = [y for _idx, y in entries]
        for idx in range(1, len(spread)):
            spread[idx] = max(spread[idx], spread[idx - 1] + min_sep)
        overflow = spread[-1] - upper
        if overflow > 0:
            spread = [y - overflow for y in spread]
        for idx in range(1, len(spread)):
            spread[idx] = max(spread[idx], spread[idx - 1] + min_sep)
        if spread[0] < lower:
            shift = lower - spread[0]
            spread = [y + shift for y in spread]
        for (row_index, _original_y), adjusted_y in zip(entries, spread):
            placed[row_index] = (side, adjusted_y)

    for row_index, row in highlights.iterrows():
        x = float(row["protein_log2_tumor_over_nat"])
        y = float(row["phosphosite_log2_tumor_over_nat"])
        label = row.get("HighlightLabel") or row["GeneSite"]
        side, text_y = placed[row_index]
        text_x = left_x if side == "left" else right_x
        ax.annotate(
            _format_site_label(label),
            xy=(x, y),
            xytext=(text_x, text_y),
            fontsize=6.5 * font_scale,
            fontweight="bold",
            color="black",
            ha="left" if side == "left" else "right",
            va="center",
            arrowprops={"arrowstyle": "-", "color": "#f3c38b", "linewidth": 1.0, "shrinkA": 0, "shrinkB": 3},
            zorder=4,
        )
    return

    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    axes_bbox = ax.bbox.padded(-2)
    point_pixels = ax.transData.transform(
        points[["protein_log2_tumor_over_nat", "phosphosite_log2_tumor_over_nat"]].to_numpy(dtype=float)
    )
    existing_bboxes: list[Any] = []

    for _row_index, row in highlights.iterrows():
        x = float(row["protein_log2_tumor_over_nat"])
        y = float(row["phosphosite_log2_tumor_over_nat"])
        label = row.get("HighlightLabel") or row["GeneSite"]
        text = _format_site_label(label)
        prefer_right = x < 0.2
        best: tuple[float, tuple[int, int], Any] | None = None
        for dx, dy in _annotation_candidates(prefer_right):
            probe = ax.annotate(
                text,
                xy=(x, y),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=6.5,
                fontweight="bold",
                color="black",
                ha="left" if dx >= 0 else "right",
                va="center",
                arrowprops={"arrowstyle": "-", "color": "#f3c38b", "linewidth": 1.0, "shrinkA": 0, "shrinkB": 3},
                zorder=4,
            )
            ax.figure.canvas.draw()
            bbox = probe.get_window_extent(renderer)
            score = _score_annotation_bbox(
                bbox,
                axes_bbox=axes_bbox,
                point_pixels=point_pixels,
                existing_bboxes=existing_bboxes,
            )
            probe.remove()
            if best is None or score < best[0]:
                best = (score, (dx, dy), bbox)
        if best is None:
            continue

        _score, (best_dx, best_dy), _bbox = best
        annotation = ax.annotate(
            text,
            xy=(x, y),
            xytext=(best_dx, best_dy),
            textcoords="offset points",
            fontsize=6.5,
            fontweight="bold",
            color="black",
            ha="left" if best_dx >= 0 else "right",
            va="center",
            arrowprops={"arrowstyle": "-", "color": "#f3c38b", "linewidth": 1.0, "shrinkA": 0, "shrinkB": 3},
            zorder=4,
        )
        ax.figure.canvas.draw()
        existing_bboxes.append(annotation.get_window_extent(renderer))


def _plot_points(
    points: pd.DataFrame,
    pathway: str,
    png_path: Path,
    pdf_path: Path,
    tiff_path: Path,
    *,
    width: float,
    height: float,
    dpi: int,
    tiff_dpi: int,
    show_title: bool,
    title_template: str,
    log2fc_label: str,
    font_scale: float,
    point_size: float,
    highlight_point_size: float,
) -> tuple[float | None, float | None, float | None]:
    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    background = points[~points["highlight"]]
    highlights = points[points["highlight"]]

    ax.scatter(
        background["protein_log2_tumor_over_nat"],
        background["phosphosite_log2_tumor_over_nat"],
        s=point_size,
        color="#1f77b4",
        alpha=0.9,
    )
    ax.scatter(
        highlights["protein_log2_tumor_over_nat"],
        highlights["phosphosite_log2_tumor_over_nat"],
        s=highlight_point_size,
        color="red",
        edgecolor="red",
        linewidth=0.2,
        zorder=3,
    )
    ax.axvline(0, color="#cfcfcf", linewidth=1.0)
    ax.axhline(0, color="#cfcfcf", linewidth=1.0)

    slope, intercept, r = _fit_trend(points)
    if slope is not None and intercept is not None:
        x_min = float(points["protein_log2_tumor_over_nat"].min())
        x_max = float(points["protein_log2_tumor_over_nat"].max())
        xs = np.linspace(x_min, x_max, 100)
        ax.plot(xs, slope * xs + intercept, color="#cfcfcf", linestyle="--", linewidth=1.2)

    ax.set_xlabel(f"Protein {log2fc_label} log2(Tumor/NATs)", fontsize=9 * font_scale)
    ax.set_ylabel(f"Phosphosite {log2fc_label} log2(Tumor/NATs)", fontsize=9 * font_scale)
    if show_title:
        ax.set_title(title_template.format(pathway=pathway), fontsize=11 * font_scale, color="#6f6f6f")
    ax.tick_params(axis="both", labelsize=8 * font_scale)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.margins(x=0.18, y=0.15)
    _annotate_highlights(ax, points, highlights, font_scale=font_scale)
    fig.tight_layout(pad=2.0)
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", pad_inches=0.2)
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.2)
    fig.savefig(tiff_path, dpi=tiff_dpi, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)
    return slope, intercept, r


def run_pathway_scatter_analysis(
    *,
    protein_diff_path: str,
    phospho_diff_path: str,
    output_dir: str,
    pathways: dict[str, str],
    highlight_paths: dict[str, str],
    fasta_path: str = DEFAULT_FASTA_PATH,
    gmt_path: str = str(DEFAULT_HALLMARK_GMT),
    width: float = 5.0,
    height: float = 4.0,
    dpi: int = 300,
    tiff_dpi: int = 600,
    show_title: bool = True,
    title_template: str = "{pathway}",
    font_scale: float = 1.0,
    point_size: float = 10.0,
    highlight_point_size: float = 18.0,
    font_family: str = "Liberation Sans",
    protein_log2fc_column: str = "Log2FC(median)",
    phospho_log2fc_column: str = "Log2FC(median)",
) -> list[PathwayScatterResult]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    gene_map = _read_gene_map(fasta_path)
    terms = _read_gmt_terms(gmt_path)
    protein_log2fc = _protein_log2fc_by_gene(
        protein_diff_path,
        gene_map,
        protein_log2fc_column,
    )

    results = []
    rc_updates: dict[str, Any] = {
        "font.family": "sans-serif",
        "font.sans-serif": [font_family, "Liberation Sans", "Arial", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
    with matplotlib.rc_context(rc_updates):
        for folder_name, pathway_term in pathways.items():
            if pathway_term not in terms:
                raise ValueError(f"Pathway term not found in GMT: {pathway_term}")
            if folder_name not in highlight_paths:
                raise ValueError(f"Missing highlight path for pathway folder: {folder_name}")

            pathway_dir = output_path / folder_name
            pathway_dir.mkdir(parents=True, exist_ok=True)
            highlight_txt = Path(highlight_paths[folder_name])
            highlight_sites = _read_highlight_sites(highlight_txt)
            points = _phospho_points(
                phospho_diff_path=phospho_diff_path,
                protein_log2fc=protein_log2fc,
                pathway_genes=terms[pathway_term],
                gene_map=gene_map,
                highlight_sites=highlight_sites,
                log2fc_column=phospho_log2fc_column,
            )

            safe_name = folder_name.replace(" ", "_").replace("-", "_").lower()
            points_tsv = pathway_dir / f"{safe_name}_scatter_points.tsv"
            highlight_tsv = pathway_dir / f"{safe_name}_highlight_points.tsv"
            missing_highlights_tsv = pathway_dir / f"{safe_name}_missing_highlights.tsv"
            png_path = pathway_dir / f"{safe_name}_protein_phosphosite_scatter.png"
            pdf_path = pathway_dir / f"{safe_name}_protein_phosphosite_scatter.pdf"
            tiff_path = pathway_dir / f"{safe_name}_protein_phosphosite_scatter.tiff"

            points.to_csv(points_tsv, sep="\t", index=False)
            points[points["highlight"]].to_csv(highlight_tsv, sep="\t", index=False)
            _missing_highlights(points, highlight_sites).to_csv(
                missing_highlights_tsv,
                sep="\t",
                index=False,
            )
            slope, intercept, r = _plot_points(
                points,
                pathway_term,
                png_path,
                pdf_path,
                tiff_path,
                width=width,
                height=height,
                dpi=dpi,
                tiff_dpi=tiff_dpi,
                show_title=show_title,
                title_template=title_template,
                log2fc_label="mean" if "mean" in phospho_log2fc_column.lower() else "median",
                font_scale=font_scale,
                point_size=point_size,
                highlight_point_size=highlight_point_size,
            )

            results.append(
                PathwayScatterResult(
                    pathway=pathway_term,
                    output_dir=pathway_dir,
                    points_tsv=points_tsv,
                    highlight_tsv=highlight_tsv,
                    missing_highlights_tsv=missing_highlights_tsv,
                    highlight_txt=highlight_txt,
                    png=png_path,
                    pdf=pdf_path,
                    tiff=tiff_path,
                    point_count=points.shape[0],
                    highlight_count=int(points["highlight"].sum()),
                    trend_slope=slope,
                    trend_intercept=intercept,
                    trend_r=r,
                )
            )

    summary = pd.DataFrame(
        [
            {
                "pathway": result.pathway,
                "output_dir": str(result.output_dir),
                "point_count": result.point_count,
                "highlight_count": result.highlight_count,
                "trend_slope": result.trend_slope,
                "trend_intercept": result.trend_intercept,
                "trend_r": result.trend_r,
                "protein_log2fc_column": protein_log2fc_column,
                "phospho_log2fc_column": phospho_log2fc_column,
                "png": str(result.png),
                "pdf": str(result.pdf),
                "tiff": str(result.tiff),
                "missing_highlights_tsv": str(result.missing_highlights_tsv),
            }
            for result in results
        ]
    )
    summary.to_csv(output_path / "pathway_scatter_summary.tsv", sep="\t", index=False)
    return results


def _split_lines(value: str) -> list[str]:
    return [item.strip() for item in value.replace("\n", ",").split(",") if item.strip()]


def load_pathway_scatter_config(config_path: str | Path) -> dict[str, Any]:
    parser = configparser.ConfigParser()
    parser.optionxform = str
    parser.read(config_path, encoding="utf-8-sig")

    pathways = dict(parser.items("pathways"))
    highlight_paths = dict(parser.items("highlights"))
    payload: dict[str, Any] = {
        "protein_diff_path": parser.get("input", "protein_diff_path"),
        "phospho_diff_path": parser.get("input", "phospho_diff_path"),
        "output_dir": parser.get("output", "output_dir"),
        "pathways": pathways,
        "highlight_paths": highlight_paths,
    }
    if parser.has_option("input", "fasta_path"):
        payload["fasta_path"] = parser.get("input", "fasta_path")
    if parser.has_option("input", "gmt_path"):
        payload["gmt_path"] = parser.get("input", "gmt_path")
    if parser.has_option("input", "protein_log2fc_column"):
        payload["protein_log2fc_column"] = parser.get("input", "protein_log2fc_column")
    if parser.has_option("input", "phospho_log2fc_column"):
        payload["phospho_log2fc_column"] = parser.get("input", "phospho_log2fc_column")
    if parser.has_section("plot"):
        for key in ["width", "height"]:
            if parser.has_option("plot", key):
                payload[key] = parser.getfloat("plot", key)
        for key in ["dpi", "tiff_dpi"]:
            if parser.has_option("plot", key):
                payload[key] = parser.getint("plot", key)
        if parser.has_option("plot", "show_title"):
            payload["show_title"] = parser.getboolean("plot", "show_title")
        if parser.has_option("plot", "title_template"):
            payload["title_template"] = parser.get("plot", "title_template")
        for key in ["font_scale", "point_size", "highlight_point_size"]:
            if parser.has_option("plot", key):
                payload[key] = parser.getfloat("plot", key)
        if parser.has_option("plot", "font_family"):
            payload["font_family"] = parser.get("plot", "font_family")
    return payload


def run_pathway_scatter_analysis_from_config(config_path: str | Path) -> list[PathwayScatterResult]:
    return run_pathway_scatter_analysis(**load_pathway_scatter_config(config_path))
