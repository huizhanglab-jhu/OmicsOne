from __future__ import annotations

import argparse
import configparser
import json
import shutil
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

from omicsone.plots.cnv_correlation import (
    load_chromosome_layout,
    plot_chromosome_panel,
    plot_corr_heatmap_panel,
    plot_gistic_panel,
    plot_local_distal_distribution_panel,
    plot_self_cis_distribution_panel,
    save_combined_three_class_split_self_figure,
    summarize_cnv_distribution_three_class,
)


def run_cnv_correlation_clean_figures(
    config_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    print_json: bool = True,
) -> dict[str, Any]:
    config_file = Path(config_path)
    config = _read_config(config_file)

    output_path = Path(
        output_dir
        or config.get("output", "clean_dir", fallback="")
        or config.get("output", "figures_dir", fallback="")
    )
    if not output_path:
        raise ValueError("Config must define [output] clean_dir, or pass --output-dir.")
    output_path.mkdir(parents=True, exist_ok=True)

    prefix = config.get("output", "prefix", fallback="cnv_correlation")
    target_type = config.get("output", "target_type", fallback="protein")
    use_simple_names = config.getboolean("output", "use_simple_names", fallback=True)
    figure_dpi = config.getint("plot", "figure_dpi", fallback=150)
    tiff_dpi = config.getint("plot", "tiff_dpi", fallback=600)
    heatmap_vmin = config.getfloat("plot", "heatmap_vmin", fallback=-0.5)
    heatmap_vmax = config.getfloat("plot", "heatmap_vmax", fallback=0.5)
    mpl.rcParams["font.family"] = config.get("plot", "font_family", fallback="Liberation Sans")

    annotated_path = _resolve_input_file(Path(config["inputs"]["annotated_correlations"]))
    count_path = _resolve_input_file(Path(config["inputs"]["cnv_distribution_counts"]))
    gistic_path = _resolve_input_file(Path(config["inputs"]["gistic_counts"]))
    chromosomes_file = Path(config["inputs"]["chromosomes_file"])

    annotated = pd.read_csv(annotated_path, sep="\t", low_memory=False)
    old_counts = pd.read_csv(count_path, sep="\t")
    gistic_counts = pd.read_csv(gistic_path, sep="\t")
    layout = load_chromosome_layout(chromosomes_file)
    arm_bounds = old_counts[
        ["cnv.arm", "cnv.arm.start", "cnv.arm.end"]
    ].rename(columns={"cnv.arm": "arm"})
    count_df = summarize_cnv_distribution_three_class(
        annotated,
        arm_bounds=arm_bounds,
        target_type=target_type,
    )
    amp_segments, deletion_segments = _gistic_segments(gistic_counts)

    image_outputs: list[Path] = []
    data_outputs: list[Path] = []

    def stem(name: str) -> str:
        if use_simple_names:
            return name
        legacy_names = {
            "chromosome": f"{prefix}_chromosome",
            "gistic": f"{prefix}_gistic_3class_split_self",
            "corr_heatmap": f"{prefix}_corr_heatmap_3class_split_self_4x3",
            "self_cis": f"{prefix}_self_cis_distribution_3class_split_self",
            "local_distal": f"{prefix}_local_distal_distribution_3class_split_self",
            "combined": f"{prefix}_combined_3class_split_self",
            "counts": f"{prefix}_cnv_distribution_3class_counts",
            "summary": f"{prefix}_3class_summary",
        }
        return legacy_names[name]

    _save_table(count_df, output_path / f"{stem('counts')}.tsv", data_outputs)
    annotated_output = output_path / f"{stem('corr_heatmap')}.tsv"
    shutil.copy2(annotated_path, annotated_output)
    data_outputs.append(annotated_output)

    chrom_df = layout.chromosomes.copy()
    chrom_df["End"] = chrom_df["Start"] + chrom_df["Total length (bp)"]
    chrom_df["Tick"] = layout.ticks
    chrom_df["Display label"] = [
        str(label).replace("\n", "") for label in layout.labels
    ]
    _save_table(chrom_df, output_path / f"{stem('chromosome')}.tsv", data_outputs)
    with plt.style.context("seaborn-v0_8-white"):
        fig = plt.figure(figsize=(8, 0.9), dpi=figure_dpi)
        ax = fig.add_axes([0.05, 0.42, 0.9, 0.23])
        plot_chromosome_panel(ax, layout=layout, stagger_dense_labels=True)
        _save_all(
            fig,
            output_path / stem("chromosome"),
            image_outputs,
            tiff_dpi=tiff_dpi,
        )

    _save_table(gistic_counts, output_path / f"{stem('gistic')}.tsv", data_outputs)
    with plt.style.context("seaborn-v0_8-white"):
        fig, ax = plt.subplots(figsize=(8, 2.5), dpi=figure_dpi)
        plot_gistic_panel(ax, amp_segments, deletion_segments, layout=layout)
        handles = [
            Rectangle((0, 0), 1, 1, color="red"),
            Rectangle((0, 0), 1, 1, color="blue"),
        ]
        ax.legend(
            handles,
            ["Amplification", "Deletion"],
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=7,
        )
        _save_all(
            fig,
            output_path / stem("gistic"),
            image_outputs,
            tiff_dpi=tiff_dpi,
        )

    with plt.style.context("seaborn-v0_8-white"):
        fig = plt.figure(figsize=(8, 6), dpi=figure_dpi)
        ax = fig.add_axes([0.08, 0.11, 0.76, 0.82])
        cax = fig.add_axes([0.88, 0.11, 0.025, 0.82])
        plot_corr_heatmap_panel(
            ax,
            annotated,
            layout=layout,
            target_type=target_type,
            vmin=heatmap_vmin,
            vmax=heatmap_vmax,
        )
        mappable = mpl.cm.ScalarMappable(
            norm=mpl.colors.Normalize(vmin=heatmap_vmin, vmax=heatmap_vmax),
            cmap="bwr",
        )
        mappable.set_array([])
        colorbar = fig.colorbar(mappable, cax=cax)
        colorbar.set_label("Spearman r", fontsize=7)
        colorbar.ax.tick_params(labelsize=7)
        _save_all(
            fig,
            output_path / stem("corr_heatmap"),
            image_outputs,
            tight=False,
            tiff_dpi=tiff_dpi,
        )

    self_df = count_df[
        ["cnv.arm", "pos.self_cis.count", "cnv.arm.start", "cnv.arm.end"]
    ].rename(columns={"pos.self_cis.count": "positive_self_cis_count"})
    _save_table(self_df, output_path / f"{stem('self_cis')}.tsv", data_outputs)
    with plt.style.context("seaborn-v0_8-white"):
        fig, ax = plt.subplots(figsize=(8, 2.0), dpi=figure_dpi)
        plot_self_cis_distribution_panel(ax, count_df, layout=layout)
        _save_all(
            fig,
            output_path / stem("self_cis"),
            image_outputs,
            tiff_dpi=tiff_dpi,
        )

    local_distal_df = count_df[
        [
            "cnv.arm",
            "pos.same_arm_local.count",
            "pos.distal_trans.count",
            "neg.same_arm_local.count",
            "neg.distal_trans.count",
            "cnv.arm.start",
            "cnv.arm.end",
        ]
    ].copy()
    _save_table(
        local_distal_df,
        output_path / f"{stem('local_distal')}.tsv",
        data_outputs,
    )
    with plt.style.context("seaborn-v0_8-white"):
        fig, ax = plt.subplots(figsize=(8, 2.5), dpi=figure_dpi)
        plot_local_distal_distribution_panel(ax, count_df, layout=layout)
        _save_all(
            fig,
            output_path / stem("local_distal"),
            image_outputs,
            tiff_dpi=tiff_dpi,
        )

    combined = save_combined_three_class_split_self_figure(
        output_path / f"{stem('combined')}.png",
        annotated,
        count_df,
        layout=layout,
        target_type=target_type,
        gistic_segments=(amp_segments, deletion_segments),
        dpi=figure_dpi,
        heatmap_vmin=heatmap_vmin,
        heatmap_vmax=heatmap_vmax,
    )
    image_outputs.extend(
        [
            combined.with_suffix(".png"),
            combined.with_suffix(".pdf"),
            combined.with_suffix(".tiff"),
        ],
    )

    manifest = pd.DataFrame(
        [
            _manifest_row("chromosome", stem("chromosome"), "Chromosome positions"),
            _manifest_row("gistic", stem("gistic"), "GISTIC amplification/deletion"),
            _manifest_row("corr_heatmap", stem("corr_heatmap"), "Correlation heatmap"),
            _manifest_row("self_cis", stem("self_cis"), "Positive self-cis counts"),
            _manifest_row("local_distal", stem("local_distal"), "Local/distal counts"),
            _manifest_row(
                "combined",
                stem("combined"),
                "Combined figure",
                data_table=f"{stem('combined')}_manifest.tsv",
            ),
        ],
    )
    _save_table(
        manifest,
        output_path / f"{stem('combined')}_manifest.tsv",
        data_outputs,
    )
    summary = count_df[
        [
            "pos.self_cis.count",
            "neg.self_cis.count",
            "pos.same_arm_local.count",
            "pos.distal_trans.count",
            "neg.same_arm_local.count",
            "neg.distal_trans.count",
            "count.sum",
        ]
    ].sum().astype(int)
    summary_df = summary.reset_index()
    summary_df.columns = ["metric", "count"]
    _save_table(summary_df, output_path / f"{stem('summary')}.tsv", data_outputs)

    result = {
        "config": str(config_file),
        "output_dir": str(output_path),
        "image_outputs": [str(path) for path in image_outputs],
        "data_outputs": [str(path) for path in data_outputs],
    }
    if print_json:
        print(json.dumps(result, indent=2), flush=True)
    return result


def _read_config(path: Path) -> configparser.ConfigParser:
    if not path.is_file():
        raise FileNotFoundError(f"Config file does not exist: {path}")
    config = configparser.ConfigParser()
    config.read(path, encoding="utf-8-sig")
    return config


def _resolve_input_file(path: Path) -> Path:
    if path.is_file():
        return path
    history_path = path.parent / "history" / path.name
    if history_path.is_file():
        return history_path
    raise FileNotFoundError(f"Input file does not exist: {path}")


def _save_table(df: pd.DataFrame, path: Path, outputs: list[Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    outputs.append(path)


def _save_all(
    fig: plt.Figure,
    base: Path,
    outputs: list[Path],
    *,
    tight: bool = True,
    tiff_dpi: int = 600,
) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    kwargs = {"bbox_inches": "tight"} if tight else {}
    png = base.with_suffix(".png")
    pdf = base.with_suffix(".pdf")
    tiff = base.with_suffix(".tiff")
    fig.savefig(png, **kwargs)
    with mpl.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42}):
        fig.savefig(pdf, **kwargs)
    fig.savefig(tiff, dpi=tiff_dpi, **kwargs)
    plt.close(fig)
    outputs.extend([png, pdf, tiff])


def _gistic_segments(
    table: pd.DataFrame,
) -> tuple[
    list[tuple[tuple[float, float], tuple[float, float]]],
    list[tuple[tuple[float, float], tuple[float, float]]],
]:
    amp_segments = [
        (
            (float(row.location), 0.0),
            (float(row.location), float(row.amplification_count)),
        )
        for row in table.itertuples(index=False)
        if float(row.amplification_count) != 0.0
    ]
    deletion_segments = [
        (
            (float(row.location), 0.0),
            (float(row.location), -float(row.deletion_count)),
        )
        for row in table.itertuples(index=False)
        if float(row.deletion_count) != 0.0
    ]
    return amp_segments, deletion_segments


def _manifest_row(
    panel: str,
    image_stem: str,
    description: str,
    *,
    data_table: str | None = None,
) -> dict[str, str]:
    return {
        "panel": panel,
        "image_stem": image_stem,
        "data_table": data_table or f"{image_stem}.tsv",
        "description": description,
    }


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate standardized clean CNV correlation figures.",
    )
    parser.add_argument("--config", required=True, help="Path to clean figure config.ini.")
    parser.add_argument("--output-dir", help="Override [output] clean_dir.")
    parser.add_argument("--quiet", action="store_true", help="Do not print JSON output.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    run_cnv_correlation_clean_figures(
        args.config,
        output_dir=args.output_dir,
        print_json=not args.quiet,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
