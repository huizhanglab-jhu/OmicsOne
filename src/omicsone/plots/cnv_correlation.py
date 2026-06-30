from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle


mpl.rcParams["font.family"] = "Liberation Sans"
mpl.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "DejaVu Sans", *mpl.rcParams["font.sans-serif"]]

TargetType = Literal["protein", "rna"]
BP_PLOT_SCALE = 1_000_000.0


@dataclass(frozen=True)
class ChromosomeLayout:
    chromosomes: pd.DataFrame
    total_bp: int
    start_map: dict[str, int]
    bounds: list[int]
    ticks: list[int]
    labels: list[str]
    y_labels: list[str]


@dataclass(frozen=True)
class CnvCorrelationFigureResult:
    output_dir: Path
    annotated_correlations_file: Path
    cnv_distribution_counts_file: Path
    corr_heatmap_file: Path
    cnv_distribution_file: Path
    chromosome_file: Path
    gistic_file: Path | None
    gistic_counts_file: Path | None
    combined_file: Path
    filtered_correlation_count: int
    annotated_correlation_count: int


def generate_cnv_correlation_figures(
    correlation_file: str | Path,
    cnv_for_corr_file: str | Path,
    target_for_corr_file: str | Path,
    fasta_file: str | Path,
    chromosomes_file: str | Path,
    cytoband_file: str | Path,
    output_dir: str | Path,
    *,
    gistic_file: str | Path | None = None,
    target_type: TargetType = "protein",
    correlation_threshold: float = 0.5,
    output_prefix: str | None = None,
    chunksize: int = 1_000_000,
    dpi: int = 300,
) -> CnvCorrelationFigureResult:
    """Generate CNV correlation panels from Spearman and *_for_corr files.

    The implementation follows the HNSCC protein-CNV notebook workflow, but
    keeps each step reusable and avoids the slow notebook row loops.
    """

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    prefix = output_prefix or f"cnv_vs_{target_type}"
    layout = load_chromosome_layout(chromosomes_file)
    cytobands = load_cytobands(cytoband_file)
    arm_bounds = build_cytoband_arm_bounds(cytobands, layout.start_map)
    gene_annotation = load_gene_annotation(fasta_file, layout, cytobands)

    corr = read_correlation_subset(
        correlation_file,
        target_type=target_type,
        threshold=correlation_threshold,
        chunksize=chunksize,
    )
    filtered_count = int(corr.shape[0])

    cnv_genes = read_matrix_gene_index(cnv_for_corr_file)
    cnv_samples = read_matrix_sample_columns(cnv_for_corr_file)
    target_genes = read_matrix_gene_index(target_for_corr_file)
    annotated = annotate_correlation_genes(
        corr,
        cnv_genes=cnv_genes,
        target_genes=target_genes,
        gene_annotation=gene_annotation,
        target_type=target_type,
    )

    count_df = summarize_cnv_distribution(
        annotated,
        arm_bounds=arm_bounds,
        target_type=target_type,
    )

    annotated_file = output_path / f"{prefix}_annotated_correlations.tsv"
    count_file = output_path / f"{prefix}_cnv_distribution_counts.tsv"
    annotated.to_csv(annotated_file, sep="\t", index=False)
    count_df.to_csv(count_file, sep="\t", index=False)

    heatmap_file = output_path / f"{prefix}_corr_heatmap.png"
    distribution_file = output_path / f"{prefix}_cnv_distribution.png"
    chromosome_file = output_path / f"{prefix}_chromosome.png"
    gistic_panel_file = output_path / f"{prefix}_gistic.png"
    gistic_counts_file = output_path / f"{prefix}_gistic_counts.tsv"
    combined_file = output_path / f"{prefix}_combined.png"

    gistic_segments: tuple[list[tuple[tuple[float, float], tuple[float, float]]], list[tuple[tuple[float, float], tuple[float, float]]]] | None = None
    if gistic_file is not None:
        gistic_summary, amp_segments, deletion_segments = summarize_gistic_segments(
            gistic_file,
            genes=cnv_genes,
            samples=cnv_samples,
            gene_annotation=gene_annotation,
        )
        gistic_summary.to_csv(gistic_counts_file, sep="\t", index=False)
        gistic_segments = (amp_segments, deletion_segments)
        save_gistic_panel(
            gistic_panel_file,
            amp_segments,
            deletion_segments,
            layout=layout,
            dpi=dpi,
        )

    save_corr_heatmap_panel(
        heatmap_file,
        annotated,
        layout=layout,
        target_type=target_type,
        dpi=dpi,
    )
    save_cnv_distribution_panel(
        distribution_file,
        count_df,
        layout=layout,
        dpi=dpi,
    )
    save_chromosome_panel(chromosome_file, layout=layout, dpi=dpi)
    save_combined_figure(
        combined_file,
        annotated,
        count_df,
        layout=layout,
        target_type=target_type,
        gistic_segments=gistic_segments,
        dpi=dpi,
    )

    return CnvCorrelationFigureResult(
        output_dir=output_path,
        annotated_correlations_file=annotated_file,
        cnv_distribution_counts_file=count_file,
        corr_heatmap_file=heatmap_file,
        cnv_distribution_file=distribution_file,
        chromosome_file=chromosome_file,
        gistic_file=gistic_panel_file if gistic_file is not None else None,
        gistic_counts_file=gistic_counts_file if gistic_file is not None else None,
        combined_file=combined_file,
        filtered_correlation_count=filtered_count,
        annotated_correlation_count=int(annotated.shape[0]),
    )


def read_correlation_subset(
    correlation_file: str | Path,
    *,
    target_type: TargetType,
    threshold: float = 0.5,
    chunksize: int = 1_000_000,
) -> pd.DataFrame:
    path = Path(correlation_file)
    if not path.is_file():
        raise FileNotFoundError(f"Correlation file does not exist: {path}")

    target_index_col = f"{target_type}.index"
    columns = ["cnv.index", target_index_col, "correlation"]
    parts: list[pd.DataFrame] = []

    reader = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=columns,
        chunksize=chunksize,
    )
    for chunk in reader:
        chunk["correlation"] = pd.to_numeric(chunk["correlation"], errors="coerce")
        filtered = chunk[chunk["correlation"].abs() > threshold]
        if not filtered.empty:
            parts.append(filtered.copy())

    if not parts:
        return pd.DataFrame(columns=columns)

    return pd.concat(parts, ignore_index=True)


def read_matrix_gene_index(matrix_file: str | Path) -> np.ndarray:
    path = Path(matrix_file)
    if not path.is_file():
        raise FileNotFoundError(f"Matrix file does not exist: {path}")

    index_df = pd.read_csv(path, sep="\t", usecols=[0], dtype=str)
    return index_df.iloc[:, 0].str.split(".").str[0].to_numpy(dtype=object)


def read_matrix_sample_columns(matrix_file: str | Path) -> list[str]:
    path = Path(matrix_file)
    if not path.is_file():
        raise FileNotFoundError(f"Matrix file does not exist: {path}")

    columns = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    return [str(column) for column in columns[1:]]


def load_chromosome_layout(chromosomes_file: str | Path) -> ChromosomeLayout:
    path = Path(chromosomes_file)
    if not path.is_file():
        raise FileNotFoundError(f"Chromosome file does not exist: {path}")

    chromosomes = pd.read_excel(path).copy()
    required = {"Chromosome", "Total length (bp)"}
    missing = required - set(chromosomes.columns)
    if missing:
        raise ValueError(f"Chromosome file is missing columns: {sorted(missing)}")

    chromosomes["Chromosome"] = chromosomes["Chromosome"].astype(str)
    lengths = chromosomes["Total length (bp)"].astype(np.int64).to_numpy()
    starts = np.concatenate(([0], np.cumsum(lengths[:-1]))).astype(np.int64)
    chromosomes["Start"] = starts

    total_bp = int(lengths.sum())
    bounds = [int(value) for value in starts] + [total_bp]
    ticks = [int((bounds[i] + bounds[i + 1]) / 2) for i in range(len(bounds) - 1)]
    labels = [
        "\n" + str(i) if isinstance(i, int) and i > 10 and i % 2 == 0 else str(i)
        for i in list(range(1, 23)) + ["X", "Y"]
    ]
    y_labels = [
        str(i) + "    " if isinstance(i, int) and i > 10 and i % 2 == 0 else str(i)
        for i in list(range(1, 23)) + ["X", "Y"]
    ]

    return ChromosomeLayout(
        chromosomes=chromosomes,
        total_bp=total_bp,
        start_map={
            str(chromosome): int(start)
            for chromosome, start in zip(chromosomes["Chromosome"], starts)
        },
        bounds=bounds,
        ticks=ticks,
        labels=labels,
        y_labels=y_labels,
    )


def load_cytobands(cytoband_file: str | Path) -> pd.DataFrame:
    path = Path(cytoband_file)
    if not path.is_file():
        raise FileNotFoundError(f"Cytoband file does not exist: {path}")

    cytobands = pd.read_csv(path, sep="\t", header=None)
    cytobands.columns = ["chrom", "start", "end", "band", "stain"]
    cytobands = cytobands[
        cytobands["chrom"].astype(str).str.match(r"^chr([0-9]+|X|Y)$", na=False)
    ].copy()
    cytobands["chromosome"] = cytobands["chrom"].str[3:]
    cytobands["arm"] = cytobands["chromosome"] + cytobands["band"].str[0]
    return cytobands


def build_cytoband_arm_bounds(
    cytobands: pd.DataFrame,
    start_map: dict[str, int],
) -> pd.DataFrame:
    arm_bounds = (
        cytobands.groupby(["chromosome", "arm"], as_index=False)
        .agg(arm_start=("start", "min"), arm_end=("end", "max"))
        .copy()
    )
    arm_bounds["chromosome_start"] = arm_bounds["chromosome"].map(start_map)
    arm_bounds = arm_bounds.dropna(subset=["chromosome_start"])
    arm_bounds["cnv.arm.start"] = (
        arm_bounds["chromosome_start"].astype(np.int64)
        + arm_bounds["arm_start"].astype(np.int64)
    )
    arm_bounds["cnv.arm.end"] = (
        arm_bounds["chromosome_start"].astype(np.int64)
        + arm_bounds["arm_end"].astype(np.int64)
    )
    order = {str(i): i for i in range(1, 23)}
    order.update({"X": 23, "Y": 24})
    arm_bounds["chromosome_order"] = arm_bounds["chromosome"].map(order)
    arm_bounds["arm_order"] = arm_bounds["arm"].str[-1].map({"p": 0, "q": 1})
    return arm_bounds.sort_values(["chromosome_order", "arm_order"]).reset_index(
        drop=True
    )


def load_gene_annotation(
    fasta_file: str | Path,
    layout: ChromosomeLayout,
    cytobands: pd.DataFrame,
) -> pd.DataFrame:
    path = Path(fasta_file)
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file does not exist: {path}")

    records = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            parsed = _parse_gencode_header(line[1:].strip())
            if parsed is not None:
                records.append(parsed)

    annotation = pd.DataFrame(
        records,
        columns=["gene_id", "gene_symbol", "chromosome", "offset"],
    )
    annotation = annotation.drop_duplicates(subset=["gene_id"], keep="first")
    annotation = annotation[annotation["chromosome"].isin(layout.start_map)].copy()
    annotation["chromosome_start"] = annotation["chromosome"].map(layout.start_map)
    annotation["location"] = (
        annotation["chromosome_start"].astype(np.int64)
        + annotation["offset"].astype(np.int64)
    )
    annotation["cytoband"] = assign_cytobands(
        annotation["chromosome"].to_numpy(dtype=object),
        annotation["offset"].to_numpy(dtype=np.int64),
        cytobands,
    )
    annotation = annotation.dropna(subset=["cytoband"]).copy()
    annotation["arm"] = annotation["chromosome"] + annotation["cytoband"].str[0]
    return annotation


def annotate_correlation_genes(
    corr: pd.DataFrame,
    *,
    cnv_genes: np.ndarray,
    target_genes: np.ndarray,
    gene_annotation: pd.DataFrame,
    target_type: TargetType,
) -> pd.DataFrame:
    result = corr.copy()
    target_index_col = f"{target_type}.index"
    target_gene_col = f"{target_type}.gene"

    result["cnv.gene"] = _take_gene_ids(cnv_genes, result["cnv.index"])
    result[target_gene_col] = _take_gene_ids(target_genes, result[target_index_col])

    annotation = gene_annotation.set_index("gene_id")
    for prefix, gene_col in (("cnv", "cnv.gene"), (target_type, target_gene_col)):
        result[f"{prefix}.gene_symbol"] = result[gene_col].map(
            annotation["gene_symbol"]
        )
        result[f"{prefix}.location"] = result[gene_col].map(annotation["location"])
        result[f"{prefix}.chr"] = result[gene_col].map(annotation["chromosome"])
        result[f"{prefix}.cytoband"] = result[gene_col].map(annotation["cytoband"])
        result[f"{prefix}.band"] = result[gene_col].map(annotation["arm"])

    required = [
        "cnv.gene",
        target_gene_col,
        "cnv.location",
        f"{target_type}.location",
        "cnv.cytoband",
        f"{target_type}.cytoband",
        "cnv.band",
    ]
    return result.dropna(subset=required).reset_index(drop=True)


def summarize_cnv_distribution(
    annotated: pd.DataFrame,
    *,
    arm_bounds: pd.DataFrame,
    target_type: TargetType,
) -> pd.DataFrame:
    target_gene_col = f"{target_type}.gene"
    target_chr_col = f"{target_type}.chr"
    target_cytoband_col = f"{target_type}.cytoband"

    comparable = annotated["cnv.gene"] != annotated[target_gene_col]
    same_arm = (
        annotated["cnv.chr"].eq(annotated[target_chr_col])
        & annotated["cnv.cytoband"].str[0].eq(annotated[target_cytoband_col].str[0])
    )
    positive = annotated["correlation"] > 0

    pos_cis = _count_by_arm(annotated, comparable & same_arm & positive)
    pos_trans = _count_by_arm(annotated, comparable & ~same_arm & positive)
    neg_cis = _count_by_arm(annotated, comparable & same_arm & ~positive)
    neg_trans = _count_by_arm(annotated, comparable & ~same_arm & ~positive)

    counts = pd.DataFrame({"cnv.arm": arm_bounds["arm"]})
    counts["pos.cis.count"] = counts["cnv.arm"].map(pos_cis)
    counts["pos.trans.count"] = counts["cnv.arm"].map(pos_trans)
    counts["neg.cis.count"] = counts["cnv.arm"].map(neg_cis)
    counts["neg.trans.count"] = counts["cnv.arm"].map(neg_trans)
    counts = counts.merge(
        arm_bounds[["arm", "cnv.arm.start", "cnv.arm.end"]],
        left_on="cnv.arm",
        right_on="arm",
        how="left",
    ).drop(columns=["arm"])
    count_cols = [
        "pos.cis.count",
        "pos.trans.count",
        "neg.cis.count",
        "neg.trans.count",
    ]
    counts[count_cols] = counts[count_cols].fillna(0).astype(int)
    counts["count.sum"] = counts[count_cols].sum(axis=1)
    return counts


def summarize_cnv_distribution_three_class(
    annotated: pd.DataFrame,
    *,
    arm_bounds: pd.DataFrame,
    target_type: TargetType,
) -> pd.DataFrame:
    target_gene_col = f"{target_type}.gene"
    target_chr_col = f"{target_type}.chr"
    target_cytoband_col = f"{target_type}.cytoband"

    self_cis = annotated["cnv.gene"].eq(annotated[target_gene_col])
    same_arm = (
        annotated["cnv.chr"].eq(annotated[target_chr_col])
        & annotated["cnv.cytoband"].str[0].eq(annotated[target_cytoband_col].str[0])
    )
    same_arm_local = ~self_cis & same_arm
    distal_trans = ~self_cis & ~same_arm
    positive = annotated["correlation"] > 0

    counts = pd.DataFrame({"cnv.arm": arm_bounds["arm"]})
    categories = {
        "pos.self_cis.count": self_cis & positive,
        "neg.self_cis.count": self_cis & ~positive,
        "pos.same_arm_local.count": same_arm_local & positive,
        "pos.distal_trans.count": distal_trans & positive,
        "neg.same_arm_local.count": same_arm_local & ~positive,
        "neg.distal_trans.count": distal_trans & ~positive,
    }
    for column, mask in categories.items():
        counts[column] = counts["cnv.arm"].map(_count_by_arm(annotated, mask))

    counts = counts.merge(
        arm_bounds[["arm", "cnv.arm.start", "cnv.arm.end"]],
        left_on="cnv.arm",
        right_on="arm",
        how="left",
    ).drop(columns=["arm"])
    count_cols = list(categories)
    counts[count_cols] = counts[count_cols].fillna(0).astype(int)
    counts["count.sum"] = counts[count_cols].sum(axis=1)
    return counts


def summarize_gistic_segments(
    gistic_file: str | Path,
    *,
    genes: np.ndarray,
    samples: list[str],
    gene_annotation: pd.DataFrame,
) -> tuple[
    pd.DataFrame,
    list[tuple[tuple[float, float], tuple[float, float]]],
    list[tuple[tuple[float, float], tuple[float, float]]],
]:
    path = Path(gistic_file)
    if not path.is_file():
        raise FileNotFoundError(f"GISTIC file does not exist: {path}")

    gistic = pd.read_csv(path, sep="\t")
    if "idx" not in gistic.columns:
        raise ValueError(f"GISTIC file must contain an 'idx' column: {path}")

    gistic["idx"] = gistic["idx"].astype(str).str.split(".").str[0]
    gistic = gistic.set_index("idx")
    sample_columns = [sample for sample in samples if sample in gistic.columns]
    if not sample_columns:
        raise ValueError("No common sample columns found between CNV and GISTIC files")

    annotation = gene_annotation.set_index("gene_id")
    selected_genes = [
        gene for gene in genes if gene in gistic.index and gene in annotation.index
    ]
    if not selected_genes:
        raise ValueError("No common annotated genes found between CNV and GISTIC files")

    gistic = gistic.loc[selected_genes, sample_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )
    gistic = gistic.join(
        annotation[["gene_symbol", "chromosome", "cytoband", "location"]],
        how="inner",
    )
    gistic = gistic.drop_duplicates()

    values = gistic[sample_columns]
    summary = gistic[["gene_symbol", "chromosome", "cytoband", "location"]].copy()
    summary.insert(0, "gene_id", gistic.index)
    summary["amplification_count"] = (values >= 1).sum(axis=1).astype(int)
    summary["deletion_count"] = (values <= -1).sum(axis=1).astype(int)

    amp_segments = _vertical_segments(
        summary["location"].to_numpy(dtype=float),
        summary["amplification_count"].to_numpy(dtype=float),
    )
    deletion_segments = _vertical_segments(
        summary["location"].to_numpy(dtype=float),
        -summary["deletion_count"].to_numpy(dtype=float),
    )
    return summary.reset_index(drop=True), amp_segments, deletion_segments


def plot_corr_heatmap_panel(
    ax: Axes,
    annotated: pd.DataFrame,
    *,
    layout: ChromosomeLayout,
    target_type: TargetType,
    vmin: float = -1,
    vmax: float = 1,
) -> Axes:
    target_location = f"{target_type}.location"
    scatter = ax.scatter(
        _scale_bp_values(annotated["cnv.location"].to_numpy(dtype=float)),
        _scale_bp_values(annotated[target_location].to_numpy(dtype=float)),
        c=annotated["correlation"].to_numpy(dtype=float),
        cmap="bwr",
        vmin=vmin,
        vmax=vmax,
        marker="_",
        s=10,
        alpha=0.05,
        linewidths=0.7,
    )
    ax.set_xlim(0, _scale_bp_value(layout.total_bp))
    ax.set_ylim(0, _scale_bp_value(layout.total_bp))
    ax.set_xticks(_scale_bp_values(layout.bounds))
    ax.set_xticklabels([])
    ax.set_yticks(_scale_bp_values(layout.bounds))
    ax.set_yticklabels([""] + layout.y_labels)
    ax.set_ylabel("")
    ax.grid(alpha=0.2)
    ax.patch.set_visible(False)
    return scatter


def plot_cnv_distribution_panel(
    ax: Axes,
    count_df: pd.DataFrame,
    *,
    layout: ChromosomeLayout,
    top_n: int = 8,
) -> Axes:
    top_arms = set(
        count_df.sort_values("count.sum", ascending=False).head(top_n)["cnv.arm"]
    )
    colors = {
        "pos.cis.count": "pink",
        "pos.trans.count": "red",
        "neg.cis.count": "skyblue",
        "neg.trans.count": "blue",
    }

    for _, row in count_df.iterrows():
        x = _scale_bp_value(row["cnv.arm.start"])
        width = _scale_bp_value(row["cnv.arm.end"] - row["cnv.arm.start"])
        pos_cis = row["pos.cis.count"]
        pos_trans = row["pos.trans.count"]
        neg_cis = row["neg.cis.count"]
        neg_trans = row["neg.trans.count"]

        _add_rect(ax, x, 0, width, pos_cis, colors["pos.cis.count"])
        _add_rect(ax, x, pos_cis, width, pos_trans, colors["pos.trans.count"])
        _add_rect(ax, x, -neg_cis, width, neg_cis, colors["neg.cis.count"])
        _add_rect(
            ax,
            x,
            -(neg_cis + neg_trans),
            width,
            neg_trans,
            colors["neg.trans.count"],
        )

        if row["cnv.arm"] in top_arms and row["count.sum"] > 0:
            ax.text(x, pos_cis + pos_trans, row["cnv.arm"], fontsize=7)

    pos_total = count_df["pos.cis.count"] + count_df["pos.trans.count"]
    neg_total = count_df["neg.cis.count"] + count_df["neg.trans.count"]
    y_max = max(1, int(pos_total.max()))
    y_min = -max(1, int(neg_total.max()))
    ax.set_ylim(y_min * 1.15, y_max * 1.25)
    ax.set_xlim(0, _scale_bp_value(layout.total_bp))
    ax.set_xticks(_scale_bp_values(layout.bounds))
    ax.set_xticklabels([])
    ax.set_ylabel("")
    ax.tick_params(direction="out")
    ax.grid(alpha=0.2)
    return ax


def plot_self_cis_distribution_panel(
    ax: Axes,
    count_df: pd.DataFrame,
    *,
    layout: ChromosomeLayout,
    top_n: int = 8,
    show_legend: bool = True,
) -> Axes:
    top_arms = set(
        count_df.sort_values("count.sum", ascending=False).head(top_n)["cnv.arm"]
    )
    color = "red"

    for _, row in count_df.iterrows():
        x = _scale_bp_value(row["cnv.arm.start"])
        width = _scale_bp_value(row["cnv.arm.end"] - row["cnv.arm.start"])
        pos_self = int(row["pos.self_cis.count"])
        _add_rect(ax, x, 0, width, pos_self, color)
        if row["cnv.arm"] in top_arms and pos_self > 0:
            ax.text(x, pos_self, row["cnv.arm"], fontsize=7)

    y_max = max(1, int(count_df["pos.self_cis.count"].max()))
    ax.set_ylim(0, y_max * 1.35)
    ax.set_xlim(0, _scale_bp_value(layout.total_bp))
    ax.set_xticks(_scale_bp_values(layout.bounds))
    ax.set_xticklabels([])
    ax.set_ylabel("")
    ax.tick_params(direction="out")
    ax.grid(alpha=0.2)

    if show_legend:
        handle = Rectangle((0, 0), 1, 1, color=color)
        ax.legend(
            [handle],
            ["+ self-cis"],
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=7,
        )
    return ax


def plot_local_distal_distribution_panel(
    ax: Axes,
    count_df: pd.DataFrame,
    *,
    layout: ChromosomeLayout,
    top_n: int = 8,
    show_legend: bool = True,
) -> Axes:
    local_distal_cols = [
        "pos.same_arm_local.count",
        "pos.distal_trans.count",
        "neg.same_arm_local.count",
        "neg.distal_trans.count",
    ]
    top_arms = set(
        count_df.assign(_local_distal_sum=count_df[local_distal_cols].sum(axis=1))
        .sort_values("_local_distal_sum", ascending=False)
        .head(top_n)["cnv.arm"]
    )
    colors = {
        "pos.same_arm_local.count": "pink",
        "pos.distal_trans.count": "red",
        "neg.same_arm_local.count": "skyblue",
        "neg.distal_trans.count": "blue",
    }
    labels = {
        "pos.same_arm_local.count": "+ same-arm local",
        "pos.distal_trans.count": "+ distal/trans",
        "neg.same_arm_local.count": "- same-arm local",
        "neg.distal_trans.count": "- distal/trans",
    }
    positive_cols = [
        "pos.same_arm_local.count",
        "pos.distal_trans.count",
    ]
    negative_cols = [
        "neg.same_arm_local.count",
        "neg.distal_trans.count",
    ]

    for _, row in count_df.iterrows():
        x = _scale_bp_value(row["cnv.arm.start"])
        width = _scale_bp_value(row["cnv.arm.end"] - row["cnv.arm.start"])
        y_positive = 0
        for column in positive_cols:
            height = int(row[column])
            _add_rect(ax, x, y_positive, width, height, colors[column])
            y_positive += height

        y_negative = 0
        for column in negative_cols:
            height = int(row[column])
            _add_rect(ax, x, -(y_negative + height), width, height, colors[column])
            y_negative += height

        if row["cnv.arm"] in top_arms and y_positive + y_negative > 0:
            ax.text(x, y_positive, row["cnv.arm"], fontsize=7)

    pos_total = count_df[positive_cols].sum(axis=1)
    neg_total = count_df[negative_cols].sum(axis=1)
    y_max = max(1, int(pos_total.max()))
    y_min = -max(1, int(neg_total.max()))
    ax.set_ylim(y_min * 1.2, y_max * 1.35)
    ax.set_xlim(0, _scale_bp_value(layout.total_bp))
    ax.set_xticks(_scale_bp_values(layout.bounds))
    ax.set_xticklabels([])
    ax.set_ylabel("")
    ax.tick_params(direction="out")
    ax.grid(alpha=0.2)

    if show_legend:
        columns = positive_cols + negative_cols
        handles = [Rectangle((0, 0), 1, 1, color=colors[column]) for column in columns]
        ax.legend(
            handles,
            [labels[column] for column in columns],
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=7,
        )
    return ax


def plot_chromosome_panel(
    ax: Axes,
    *,
    layout: ChromosomeLayout,
    stagger_dense_labels: bool = False,
) -> Axes:
    color_list = ["black" if i % 2 == 0 else "white" for i in range(24)]
    cmap = mpl.colors.ListedColormap(color_list)
    scaled_bounds = _scale_bp_values(layout.bounds)
    norm = mpl.colors.BoundaryNorm(scaled_bounds, cmap.N)
    mpl.colorbar.ColorbarBase(
        ax,
        cmap=cmap,
        norm=norm,
        ticks=_scale_bp_values(layout.ticks),
        spacing="proportional",
        orientation="horizontal",
    )
    if stagger_dense_labels:
        ax.set_xticklabels([])
        transform = ax.get_xaxis_transform()
        for tick, label in zip(_scale_bp_values(layout.ticks), layout.labels):
            plain_label = str(label).replace("\n", "")
            y = (
                -1.35
                if plain_label.isdigit()
                and int(plain_label) >= 12
                and int(plain_label) % 2 == 0
                else -0.72
            )
            ax.text(
                tick,
                y,
                plain_label,
                transform=transform,
                ha="center",
                va="top",
                fontsize=mpl.rcParams["xtick.labelsize"],
                clip_on=False,
            )
    else:
        ax.set_xticklabels(layout.labels, rotation=0)
        plt.setp(ax.xaxis.get_majorticklabels(), ha="center")
    return ax


def plot_gistic_panel(
    ax: Axes,
    amp_segments: list[tuple[tuple[float, float], tuple[float, float]]],
    deletion_segments: list[tuple[tuple[float, float], tuple[float, float]]],
    *,
    layout: ChromosomeLayout,
) -> Axes:
    amp_x, amp_y = _segments_to_bar_arrays(amp_segments)
    deletion_x, deletion_y = _segments_to_bar_arrays(deletion_segments)
    bar_width = max(_scale_bp_value(layout.total_bp) / 6000, 0.2)
    ax.bar(amp_x, amp_y, width=bar_width, color="red", linewidth=0)
    ax.bar(deletion_x, deletion_y, width=bar_width, color="blue", linewidth=0)
    ax.set_xlim(0, _scale_bp_value(layout.total_bp))
    ax.set_ylim(-100, 100)
    ax.set_xticks(_scale_bp_values(layout.bounds))
    ax.set_xticklabels([])
    ax.set_ylabel("")
    ax.grid(alpha=0.2)
    ax.autoscale_view()
    return ax


def save_corr_heatmap_panel(
    path: str | Path,
    annotated: pd.DataFrame,
    *,
    layout: ChromosomeLayout,
    target_type: TargetType,
    dpi: int = 300,
) -> Path:
    output_path = Path(path)
    fig, ax = plt.subplots(figsize=(6, 5), dpi=dpi)
    scatter = plot_corr_heatmap_panel(
        ax,
        annotated,
        layout=layout,
        target_type=target_type,
    )
    fig.colorbar(scatter, ax=ax, fraction=0.03, pad=0.02, label="Spearman r")
    _save_png_and_editable_pdf(fig, output_path)
    plt.close(fig)
    return output_path


def save_cnv_distribution_panel(
    path: str | Path,
    count_df: pd.DataFrame,
    *,
    layout: ChromosomeLayout,
    dpi: int = 300,
) -> Path:
    output_path = Path(path)
    fig, ax = plt.subplots(figsize=(8, 2.5), dpi=dpi)
    plot_cnv_distribution_panel(ax, count_df, layout=layout)
    _save_png_and_editable_pdf(fig, output_path)
    plt.close(fig)
    return output_path


def save_chromosome_panel(
    path: str | Path,
    *,
    layout: ChromosomeLayout,
    dpi: int = 300,
) -> Path:
    output_path = Path(path)
    fig = plt.figure(figsize=(8, 0.8), dpi=dpi)
    ax = fig.add_axes([0.05, 0.35, 0.9, 0.25])
    plot_chromosome_panel(ax, layout=layout)
    _save_png_and_editable_pdf(fig, output_path)
    plt.close(fig)
    return output_path


def save_gistic_panel(
    path: str | Path,
    amp_segments: list[tuple[tuple[float, float], tuple[float, float]]],
    deletion_segments: list[tuple[tuple[float, float], tuple[float, float]]],
    *,
    layout: ChromosomeLayout,
    dpi: int = 300,
) -> Path:
    output_path = Path(path)
    fig, ax = plt.subplots(figsize=(8, 2.5), dpi=dpi)
    plot_gistic_panel(
        ax,
        amp_segments,
        deletion_segments,
        layout=layout,
    )
    _save_png_and_editable_pdf(fig, output_path)
    plt.close(fig)
    return output_path


def save_combined_figure(
    path: str | Path,
    annotated: pd.DataFrame,
    count_df: pd.DataFrame,
    *,
    layout: ChromosomeLayout,
    target_type: TargetType,
    gistic_segments: tuple[
        list[tuple[tuple[float, float], tuple[float, float]]],
        list[tuple[tuple[float, float], tuple[float, float]]],
    ]
    | None = None,
    dpi: int = 300,
) -> Path:
    output_path = Path(path)
    with plt.style.context("seaborn-v0_8-white"):
        fig = plt.figure(figsize=(5, 10), dpi=dpi)
        ax_chromosome = fig.add_axes([0.05, 0.1, 0.9, 0.02])
        ax_distribution = fig.add_axes([0.05, 0.15, 0.9, 0.1])
        if gistic_segments is None:
            ax_heatmap = fig.add_axes([0.05, 0.3, 0.9, 0.4])
        else:
            ax_gistic = fig.add_axes([0.05, 0.28, 0.9, 0.1])
            ax_heatmap = fig.add_axes([0.05, 0.4, 0.9, 0.4])

        plot_chromosome_panel(ax_chromosome, layout=layout)
        plot_cnv_distribution_panel(ax_distribution, count_df, layout=layout)
        ax_distribution.set_ylabel("")
        if gistic_segments is not None:
            plot_gistic_panel(
                ax_gistic,
                gistic_segments[0],
                gistic_segments[1],
                layout=layout,
            )
            ax_gistic.set_ylabel("")
        plot_corr_heatmap_panel(
            ax_heatmap,
            annotated,
            layout=layout,
            target_type=target_type,
        )
        ax_heatmap.set_ylabel("")
        _save_png_and_editable_pdf(fig, output_path)
        plt.close(fig)
    return output_path


def save_combined_three_class_split_self_figure(
    path: str | Path,
    annotated: pd.DataFrame,
    count_df: pd.DataFrame,
    *,
    layout: ChromosomeLayout,
    target_type: TargetType,
    gistic_segments: tuple[
        list[tuple[tuple[float, float], tuple[float, float]]],
        list[tuple[tuple[float, float], tuple[float, float]]],
    ]
    | None = None,
    dpi: int = 300,
    heatmap_vmin: float = -0.5,
    heatmap_vmax: float = 0.5,
) -> Path:
    output_path = Path(path)
    with plt.style.context("seaborn-v0_8-white"):
        fig = plt.figure(figsize=(6.8, 12.5), dpi=dpi)
        left = 0.08
        width = 0.68

        ax_chromosome = fig.add_axes([left, 0.06, width, 0.018])
        ax_local_distal = fig.add_axes([left, 0.12, width, 0.115])
        ax_self = fig.add_axes([left, 0.28, width, 0.095])
        if gistic_segments is None:
            ax_heatmap = fig.add_axes([left, 0.45, width, 0.39])
            ax_heatmap_colorbar = fig.add_axes([0.79, 0.45, 0.025, 0.39])
        else:
            ax_gistic = fig.add_axes([left, 0.41, width, 0.095])
            ax_heatmap = fig.add_axes([left, 0.56, width, 0.36])
            ax_heatmap_colorbar = fig.add_axes([0.79, 0.56, 0.025, 0.36])

        plot_chromosome_panel(
            ax_chromosome,
            layout=layout,
            stagger_dense_labels=True,
        )
        plot_local_distal_distribution_panel(
            ax_local_distal,
            count_df,
            layout=layout,
        )
        ax_local_distal.set_ylabel("")
        plot_self_cis_distribution_panel(ax_self, count_df, layout=layout)
        ax_self.set_ylabel("")

        if gistic_segments is not None:
            plot_gistic_panel(
                ax_gistic,
                gistic_segments[0],
                gistic_segments[1],
                layout=layout,
            )
            ax_gistic.set_ylabel("")
            gistic_handles = [
                Rectangle((0, 0), 1, 1, color="red"),
                Rectangle((0, 0), 1, 1, color="blue"),
            ]
            ax_gistic.legend(
                gistic_handles,
                ["Amplification", "Deletion"],
                loc="center left",
                bbox_to_anchor=(1.01, 0.5),
                frameon=False,
                fontsize=7,
            )

        plot_corr_heatmap_panel(
            ax_heatmap,
            annotated,
            layout=layout,
            target_type=target_type,
            vmin=heatmap_vmin,
            vmax=heatmap_vmax,
        )
        ax_heatmap.set_ylabel("")
        colorbar_mappable = mpl.cm.ScalarMappable(
            norm=mpl.colors.Normalize(vmin=heatmap_vmin, vmax=heatmap_vmax),
            cmap="bwr",
        )
        colorbar_mappable.set_array([])
        colorbar = fig.colorbar(colorbar_mappable, cax=ax_heatmap_colorbar)
        colorbar.set_label("Spearman r", fontsize=7)
        colorbar.ax.tick_params(labelsize=7)

        _save_png_and_editable_pdf(fig, output_path)
        plt.close(fig)
    return output_path


def assign_cytobands(
    chromosomes: np.ndarray,
    positions: np.ndarray,
    cytobands: pd.DataFrame,
) -> np.ndarray:
    lookup = {}
    for chromosome, table in cytobands.groupby("chromosome"):
        ordered = table.sort_values("start")
        lookup[str(chromosome)] = (
            ordered["start"].to_numpy(dtype=np.int64),
            ordered["end"].to_numpy(dtype=np.int64),
            ordered["band"].to_numpy(dtype=object),
        )

    result = np.empty(len(chromosomes), dtype=object)
    result[:] = None
    for i, (chromosome, position) in enumerate(zip(chromosomes, positions)):
        entry = lookup.get(str(chromosome))
        if entry is None:
            continue
        starts, ends, bands = entry
        idx = np.searchsorted(ends, int(position), side="left")
        if idx < len(starts) and starts[idx] <= int(position) <= ends[idx]:
            result[i] = bands[idx]
    return result


def _parse_gencode_header(description: str) -> tuple[str, str, str, int] | None:
    items = description.split("|")
    try:
        gene_id = next(item for item in items if item.startswith("ENSG")).split(".")[0]
        gene_symbol = next(item for item in items if item.startswith("GN=")).split(
            "=",
            1,
        )[1]
        location = next(item for item in items if item.startswith("chr"))
        chromosome, offset = location.split("-", 1)[0].split(":", 1)
        chromosome = chromosome[3:]
        return gene_id, gene_symbol, chromosome, int(offset)
    except (StopIteration, ValueError, IndexError):
        return None


def _take_gene_ids(genes: np.ndarray, indices: pd.Series) -> np.ndarray:
    numeric = pd.to_numeric(indices, errors="coerce").to_numpy(dtype=float)
    result = np.empty(len(numeric), dtype=object)
    result[:] = None
    valid = np.isfinite(numeric)
    int_indices = numeric[valid].astype(np.int64)
    valid_positions = np.flatnonzero(valid)
    in_range = (int_indices >= 0) & (int_indices < len(genes))
    result[valid_positions[in_range]] = genes[int_indices[in_range]]
    return result


def _count_by_arm(annotated: pd.DataFrame, mask: pd.Series) -> pd.Series:
    return annotated.loc[mask].groupby("cnv.band").size()


def _add_rect(
    ax: Axes,
    startx: float,
    starty: float,
    width: float,
    height: float,
    color: str,
) -> None:
    if height > 0:
        ax.add_patch(Rectangle((startx, starty), width, height, color=color))


def _save_png_and_editable_pdf(fig: plt.Figure, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")

    pdf_path = output_path.with_suffix(".pdf")
    with mpl.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42}):
        fig.savefig(pdf_path, bbox_inches="tight")

    tiff_path = output_path.with_suffix(".tiff")
    fig.savefig(tiff_path, bbox_inches="tight", dpi=600)


def _scale_bp_value(value: float | int) -> float:
    return float(value) / BP_PLOT_SCALE


def _scale_bp_values(values: list[int] | np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=float) / BP_PLOT_SCALE


def _segments_to_bar_arrays(
    segments: list[tuple[tuple[float, float], tuple[float, float]]],
) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray([segment[0][0] for segment in segments], dtype=float) / BP_PLOT_SCALE
    y = np.asarray([segment[1][1] for segment in segments], dtype=float)
    return x, y


def _vertical_segments(
    positions: np.ndarray,
    values: np.ndarray,
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    return [
        ((float(position), 0.0), (float(position), float(value)))
        for position, value in zip(positions, values)
    ]
