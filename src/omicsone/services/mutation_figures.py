from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "DejaVu Sans"]
matplotlib.rcParams["font.weight"] = "regular"
matplotlib.rcParams["axes.titleweight"] = "regular"
matplotlib.rcParams["axes.labelweight"] = "regular"
matplotlib.rcParams["font.size"] = 16
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


DEFAULT_MAF_PATH = (
    r"C:\Users\yhu39\OneDrive - Johns Hopkins\project\2025_HN-Lung\google"
    r"\Head and Neck & Lung\HNSCC\HNSCC_somatic_mutation.maf"
)
DEFAULT_MUTATION_EXCEL_PATH = (
    r"C:\Users\yhu39\OneDrive - Johns Hopkins\project\2025_HN-Lung\google"
    r"\Head and Neck & Lung\Processing\HNSCC_somatic_mutation_gene_level_revised.xlsx"
)
DEFAULT_OUTPUT_DIR = r"E:\lab\HSinI\runs\20260501_hnscc_mutations"
DEFAULT_COHORT = "hnscc"
DEFAULT_MUTATION_BINARY_PATH = (
    r"F:\lab\HsinI\Head and Neck & Lung\HNSCC"
    r"\HNSCC_somatic_mutation_gene_level_binary.txt"
)
DEFAULT_MUTATION_META_PATH = (
    r"F:\lab\HsinI\Head and Neck & Lung\HNSCC\HNSCC_meta.txt"
)
DEFAULT_V2_MAF_PATH = (
    r"F:\lab\HsinI\Head and Neck & Lung\HNSCC\HNSCC_somatic_mutation.maf"
)
DEFAULT_V2_OUTPUT_DIR = r"E:\lab\HSinI\runs\20260501_hnscc_mutations_v2"
DEFAULT_V2_COHORT = "hnscc_v2"
RESOURCE_DIR = Path(__file__).resolve().parents[1] / "resources"
DEFAULT_HUMAN_GENE_SYMBOL_MAP_PATH = (
    RESOURCE_DIR / "gene_symbol_maps" / "gencode_v42_human_gene_symbol_map.tsv"
)
DEFAULT_MOUSE_GENE_SYMBOL_MAP_PATH = (
    RESOURCE_DIR / "gene_symbol_maps" / "gencode_m31_mouse_gene_symbol_map.tsv"
)
GENE_COLUMN_CANDIDATES = (
    "GENE NAME",
    "gene name",
    "Gene Name",
    "gene_name",
    "Gene",
    "gene",
    "Hugo_Symbol",
    "HUGO_SYMBOL",
    "Symbol",
)

GRADE_COLOR_MAP = {
    "G1 Well differentiated": "#4daf4a",
    "G2 Moderately differentiated": "#377eb8",
    "G3 Poorly differentiated": "#ff7f00",
}

STAGE_COLOR_MAP = {
    "Stage I": "#4daf4a",
    "Stage II": "#377eb8",
    "Stage III": "#ff7f00",
    "Stage IV": "#984ea3",
}

CUSTOM_MUTATION_COLORS = [
    "white",
    "#d3d3d3",
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#ff9896",
    "#98df8a",
    "#c49c94",
    "#f7b6d2",
    "#9edae5",
    "#2b83ba",
    "#ffcc00",
    "#ff5511",
]

TARGET_HEATMAP_WIDTH_RATIO = 2.0
TARGET_HEATMAP_HEIGHT_RATIO = 1.25
TARGET_HEATMAP_ASPECT = TARGET_HEATMAP_WIDTH_RATIO / TARGET_HEATMAP_HEIGHT_RATIO
BASE_HEATMAP_SAMPLE_COUNT = 107
BASE_HEATMAP_GENE_COUNT = 14
BASE_HEATMAP_WIDTH_INCH = 10.0
BASE_HEATMAP_HEIGHT_INCH = BASE_HEATMAP_WIDTH_INCH / TARGET_HEATMAP_ASPECT
HEATMAP_WIDTH_PER_SAMPLE = BASE_HEATMAP_WIDTH_INCH / BASE_HEATMAP_SAMPLE_COUNT
HEATMAP_HEIGHT_PER_GENE = BASE_HEATMAP_HEIGHT_INCH / BASE_HEATMAP_GENE_COUNT
MIN_HEATMAP_WIDTH_INCH = 6.0
MIN_HEATMAP_HEIGHT_INCH = 4.0
BASE_MUTATION_TYPE_WIDTH_INCH = 12.0
BASE_MUTATION_TYPE_HEIGHT_INCH = 12.0
MUTATION_TYPE_HEIGHT_PER_GENE = BASE_MUTATION_TYPE_HEIGHT_INCH / BASE_HEATMAP_GENE_COUNT
MIN_MUTATION_TYPE_HEIGHT_INCH = 6.0


@dataclass
class MutationFigureResult:
    heatmap_pdf: Path
    mutation_type_pdf: Path
    result_log: Path
    gene_summary_tsv: Path
    binary_matrix_tsv: Path
    mutation_type_matrix_tsv: Path
    encoded_matrix_tsv: Path
    sample_annotations_tsv: Path
    sample_annotation_colors_tsv: Path
    mutation_color_table_tsv: Path
    filtered_maf_tsv: Path
    gene_count: int
    sample_count: int
    cohort: str
    total_maf_rows: int
    filtered_maf_rows: int
    found_mutations: list[str]
    heatmap_width_inch: float
    heatmap_height_inch: float
    heatmap_aspect: float
    mutation_type_width_inch: float
    mutation_type_height_inch: float
    mutation_type_aspect: float
    sample_gene_ratio: float


@dataclass
class MutationFigureV2Result(MutationFigureResult):
    mutation_binary_path: Path
    meta_path: Path
    gene_symbol_map_path: Path
    species: str


def _ensure_file(path: str | Path, label: str) -> Path:
    resolved = Path(path).expanduser()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} does not exist: {resolved}")
    return resolved


def _normalize_cohort(cohort: str | None) -> str:
    value = (cohort or DEFAULT_COHORT).strip().lower()
    if not value:
        value = DEFAULT_COHORT
    return value


def _normalize_species(species: str | None) -> str:
    value = (species or "human").strip().lower()
    if value not in {"human", "mouse"}:
        raise ValueError("species must be either 'human' or 'mouse'")
    return value


def _default_gene_symbol_map_path(species: str) -> Path:
    if species == "human":
        return DEFAULT_HUMAN_GENE_SYMBOL_MAP_PATH
    return DEFAULT_MOUSE_GENE_SYMBOL_MAP_PATH


def _resolve_gene_column(mutation_df: pd.DataFrame) -> str:
    normalized_columns = {str(column).strip().lower(): column for column in mutation_df.columns}
    for candidate in GENE_COLUMN_CANDIDATES:
        match = normalized_columns.get(candidate.lower())
        if match is not None:
            return match

    raise ValueError(
        "Mutation Excel file must include a gene-name column. "
        f"Tried: {', '.join(GENE_COLUMN_CANDIDATES)}"
    )


def _metadata_columns(mutation_df: pd.DataFrame, gene_column: str) -> list[str]:
    candidates = ["NUM_MUT", gene_column, "idx"]
    return [column for column in candidates if column in mutation_df.columns]


def _normalize_gene_name(value: object) -> object:
    if pd.isna(value):
        return value

    gene_name = str(value).strip().split(",", maxsplit=1)[0].strip()
    if gene_name == "MUC16(CA125?)":
        return "MUC16"
    return gene_name


def _normalize_ensembl_gene_id(value: object) -> str:
    return str(value).strip().split(".", maxsplit=1)[0]


def _load_gene_symbol_map(path: str | Path) -> dict[str, str]:
    mapping_path = _ensure_file(path, "Gene symbol map file")
    mapping_df = pd.read_csv(mapping_path, sep="\t")
    required_columns = {"gene_id", "gene_symbol"}
    if not required_columns.issubset(mapping_df.columns):
        raise ValueError(
            "Gene symbol map must be a TSV with columns: gene_id, gene_symbol"
        )

    mapping_df = mapping_df.dropna(subset=["gene_id", "gene_symbol"])
    return {
        _normalize_ensembl_gene_id(row["gene_id"]): str(row["gene_symbol"]).strip()
        for _, row in mapping_df.iterrows()
    }


def _build_column_annotations(
    mutation_df: pd.DataFrame,
    gene_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta_df = mutation_df.iloc[:2, :]
    meta_df = meta_df.drop(_metadata_columns(mutation_df, gene_column), axis=1)
    meta_df = meta_df.dropna(axis=1, how="all")
    meta_df = meta_df.loc[:, ~meta_df.columns.astype(str).str.startswith("Unnamed:")]

    rows = []
    for sample in meta_df.columns:
        grade, stage = meta_df[sample].values
        rows.append([sample, grade, stage])

    col_df = pd.DataFrame(rows, columns=["Sample", "Grade", "Stage"]).set_index("Sample")
    col_color_df = col_df.copy(deep=True)
    col_color_df["Grade"] = col_df["Grade"].map(GRADE_COLOR_MAP)
    col_color_df["Stage"] = col_df["Stage"].map(STAGE_COLOR_MAP)
    return col_df, col_color_df


def _build_mutation_color_table(maf_df: pd.DataFrame) -> pd.DataFrame:
    variant_classes = sorted(maf_df["Variant_Classification"].dropna().unique())
    if len(variant_classes) + 2 > len(CUSTOM_MUTATION_COLORS):
        raise ValueError(
            "Not enough configured colors for mutation classes: "
            f"{len(variant_classes)} classes found."
        )

    rows = [[0, "No Mutation", CUSTOM_MUTATION_COLORS[0]]]
    for index, mutation in enumerate(variant_classes):
        rows.append([index + 1, mutation, CUSTOM_MUTATION_COLORS[index + 1]])
    rows.append([len(variant_classes) + 1, "Multiple_Mutations", CUSTOM_MUTATION_COLORS[len(variant_classes) + 1]])
    return pd.DataFrame(rows, columns=["idx", "mut", "color"])


def _build_maf_map(maf_df: pd.DataFrame, target_genes: set[str]) -> tuple[dict[str, dict[str, str]], pd.DataFrame]:
    filtered_maf = maf_df[maf_df["Hugo_Symbol"].isin(target_genes)]
    maf_map: dict[str, dict[str, list[str] | str]] = {}

    for _, row in filtered_maf.iterrows():
        gene = row["Hugo_Symbol"]
        sample = row["Tumor_Sample_Barcode"]
        mutation_type = row["Variant_Classification"]

        maf_map.setdefault(sample, {})
        maf_map[sample].setdefault(gene, [])
        maf_map[sample][gene].append(mutation_type)

    normalized: dict[str, dict[str, str]] = {}
    for sample, gene_map in maf_map.items():
        normalized[sample] = {}
        for gene, mutation_types in gene_map.items():
            values = list(mutation_types)
            normalized[sample][gene] = values[0] if len(values) == 1 else "Multiple_Mutations"

    return normalized, filtered_maf


def _prepare_mutation_tables(
    mutation_df: pd.DataFrame,
    maf_df: pd.DataFrame,
    mutation_threshold: float,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    list[str],
]:
    gene_column = _resolve_gene_column(mutation_df)
    col_df, col_color_df = _build_column_annotations(mutation_df, gene_column)

    filtered = mutation_df.loc[2:, :]
    filtered = filtered[filtered["NUM_MUT"].astype(float) >= mutation_threshold]
    filtered = filtered.copy()
    filtered[gene_column] = filtered[gene_column].map(_normalize_gene_name)

    target_genes = set(filtered[gene_column].dropna())
    maf_map, filtered_maf = _build_maf_map(maf_df, target_genes)
    mutation_color_df = _build_mutation_color_table(maf_df)
    mutation_color_map = {
        row["mut"]: (int(row["idx"]), row["color"])
        for _, row in mutation_color_df.iterrows()
    }

    metadata_cols = set(_metadata_columns(filtered, gene_column))
    sample_cols = [
        column
        for column in filtered.columns
        if column not in metadata_cols and column in col_df.index
    ]
    binary_df = filtered.set_index(gene_column)[sample_cols]

    pairs = []
    for sample in binary_df.columns:
        signature = [int(value) for value in list(binary_df[sample])]
        pairs.append((sample, signature))
    pairs = sorted(pairs, key=lambda item: item[1], reverse=True)
    sorted_binary_df = binary_df[[sample for sample, _ in pairs]]

    encoded_df = sorted_binary_df.copy(deep=True)
    typed_df = sorted_binary_df.copy(deep=True).astype(object)
    found_mutations: set[str] = set()

    for gene, row in sorted_binary_df.iterrows():
        for sample in sorted_binary_df.columns:
            if int(row[sample]) == 1:
                mutation_type = maf_map.get(sample, {}).get(gene, "No Mutation")
                found_mutations.add(mutation_type)
                encoded_df.loc[gene, sample] = mutation_color_map[mutation_type][0]
                typed_df.loc[gene, sample] = mutation_type

    encoded_df = encoded_df.astype(float)
    source_num_mut = pd.to_numeric(
        filtered.set_index(gene_column)["NUM_MUT"],
        errors="coerce",
    )
    maf_row_counts = filtered_maf["Hugo_Symbol"].value_counts()
    gene_summary_df = pd.DataFrame(
        {
            "gene": sorted_binary_df.index,
            "source_num_mut": [source_num_mut.get(gene, np.nan) for gene in sorted_binary_df.index],
            "mutated_sample_count": sorted_binary_df.astype(int).sum(axis=1).to_numpy(),
            "mutation_ratio": (sorted_binary_df.astype(int).sum(axis=1) / sorted_binary_df.shape[1]).to_numpy(),
            "maf_row_count": [int(maf_row_counts.get(gene, 0)) for gene in sorted_binary_df.index],
        }
    )
    sorted_col_colors = col_color_df.loc[encoded_df.columns]
    sorted_col_df = col_df.loc[encoded_df.columns]
    return (
        encoded_df,
        typed_df,
        sorted_binary_df,
        sorted_col_df,
        sorted_col_colors,
        mutation_color_df,
        filtered_maf,
        gene_summary_df,
        sorted(found_mutations),
    )


def _build_column_annotations_from_meta(
    meta_df: pd.DataFrame,
    sample_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for sample in sample_columns:
        if sample not in meta_df.index:
            continue
        row = meta_df.loc[sample]
        try:
            stage = row[("Stage", "ORD")]
            grade = row[("Histologic_Grade", "ORD")]
        except KeyError as exc:
            raise ValueError(
                "Meta file must contain ('Stage', 'ORD') and "
                "('Histologic_Grade', 'ORD') columns."
            ) from exc
        rows.append([sample, grade, stage])

    col_df = pd.DataFrame(rows, columns=["Sample", "Grade", "Stage"]).set_index("Sample")
    col_color_df = col_df.copy(deep=True)
    col_color_df["Grade"] = col_df["Grade"].map(GRADE_COLOR_MAP)
    col_color_df["Stage"] = col_df["Stage"].map(STAGE_COLOR_MAP)
    return col_df, col_color_df


def _prepare_mutation_tables_from_binary(
    mutation_binary_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    maf_df: pd.DataFrame,
    gene_symbol_map: dict[str, str],
    mutation_threshold: float,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    list[str],
]:
    normalized_binary_df = mutation_binary_df.copy()
    normalized_binary_df.index = [
        _normalize_ensembl_gene_id(index) for index in normalized_binary_df.index
    ]
    normalized_binary_df = normalized_binary_df.apply(pd.to_numeric, errors="coerce").fillna(0)

    sample_cols = [column for column in normalized_binary_df.columns if column in meta_df.index]
    if not sample_cols:
        raise ValueError("No overlapping samples found between mutation binary file and meta file.")
    normalized_binary_df = normalized_binary_df[sample_cols]

    mutation_ratio = normalized_binary_df.sum(axis=1) / len(sample_cols)
    filtered_gene_ids = mutation_ratio[mutation_ratio >= mutation_threshold].sort_values(ascending=False).index

    rows = []
    missing_gene_symbols = []
    for gene_id in filtered_gene_ids:
        gene_symbol = gene_symbol_map.get(gene_id)
        if not gene_symbol:
            missing_gene_symbols.append(gene_id)
            continue
        rows.append((gene_id, gene_symbol))

    if not rows:
        raise ValueError(
            "No genes passed the mutation threshold after gene-symbol mapping. "
            f"Missing mapped symbols for {len(missing_gene_symbols)} genes."
        )

    selected_gene_ids = [gene_id for gene_id, _ in rows]
    selected_gene_symbols = [gene_symbol for _, gene_symbol in rows]
    gene_id_by_symbol = {
        gene_symbol: gene_id for gene_id, gene_symbol in zip(selected_gene_ids, selected_gene_symbols)
    }
    binary_df = normalized_binary_df.loc[selected_gene_ids].copy()
    binary_df.index = selected_gene_symbols

    col_df, col_color_df = _build_column_annotations_from_meta(meta_df, sample_cols)
    binary_df = binary_df[col_df.index]

    target_genes = set(binary_df.index)
    maf_map, filtered_maf = _build_maf_map(maf_df, target_genes)
    mutation_color_df = _build_mutation_color_table(maf_df)
    mutation_color_map = {
        row["mut"]: (int(row["idx"]), row["color"])
        for _, row in mutation_color_df.iterrows()
    }

    pairs = []
    for sample in binary_df.columns:
        signature = [int(value) for value in list(binary_df[sample])]
        pairs.append((sample, signature))
    pairs = sorted(pairs, key=lambda item: item[1], reverse=True)
    sorted_binary_df = binary_df[[sample for sample, _ in pairs]]

    encoded_df = sorted_binary_df.copy(deep=True)
    typed_df = sorted_binary_df.copy(deep=True).astype(object)
    found_mutations: set[str] = set()

    for gene, row in sorted_binary_df.iterrows():
        for sample in sorted_binary_df.columns:
            if int(row[sample]) == 1:
                mutation_type = maf_map.get(sample, {}).get(gene, "No Mutation")
                found_mutations.add(mutation_type)
                encoded_df.loc[gene, sample] = mutation_color_map[mutation_type][0]
                typed_df.loc[gene, sample] = mutation_type

    encoded_df = encoded_df.astype(float)
    maf_row_counts = filtered_maf["Hugo_Symbol"].value_counts()
    gene_summary_df = pd.DataFrame(
        {
            "gene": sorted_binary_df.index,
            "gene_id": [gene_id_by_symbol.get(gene, "") for gene in sorted_binary_df.index],
            "mutated_sample_count": sorted_binary_df.astype(int).sum(axis=1).to_numpy(),
            "mutation_ratio": (sorted_binary_df.astype(int).sum(axis=1) / sorted_binary_df.shape[1]).to_numpy(),
            "maf_row_count": [int(maf_row_counts.get(gene, 0)) for gene in sorted_binary_df.index],
        }
    )
    sorted_col_colors = col_color_df.loc[encoded_df.columns]
    sorted_col_df = col_df.loc[encoded_df.columns]
    return (
        encoded_df,
        typed_df,
        sorted_binary_df,
        sorted_col_df,
        sorted_col_colors,
        mutation_color_df,
        filtered_maf,
        gene_summary_df,
        sorted(found_mutations),
    )


def _calculate_heatmap_figure_size(sample_count: int, gene_count: int) -> tuple[float, float]:
    width = max(MIN_HEATMAP_WIDTH_INCH, sample_count * HEATMAP_WIDTH_PER_SAMPLE)
    height = max(MIN_HEATMAP_HEIGHT_INCH, gene_count * HEATMAP_HEIGHT_PER_GENE)
    return width, height


def _calculate_mutation_type_figure_size(gene_count: int) -> tuple[float, float]:
    height = max(MIN_MUTATION_TYPE_HEIGHT_INCH, gene_count * MUTATION_TYPE_HEIGHT_PER_GENE)
    return BASE_MUTATION_TYPE_WIDTH_INCH, height


def _plot_mutation_heatmap(
    encoded_df: pd.DataFrame,
    col_color_df: pd.DataFrame,
    mutation_color_df: pd.DataFrame,
    found_mutations: list[str],
    output_path: Path,
    figure_size: tuple[float, float],
) -> None:
    cmap = ListedColormap(CUSTOM_MUTATION_COLORS)
    grid = sns.clustermap(
        encoded_df,
        cmap=cmap,
        annot=False,
        col_cluster=False,
        row_cluster=False,
        col_colors=col_color_df,
        figsize=figure_size,
        linewidths=0.5,
        linecolor="grey",
        cbar_pos=None,
    )

    row_nonzero_percentage = (encoded_df != 0).sum(axis=1) / encoded_df.shape[1] * 100
    heatmap_ax = grid.ax_heatmap
    for spine in heatmap_ax.spines.values():
        spine.set_visible(True)

    fig = grid.fig
    bar_ax = fig.add_axes(
        [
            heatmap_ax.get_position().x1 + 0.01,
            heatmap_ax.get_position().y0 - 0.025,
            0.1,
            heatmap_ax.get_position().height * 1.09,
        ]
    )

    row_order = (
        grid.dendrogram_row.reordered_ind
        if grid.dendrogram_row is not None
        else np.arange(encoded_df.shape[0])[::-1]
    )
    bars = bar_ax.barh(
        np.arange(len(row_order)),
        row_nonzero_percentage.iloc[row_order],
        color="skyblue",
        edgecolor="black",
    )
    bar_ax.set_yticks([])
    bar_ax.set_xticks([0, 50, 100])
    bar_ax.set_xlabel("% samples \n with mutations")
    bar_ax.set_facecolor("none")
    bar_ax.patch.set_alpha(0)
    for spine in bar_ax.spines.values():
        spine.set_visible(False)
    for bar in bars:
        width = bar.get_width()
        y = bar.get_y() + bar.get_height() / 2
        bar_ax.text(width, y, f"{width:.0f}", va="center", ha="left", color="black", fontsize=16)

    col_colors_ax = grid.ax_col_colors
    col_colors_ax.set_position(
        [
            col_colors_ax.get_position().x0,
            col_colors_ax.get_position().y0,
            col_colors_ax.get_position().width,
            col_colors_ax.get_position().height * 2,
        ]
    )

    mutation_color_map = dict(zip(mutation_color_df["mut"], mutation_color_df["color"]))
    mutation_legend = [
        Patch(facecolor=mutation_color_map[mutation], label=mutation)
        for mutation in mutation_color_map
        if mutation in found_mutations
    ]
    stage_legend = [
        Patch(facecolor=color, label=stage.replace("Stage ", ""))
        for stage, color in STAGE_COLOR_MAP.items()
    ]
    grade_legend = [
        Patch(facecolor="#4daf4a", label="I"),
        Patch(facecolor="#377eb8", label="II"),
        Patch(facecolor="#ff7f00", label="III"),
    ]

    legend_stage = plt.legend(
        handles=stage_legend,
        title="Stage",
        loc="upper left",
        bbox_to_anchor=(1.15, 0.87),
        bbox_transform=fig.transFigure,
    )
    legend_grade = plt.legend(
        handles=grade_legend,
        title="Grade",
        loc="upper left",
        bbox_to_anchor=(1.05, 0.87),
        bbox_transform=fig.transFigure,
    )
    legend_mutation = plt.legend(
        handles=mutation_legend,
        title="Mutation",
        loc="upper left",
        bbox_to_anchor=(1.05, 0.6),
        bbox_transform=fig.transFigure,
    )
    plt.gca().add_artist(legend_grade)
    plt.gca().add_artist(legend_stage)
    plt.gca().add_artist(legend_mutation)

    grid.ax_heatmap.set_xticks([])
    grid.ax_heatmap.set_yticks(np.arange(encoded_df.shape[0]) + 0.5)
    grid.ax_heatmap.set_yticklabels(encoded_df.index.tolist(), rotation=0)
    grid.ax_heatmap.yaxis.tick_left()
    grid.ax_heatmap.yaxis.set_label_position("left")

    fig.savefig(output_path, format="pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _write_result_log(
    output_path: Path,
    cohort: str,
    mutation_excel_path: Path,
    maf_path: Path,
    heatmap_pdf: Path,
    mutation_type_pdf: Path,
    gene_summary_tsv: Path,
    binary_matrix_tsv: Path,
    mutation_type_matrix_tsv: Path,
    encoded_matrix_tsv: Path,
    sample_annotations_tsv: Path,
    sample_annotation_colors_tsv: Path,
    mutation_color_table_tsv: Path,
    filtered_maf_tsv: Path,
    sample_count: int,
    gene_count: int,
    heatmap_width: float,
    heatmap_height: float,
    mutation_type_width: float,
    mutation_type_height: float,
    total_maf_rows: int,
    filtered_maf_rows: int,
    found_mutations: list[str],
) -> None:
    sample_gene_ratio = sample_count / gene_count if gene_count else float("inf")
    figure_aspect = heatmap_width / heatmap_height if heatmap_height else float("inf")
    mutation_type_aspect = mutation_type_width / mutation_type_height if mutation_type_height else float("inf")
    lines = [
        "OmicsOne mutation figure result",
        "",
        f"mutation_excel_path = {mutation_excel_path}",
        f"maf_path = {maf_path}",
        f"heatmap_pdf = {heatmap_pdf}",
        f"mutation_type_pdf = {mutation_type_pdf}",
        f"gene_summary_tsv = {gene_summary_tsv}",
        f"binary_matrix_tsv = {binary_matrix_tsv}",
        f"mutation_type_matrix_tsv = {mutation_type_matrix_tsv}",
        f"encoded_matrix_tsv = {encoded_matrix_tsv}",
        f"sample_annotations_tsv = {sample_annotations_tsv}",
        f"sample_annotation_colors_tsv = {sample_annotation_colors_tsv}",
        f"mutation_color_table_tsv = {mutation_color_table_tsv}",
        f"filtered_maf_tsv = {filtered_maf_tsv}",
        "",
        f"cohort = {cohort}",
        f"sample_count = {sample_count}",
        f"gene_count = {gene_count}",
        f"sample_gene_ratio = {sample_count}:{gene_count} ({sample_gene_ratio:.4f}:1)",
        "",
        f"target_heatmap_width_height_ratio = {TARGET_HEATMAP_WIDTH_RATIO}:{TARGET_HEATMAP_HEIGHT_RATIO} ({TARGET_HEATMAP_ASPECT:.4f}:1)",
        f"heatmap_width_inch = {heatmap_width:.4f}",
        f"heatmap_height_inch = {heatmap_height:.4f}",
        f"actual_heatmap_width_height_ratio = {figure_aspect:.4f}:1",
        f"width_per_sample_inch = {HEATMAP_WIDTH_PER_SAMPLE:.6f}",
        f"height_per_gene_inch = {HEATMAP_HEIGHT_PER_GENE:.6f}",
        "",
        f"mutation_type_width_inch = {mutation_type_width:.4f}",
        f"mutation_type_height_inch = {mutation_type_height:.4f}",
        f"mutation_type_width_height_ratio = {mutation_type_aspect:.4f}:1",
        f"mutation_type_height_per_gene_inch = {MUTATION_TYPE_HEIGHT_PER_GENE:.6f}",
        "",
        f"total_maf_rows = {total_maf_rows}",
        f"filtered_maf_rows = {filtered_maf_rows}",
        f"found_mutations = {', '.join(found_mutations)}",
        "",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")


def _write_v2_result_log(
    output_path: Path,
    cohort: str,
    mutation_binary_path: Path,
    meta_path: Path,
    maf_path: Path,
    gene_symbol_map_path: Path,
    species: str,
    heatmap_pdf: Path,
    mutation_type_pdf: Path,
    gene_summary_tsv: Path,
    binary_matrix_tsv: Path,
    mutation_type_matrix_tsv: Path,
    encoded_matrix_tsv: Path,
    sample_annotations_tsv: Path,
    sample_annotation_colors_tsv: Path,
    mutation_color_table_tsv: Path,
    filtered_maf_tsv: Path,
    sample_count: int,
    gene_count: int,
    heatmap_width: float,
    heatmap_height: float,
    mutation_type_width: float,
    mutation_type_height: float,
    total_maf_rows: int,
    filtered_maf_rows: int,
    found_mutations: list[str],
) -> None:
    sample_gene_ratio = sample_count / gene_count if gene_count else float("inf")
    figure_aspect = heatmap_width / heatmap_height if heatmap_height else float("inf")
    mutation_type_aspect = mutation_type_width / mutation_type_height if mutation_type_height else float("inf")
    lines = [
        "OmicsOne mutation figure result v2",
        "",
        f"mutation_binary_path = {mutation_binary_path}",
        f"meta_path = {meta_path}",
        f"maf_path = {maf_path}",
        f"gene_symbol_map_path = {gene_symbol_map_path}",
        f"species = {species}",
        f"heatmap_pdf = {heatmap_pdf}",
        f"mutation_type_pdf = {mutation_type_pdf}",
        f"gene_summary_tsv = {gene_summary_tsv}",
        f"binary_matrix_tsv = {binary_matrix_tsv}",
        f"mutation_type_matrix_tsv = {mutation_type_matrix_tsv}",
        f"encoded_matrix_tsv = {encoded_matrix_tsv}",
        f"sample_annotations_tsv = {sample_annotations_tsv}",
        f"sample_annotation_colors_tsv = {sample_annotation_colors_tsv}",
        f"mutation_color_table_tsv = {mutation_color_table_tsv}",
        f"filtered_maf_tsv = {filtered_maf_tsv}",
        "",
        f"cohort = {cohort}",
        f"sample_count = {sample_count}",
        f"gene_count = {gene_count}",
        f"sample_gene_ratio = {sample_count}:{gene_count} ({sample_gene_ratio:.4f}:1)",
        "",
        f"target_heatmap_width_height_ratio = {TARGET_HEATMAP_WIDTH_RATIO}:{TARGET_HEATMAP_HEIGHT_RATIO} ({TARGET_HEATMAP_ASPECT:.4f}:1)",
        f"heatmap_width_inch = {heatmap_width:.4f}",
        f"heatmap_height_inch = {heatmap_height:.4f}",
        f"actual_heatmap_width_height_ratio = {figure_aspect:.4f}:1",
        f"width_per_sample_inch = {HEATMAP_WIDTH_PER_SAMPLE:.6f}",
        f"height_per_gene_inch = {HEATMAP_HEIGHT_PER_GENE:.6f}",
        "",
        f"mutation_type_width_inch = {mutation_type_width:.4f}",
        f"mutation_type_height_inch = {mutation_type_height:.4f}",
        f"mutation_type_width_height_ratio = {mutation_type_aspect:.4f}:1",
        f"mutation_type_height_per_gene_inch = {MUTATION_TYPE_HEIGHT_PER_GENE:.6f}",
        "",
        f"total_maf_rows = {total_maf_rows}",
        f"filtered_maf_rows = {filtered_maf_rows}",
        f"found_mutations = {', '.join(found_mutations)}",
        "",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")


def _write_intermediate_tables(
    output_dir: Path,
    cohort: str,
    encoded_df: pd.DataFrame,
    typed_df: pd.DataFrame,
    binary_df: pd.DataFrame,
    col_df: pd.DataFrame,
    col_color_df: pd.DataFrame,
    mutation_color_df: pd.DataFrame,
    filtered_maf: pd.DataFrame,
    gene_summary_df: pd.DataFrame,
) -> dict[str, Path]:
    paths = {
        "gene_summary_tsv": output_dir / f"{cohort}_plotted_gene_summary.tsv",
        "binary_matrix_tsv": output_dir / f"{cohort}_plotted_binary_matrix.tsv",
        "mutation_type_matrix_tsv": output_dir / f"{cohort}_plotted_mutation_type_matrix.tsv",
        "encoded_matrix_tsv": output_dir / f"{cohort}_plotted_encoded_matrix.tsv",
        "sample_annotations_tsv": output_dir / f"{cohort}_plotted_sample_annotations.tsv",
        "sample_annotation_colors_tsv": output_dir / f"{cohort}_plotted_sample_annotation_colors.tsv",
        "mutation_color_table_tsv": output_dir / f"{cohort}_mutation_color_table.tsv",
        "filtered_maf_tsv": output_dir / f"{cohort}_filtered_maf_for_plotted_genes.tsv",
    }
    gene_summary_df.to_csv(paths["gene_summary_tsv"], sep="\t", index=False)
    binary_df.to_csv(paths["binary_matrix_tsv"], sep="\t")
    typed_df.to_csv(paths["mutation_type_matrix_tsv"], sep="\t")
    encoded_df.to_csv(paths["encoded_matrix_tsv"], sep="\t")
    col_df.to_csv(paths["sample_annotations_tsv"], sep="\t")
    col_color_df.to_csv(paths["sample_annotation_colors_tsv"], sep="\t")
    mutation_color_df.to_csv(paths["mutation_color_table_tsv"], sep="\t", index=False)
    filtered_maf.to_csv(paths["filtered_maf_tsv"], sep="\t", index=False)
    return paths


def _type_distribution(row: pd.Series) -> pd.Series:
    non_zero = row[row != 0]
    if len(non_zero) == 0:
        return pd.Series(dtype=float)
    return non_zero.value_counts(normalize=True) * 100


def _plot_mutation_type_distribution(
    typed_df: pd.DataFrame,
    mutation_color_df: pd.DataFrame,
    output_path: Path,
    figure_size: tuple[float, float],
) -> None:
    reversed_df = typed_df.iloc[::-1].copy()
    distribution = reversed_df.apply(_type_distribution, axis=1).fillna(0)
    distribution.columns = distribution.columns.astype(str)

    color_map = dict(zip(mutation_color_df["mut"], mutation_color_df["color"]))
    fig, ax = plt.subplots(figsize=figure_size)
    distribution.plot(kind="barh", stacked=True, color=color_map, ax=ax)
    ax.set_xlim(0, 100)
    ax.set_ylabel("Gene")
    ax.set_xlabel("Percentage (%)")
    ax.set_title("Percentage of Mutations for Genes")
    ax.legend(title="Mutations", bbox_to_anchor=(1.05, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, format="pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def generate_hnsc_mutation_figures(
    mutation_excel_path: str | Path = DEFAULT_MUTATION_EXCEL_PATH,
    maf_path: str | Path = DEFAULT_MAF_PATH,
    output_dir: str | Path = DEFAULT_OUTPUT_DIR,
    mutation_threshold: float = 0.15,
    cohort: str = DEFAULT_COHORT,
    heatmap_filename: Optional[str] = None,
    mutation_type_filename: Optional[str] = None,
    output_prefix: Optional[str] = None,
) -> MutationFigureResult:
    mutation_excel = _ensure_file(mutation_excel_path, "Mutation Excel file")
    maf_file = _ensure_file(maf_path, "MAF file")
    out_dir = Path(output_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    cohort_name = _normalize_cohort(cohort)

    if heatmap_filename is None:
        heatmap_filename = f"{cohort_name}_mutation_heatmap.pdf"
    if mutation_type_filename is None:
        mutation_type_filename = f"{cohort_name}_mutation_type_distribution.pdf"

    if output_prefix:
        heatmap_filename = f"{output_prefix}_{heatmap_filename}"
        mutation_type_filename = f"{output_prefix}_{mutation_type_filename}"

    heatmap_pdf = out_dir / heatmap_filename
    mutation_type_pdf = out_dir / mutation_type_filename
    result_log = out_dir / "result.log"

    mutation_df = pd.read_excel(mutation_excel)
    maf_df = pd.read_csv(maf_file, sep="\t", comment="#")

    (
        encoded_df,
        typed_df,
        binary_df,
        col_df,
        col_color_df,
        mutation_color_df,
        filtered_maf,
        gene_summary_df,
        found_mutations,
    ) = _prepare_mutation_tables(
        mutation_df,
        maf_df,
        mutation_threshold,
    )
    sample_count = encoded_df.shape[1]
    gene_count = encoded_df.shape[0]
    heatmap_width, heatmap_height = _calculate_heatmap_figure_size(sample_count, gene_count)
    mutation_type_width, mutation_type_height = _calculate_mutation_type_figure_size(gene_count)
    intermediate_paths = _write_intermediate_tables(
        out_dir,
        cohort_name,
        encoded_df,
        typed_df,
        binary_df,
        col_df,
        col_color_df,
        mutation_color_df,
        filtered_maf,
        gene_summary_df,
    )

    _plot_mutation_heatmap(
        encoded_df,
        col_color_df,
        mutation_color_df,
        found_mutations,
        heatmap_pdf,
        (heatmap_width, heatmap_height),
    )
    _plot_mutation_type_distribution(
        typed_df,
        mutation_color_df,
        mutation_type_pdf,
        (mutation_type_width, mutation_type_height),
    )
    _write_result_log(
        result_log,
        cohort_name,
        mutation_excel,
        maf_file,
        heatmap_pdf,
        mutation_type_pdf,
        intermediate_paths["gene_summary_tsv"],
        intermediate_paths["binary_matrix_tsv"],
        intermediate_paths["mutation_type_matrix_tsv"],
        intermediate_paths["encoded_matrix_tsv"],
        intermediate_paths["sample_annotations_tsv"],
        intermediate_paths["sample_annotation_colors_tsv"],
        intermediate_paths["mutation_color_table_tsv"],
        intermediate_paths["filtered_maf_tsv"],
        sample_count,
        gene_count,
        heatmap_width,
        heatmap_height,
        mutation_type_width,
        mutation_type_height,
        maf_df.shape[0],
        filtered_maf.shape[0],
        found_mutations,
    )

    return MutationFigureResult(
        heatmap_pdf=heatmap_pdf,
        mutation_type_pdf=mutation_type_pdf,
        result_log=result_log,
        gene_summary_tsv=intermediate_paths["gene_summary_tsv"],
        binary_matrix_tsv=intermediate_paths["binary_matrix_tsv"],
        mutation_type_matrix_tsv=intermediate_paths["mutation_type_matrix_tsv"],
        encoded_matrix_tsv=intermediate_paths["encoded_matrix_tsv"],
        sample_annotations_tsv=intermediate_paths["sample_annotations_tsv"],
        sample_annotation_colors_tsv=intermediate_paths["sample_annotation_colors_tsv"],
        mutation_color_table_tsv=intermediate_paths["mutation_color_table_tsv"],
        filtered_maf_tsv=intermediate_paths["filtered_maf_tsv"],
        gene_count=gene_count,
        sample_count=sample_count,
        cohort=cohort_name,
        total_maf_rows=maf_df.shape[0],
        filtered_maf_rows=filtered_maf.shape[0],
        found_mutations=found_mutations,
        heatmap_width_inch=heatmap_width,
        heatmap_height_inch=heatmap_height,
        heatmap_aspect=heatmap_width / heatmap_height,
        mutation_type_width_inch=mutation_type_width,
        mutation_type_height_inch=mutation_type_height,
        mutation_type_aspect=mutation_type_width / mutation_type_height,
        sample_gene_ratio=sample_count / gene_count,
    )


def generate_mutation_figures_from_binary(
    mutation_binary_path: str | Path = DEFAULT_MUTATION_BINARY_PATH,
    meta_path: str | Path = DEFAULT_MUTATION_META_PATH,
    maf_path: str | Path = DEFAULT_V2_MAF_PATH,
    output_dir: str | Path = DEFAULT_V2_OUTPUT_DIR,
    mutation_threshold: float = 0.15,
    cohort: str = DEFAULT_V2_COHORT,
    species: str = "human",
    gene_symbol_map_path: Optional[str | Path] = None,
    heatmap_filename: Optional[str] = None,
    mutation_type_filename: Optional[str] = None,
    output_prefix: Optional[str] = None,
) -> MutationFigureV2Result:
    mutation_binary = _ensure_file(mutation_binary_path, "Mutation binary file")
    meta_file = _ensure_file(meta_path, "Meta file")
    maf_file = _ensure_file(maf_path, "MAF file")
    species_name = _normalize_species(species)
    resolved_gene_symbol_map = _ensure_file(
        gene_symbol_map_path or _default_gene_symbol_map_path(species_name),
        "Gene symbol map file",
    )
    gene_symbol_map = _load_gene_symbol_map(resolved_gene_symbol_map)

    out_dir = Path(output_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    cohort_name = _normalize_cohort(cohort)

    if heatmap_filename is None:
        heatmap_filename = f"{cohort_name}_mutation_heatmap.pdf"
    if mutation_type_filename is None:
        mutation_type_filename = f"{cohort_name}_mutation_type_distribution.pdf"

    if output_prefix:
        heatmap_filename = f"{output_prefix}_{heatmap_filename}"
        mutation_type_filename = f"{output_prefix}_{mutation_type_filename}"

    heatmap_pdf = out_dir / heatmap_filename
    mutation_type_pdf = out_dir / mutation_type_filename
    result_log = out_dir / "result.log"

    mutation_binary_df = pd.read_csv(mutation_binary, sep="\t", index_col=0)
    meta_df = pd.read_csv(meta_file, sep="\t", index_col=0, header=[0, 1])
    maf_df = pd.read_csv(maf_file, sep="\t", comment="#")

    (
        encoded_df,
        typed_df,
        binary_df,
        col_df,
        col_color_df,
        mutation_color_df,
        filtered_maf,
        gene_summary_df,
        found_mutations,
    ) = _prepare_mutation_tables_from_binary(
        mutation_binary_df,
        meta_df,
        maf_df,
        gene_symbol_map,
        mutation_threshold,
    )
    sample_count = encoded_df.shape[1]
    gene_count = encoded_df.shape[0]
    heatmap_width, heatmap_height = _calculate_heatmap_figure_size(sample_count, gene_count)
    mutation_type_width, mutation_type_height = _calculate_mutation_type_figure_size(gene_count)
    intermediate_paths = _write_intermediate_tables(
        out_dir,
        cohort_name,
        encoded_df,
        typed_df,
        binary_df,
        col_df,
        col_color_df,
        mutation_color_df,
        filtered_maf,
        gene_summary_df,
    )

    _plot_mutation_heatmap(
        encoded_df,
        col_color_df,
        mutation_color_df,
        found_mutations,
        heatmap_pdf,
        (heatmap_width, heatmap_height),
    )
    _plot_mutation_type_distribution(
        typed_df,
        mutation_color_df,
        mutation_type_pdf,
        (mutation_type_width, mutation_type_height),
    )
    _write_v2_result_log(
        result_log,
        cohort_name,
        mutation_binary,
        meta_file,
        maf_file,
        resolved_gene_symbol_map,
        species_name,
        heatmap_pdf,
        mutation_type_pdf,
        intermediate_paths["gene_summary_tsv"],
        intermediate_paths["binary_matrix_tsv"],
        intermediate_paths["mutation_type_matrix_tsv"],
        intermediate_paths["encoded_matrix_tsv"],
        intermediate_paths["sample_annotations_tsv"],
        intermediate_paths["sample_annotation_colors_tsv"],
        intermediate_paths["mutation_color_table_tsv"],
        intermediate_paths["filtered_maf_tsv"],
        sample_count,
        gene_count,
        heatmap_width,
        heatmap_height,
        mutation_type_width,
        mutation_type_height,
        maf_df.shape[0],
        filtered_maf.shape[0],
        found_mutations,
    )

    return MutationFigureV2Result(
        heatmap_pdf=heatmap_pdf,
        mutation_type_pdf=mutation_type_pdf,
        result_log=result_log,
        gene_summary_tsv=intermediate_paths["gene_summary_tsv"],
        binary_matrix_tsv=intermediate_paths["binary_matrix_tsv"],
        mutation_type_matrix_tsv=intermediate_paths["mutation_type_matrix_tsv"],
        encoded_matrix_tsv=intermediate_paths["encoded_matrix_tsv"],
        sample_annotations_tsv=intermediate_paths["sample_annotations_tsv"],
        sample_annotation_colors_tsv=intermediate_paths["sample_annotation_colors_tsv"],
        mutation_color_table_tsv=intermediate_paths["mutation_color_table_tsv"],
        filtered_maf_tsv=intermediate_paths["filtered_maf_tsv"],
        gene_count=gene_count,
        sample_count=sample_count,
        cohort=cohort_name,
        total_maf_rows=maf_df.shape[0],
        filtered_maf_rows=filtered_maf.shape[0],
        found_mutations=found_mutations,
        heatmap_width_inch=heatmap_width,
        heatmap_height_inch=heatmap_height,
        heatmap_aspect=heatmap_width / heatmap_height,
        mutation_type_width_inch=mutation_type_width,
        mutation_type_height_inch=mutation_type_height,
        mutation_type_aspect=mutation_type_width / mutation_type_height,
        sample_gene_ratio=sample_count / gene_count,
        mutation_binary_path=mutation_binary,
        meta_path=meta_file,
        gene_symbol_map_path=resolved_gene_symbol_map,
        species=species_name,
    )
