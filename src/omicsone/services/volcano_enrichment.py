from __future__ import annotations

import json
import math
import os
import configparser
import warnings
from html import escape
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import gseapy as gp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pyteomics import fasta
from scipy.stats import hypergeom, mannwhitneyu, ttest_ind, ttest_rel, wilcoxon
from statsmodels.stats.multitest import multipletests


DEFAULT_FASTA_PATH = (
    r"F:\lab\HsinI\Head and Neck & Lung\fasta"
    r"\GENCODE.V42.basic.CHR.combined_contaminants.gpquest3.fasta"
)
DEFAULT_BASE_OUTPUT_DIR = r"E:\lab\HSinI\runs"
DEFAULT_API_URL = "http://127.0.0.1:8001/api/v1/diff/volcano/enrichment"
DEFAULT_GSEAPY_CACHE_DIR = Path.home() / ".cache" / "gseapy"
DEFAULT_ENRICHMENT_WIDTH_RATIO = 2.0
DEFAULT_ENRICHMENT_HEIGHT_RATIO = 1.5
DEFAULT_ENRICHMENT_SIZE_SCALE = 4.0
ENRICHMENT_BACKGROUND_MODES = {"gene_list", "online", "notebook", "count"}

DEFAULT_SKIP_PATHWAYS = [
    "Phagosome",
    "Human papillomavirus infection",
    "Pertussis",
    "Malaria",
    "Arrhythmogenic right ventricular cardiomyopathy",
    "Staphylococcus aureus infection",
    "Regulation of actin cytoskeleton",
]

DEFAULT_PRESETS: dict[str, dict[str, Any]] = {
    "HNSCC_RNA": {
        "cohort": "HNSCC",
        "omics": "RNA",
        "normal_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\HNSCC"
            r"\HNSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt"
        ),
        "tumor_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\HNSCC"
            r"\HNSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"
        ),
        "output_dir": rf"{DEFAULT_BASE_OUTPUT_DIR}\20260503_HNSCC_RNA_volcano",
        "title": "HNSCC RNA differential expression analysis",
        "enrichment_title": "HNSCC RNA MSigDB hallmark Pathways Enrichment Analysis",
        "enrichment_min_x": "-25",
        "enrichment_max_x": "60",
    },
    "HNSCC_Protein": {
        "cohort": "HNSCC",
        "omics": "Protein",
        "normal_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\HNSCC"
            r"\HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt"
        ),
        "tumor_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\HNSCC"
            r"\HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
        ),
        "output_dir": rf"{DEFAULT_BASE_OUTPUT_DIR}\20260503_HNSCC_Protein_volcano",
        "title": "HNSCC protein differential expression analysis",
        "enrichment_title": "HNSCC protein MSigDB hallmark Pathways Enrichment Analysis",
        "enrichment_min_x": "-25",
        "enrichment_max_x": "40",
    },
    "HNSCC_Phospho": {
        "cohort": "HNSCC",
        "omics": "Phospho",
        "normal_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\HNSCC"
            r"\HNSCC_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt"
        ),
        "tumor_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\HNSCC"
            r"\HNSCC_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt"
        ),
        "output_dir": rf"{DEFAULT_BASE_OUTPUT_DIR}\20260517_HNSCC_phospho",
        "title": "IGP-based differential expression analysis",
        "enrichment_title": "HNSCC phosphosites MSigDB hallmark Pathways Enrichment Analysis",
        "enrichment_min_x": "-25",
        "enrichment_max_x": "40",
        "strip_feature_version": False,
        "output_prefix": "hnscc_phospho",
    },
    "LSCC_RNA": {
        "cohort": "LSCC",
        "omics": "RNA",
        "normal_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\LSCC"
            r"\LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt"
        ),
        "tumor_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\LSCC"
            r"\LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"
        ),
        "output_dir": rf"{DEFAULT_BASE_OUTPUT_DIR}\20260503_LSCC_RNA_volcano",
        "title": "LSCC RNA differential expression analysis",
        "enrichment_title": "LSCC RNA MSigDB hallmark Pathways Enrichment Analysis",
        "enrichment_min_x": "-25",
        "enrichment_max_x": "60",
    },
    "LSCC_Protein": {
        "cohort": "LSCC",
        "omics": "Protein",
        "normal_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\LSCC"
            r"\LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt"
        ),
        "tumor_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\LSCC"
            r"\LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
        ),
        "output_dir": rf"{DEFAULT_BASE_OUTPUT_DIR}\20260503_LSCC_Protein_volcano",
        "title": "LSCC protein differential expression analysis",
        "enrichment_title": "LSCC protein MSigDB hallmark Pathways Enrichment Analysis",
        "enrichment_min_x": "-25",
        "enrichment_max_x": "40",
    },
    "LSCC_Phospho": {
        "cohort": "LSCC",
        "omics": "Phospho",
        "normal_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\LSCC"
            r"\LSCC_phospho_site_abundance_log2_reference_intensity_normalized_Normal.txt"
        ),
        "tumor_path": (
            r"F:\lab\HsinI\Head and Neck & Lung\LSCC"
            r"\LSCC_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt"
        ),
        "output_dir": rf"{DEFAULT_BASE_OUTPUT_DIR}\20260517_LSCC_phospho",
        "title": "LSCC phosphosite differential abundance analysis",
        "enrichment_title": "LSCC phosphosites MSigDB hallmark Pathways Enrichment Analysis",
        "enrichment_min_x": "-25",
        "enrichment_max_x": "40",
        "strip_feature_version": False,
        "output_prefix": "lscc_phospho",
    },
}


@dataclass(frozen=True)
class VolcanoEnrichmentResult:
    cohort: str
    omics: str
    output_dir: Path
    diff_tsv: Path
    combined_matrix_tsv: Path
    up_genes_tsv: Path
    down_genes_tsv: Path
    total_genes_tsv: Path
    up_enrichment_tsv: Path
    down_enrichment_tsv: Path
    enrichment_plot_tsv: Path
    volcano_png: Path
    volcano_pdf: Path
    volcano_tiff: Path
    enrichment_png: Path
    enrichment_pdf: Path
    enrichment_tiff: Path
    report_html: Path
    result_log: Path
    n8n_js: Path
    feature_count: int
    diff_feature_count: int
    up_count: int
    down_count: int
    pure_up_gene_count: int
    pure_down_gene_count: int
    total_gene_count: int
    up_enrichment_count: int
    down_enrichment_count: int
    method: str
    fdr_cutoff: float
    log2fc_cutoff: float
    gene_sets: list[str]


def resolve_volcano_preset(name: str) -> dict[str, Any]:
    try:
        preset = DEFAULT_PRESETS[name]
    except KeyError as exc:
        raise ValueError(f"Unknown volcano preset: {name}") from exc

    resolved: dict[str, Any] = {
        "job_name": name,
        "cohort": preset["cohort"],
        "omics": preset["omics"],
        "normal_path": preset["normal_path"],
        "tumor_path": preset["tumor_path"],
        "output_dir": preset["output_dir"],
        "fasta_path": DEFAULT_FASTA_PATH,
        "title": preset["title"],
        "enrichment_title": preset["enrichment_title"],
        "enrichment_min_x": float(preset["enrichment_min_x"]),
        "enrichment_max_x": float(preset["enrichment_max_x"]),
        "strip_feature_version": bool(preset.get("strip_feature_version", True)),
    }
    if preset.get("output_prefix"):
        resolved["output_prefix"] = preset["output_prefix"]
    return resolved


@lru_cache(maxsize=4)
def _read_gene_map(fasta_path: str) -> dict[str, str]:
    gene_map: dict[str, str] = {}
    for description, _sequence in fasta.read(fasta_path):
        items = description.split("|")
        gene_id = None
        gene_name = None

        for item in items:
            if item.startswith("GI="):
                gene_id = item.split("=", maxsplit=1)[1].split(".", maxsplit=1)[0]
                break
            if item.startswith("ENSG"):
                gene_id = item.split(".", maxsplit=1)[0]

        for item in items:
            if item.startswith("GN="):
                gene_name = item.split("=", maxsplit=1)[1].split(" ", maxsplit=1)[0]
                break

        if gene_id and gene_name:
            gene_map[gene_id] = gene_name

    return gene_map


def _feature_gene_id(feature_id: str) -> str:
    return str(feature_id).split("|", maxsplit=1)[0].split(".", maxsplit=1)[0]


def _strip_ensembl_version(feature_id: str) -> str:
    return str(feature_id).split(".", maxsplit=1)[0]


def _read_group_matrix(path: Path, suffix: str, strip_feature_version: bool) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input matrix does not exist: {path}")

    df = pd.read_csv(path, sep="\t")
    if "idx" in df.columns:
        df = df.set_index("idx")
    elif df.columns.size:
        df = df.set_index(df.columns[0])
    else:
        raise ValueError(f"Input matrix has no columns: {path}")

    if strip_feature_version:
        df.index = [_strip_ensembl_version(index) for index in df.index]
    else:
        df.index = [str(index) for index in df.index]
    df.columns = [f"{column}{suffix}" for column in df.columns.values]
    return df


def _map_gene_ids(gene_ids: list[str], gene_map: dict[str, str]) -> list[str]:
    mapped = []
    for gene_id in gene_ids:
        key = _feature_gene_id(str(gene_id))
        gene = gene_map.get(key)
        if gene:
            mapped.append(gene)
    return sorted(set(mapped))


def _save_gene_list(path: Path, genes: list[str]) -> None:
    pd.Series(genes, name="gene").to_csv(path, sep="\t", index=False)


@lru_cache(maxsize=16)
def _read_gmt(path: str) -> tuple[tuple[str, str, frozenset[str]], ...]:
    gene_sets = []
    gmt_path = Path(path)
    with gmt_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term, description, *genes = parts
            gene_sets.append((gmt_path.stem, term, frozenset(gene for gene in genes if gene)))
    return tuple(gene_sets)


def _local_count_enrichr(
    genes: list[str],
    total_genes: list[str],
    gene_set_paths: list[str],
    fdr_cutoff: float,
) -> pd.DataFrame:
    query = set(genes)
    background_size = len(set(total_genes))
    query_size = len(query)
    rows = []

    for gene_set_path in gene_set_paths:
        for gene_set_name, term, term_genes in _read_gmt(gene_set_path):
            overlap_genes = sorted(query & set(term_genes))
            overlap_size = len(overlap_genes)
            term_size = len(term_genes)
            if overlap_size == 0 or background_size == 0 or query_size == 0:
                continue

            pvalue = float(hypergeom.sf(overlap_size - 1, background_size, term_size, query_size))
            non_overlap_query = query_size - overlap_size
            non_overlap_term = term_size - overlap_size
            background_not_query_or_term = max(
                background_size - query_size - term_size + overlap_size,
                0,
            )
            odds_ratio = (
                (overlap_size * background_not_query_or_term)
                / max(non_overlap_query * non_overlap_term, np.finfo(float).tiny)
            )
            rows.append(
                {
                    "Gene_set": gene_set_name,
                    "Term": term,
                    "Overlap": f"{overlap_size}/{term_size}",
                    "P-value": pvalue,
                    "Old P-value": pvalue,
                    "Old Adjusted P-value": pvalue,
                    "Odds Ratio": odds_ratio,
                    "Combined Score": -math.log(max(pvalue, np.finfo(float).tiny)) * odds_ratio,
                    "Genes": ";".join(overlap_genes),
                }
            )

    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(
            columns=[
                "Gene_set",
                "Term",
                "Overlap",
                "P-value",
                "Adjusted P-value",
                "Old P-value",
                "Old Adjusted P-value",
                "Odds Ratio",
                "Combined Score",
                "Genes",
                "-Log10(Adj.P)",
            ]
        )

    df["Adjusted P-value"] = multipletests(df["P-value"], method="fdr_bh")[1]
    df = df[df["Adjusted P-value"] < fdr_cutoff].copy()
    if df.empty:
        df["-Log10(Adj.P)"] = []
        return df

    df["-Log10(Adj.P)"] = df["Adjusted P-value"].apply(lambda value: -math.log10(value))
    columns = [
        "Gene_set",
        "Term",
        "Overlap",
        "P-value",
        "Adjusted P-value",
        "Old P-value",
        "Old Adjusted P-value",
        "Odds Ratio",
        "Combined Score",
        "Genes",
        "-Log10(Adj.P)",
    ]
    return df.loc[:, columns].sort_values("-Log10(Adj.P)", ascending=False)


def _decide_sig(fdr: float, log2fc: float, log2fc_cutoff: float, fdr_cutoff: float) -> str | None:
    if log2fc > log2fc_cutoff and fdr < fdr_cutoff:
        return "S-U"
    if log2fc < -1 * log2fc_cutoff and fdr < fdr_cutoff:
        return "S-D"
    if log2fc > 0 and fdr < fdr_cutoff:
        return "U"
    if log2fc < 0 and fdr < fdr_cutoff:
        return "D"
    return None


def _compare_two_groups_fast(
    df: pd.DataFrame,
    group_a: list[str],
    group_b: list[str],
    method: str,
    fdr_cutoff: float,
    log2fc_cutoff: float,
    max_miss_ratio_global: float,
    max_miss_ratio_group: float,
    min_sample_size: int,
) -> pd.DataFrame:
    methods = {
        "T-test(Unpaired)": ttest_ind,
        "T-test(Paired)": ttest_rel,
        "Wilcoxon(Unpaired)": mannwhitneyu,
        "Wilcoxon(Paired)": wilcoxon,
    }
    if method not in methods:
        raise ValueError(f"Unsupported differential method: {method}")

    sample_columns = group_a + group_b
    matrix = df.loc[:, sample_columns].apply(pd.to_numeric, errors="coerce")
    a_values = matrix.loc[:, group_a].to_numpy(dtype=float)
    b_values = matrix.loc[:, group_b].to_numpy(dtype=float)

    miss_global = np.isnan(matrix.to_numpy(dtype=float)).mean(axis=1)
    valid = miss_global <= max_miss_ratio_global

    if method in {"T-test(Paired)", "Wilcoxon(Paired)"}:
        if len(group_a) != len(group_b):
            raise ValueError(f"{method} requires groups with the same sample count.")
        pair_mask = ~np.isnan(a_values) & ~np.isnan(b_values)
        valid_pairs = pair_mask.sum(axis=1)
        valid &= valid_pairs >= min_sample_size
        diff_values = np.where(pair_mask, a_values - b_values, np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            log2fc_median = np.nanmedian(diff_values, axis=1)
            log2fc_mean = np.nanmean(diff_values, axis=1)
    else:
        a_non_missing = (~np.isnan(a_values)).sum(axis=1)
        b_non_missing = (~np.isnan(b_values)).sum(axis=1)
        valid &= (1 - a_non_missing / len(group_a)) <= max_miss_ratio_group
        valid &= (1 - b_non_missing / len(group_b)) <= max_miss_ratio_group
        valid &= a_non_missing >= min_sample_size
        valid &= b_non_missing >= min_sample_size
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            log2fc_median = np.nanmedian(a_values, axis=1) - np.nanmedian(b_values, axis=1)
            log2fc_mean = np.nanmean(a_values, axis=1) - np.nanmean(b_values, axis=1)

    rows = []
    test_fn = methods[method]
    stat_column = f"{method}(Stats)"
    pvalue_column = f"{method}(P-value)"

    for row_pos in np.flatnonzero(valid):
        if method in {"T-test(Paired)", "Wilcoxon(Paired)"}:
            pair_valid = ~np.isnan(a_values[row_pos]) & ~np.isnan(b_values[row_pos])
            x = a_values[row_pos, pair_valid]
            y = b_values[row_pos, pair_valid]
        else:
            x = a_values[row_pos, ~np.isnan(a_values[row_pos])]
            y = b_values[row_pos, ~np.isnan(b_values[row_pos])]

        try:
            if method == "T-test(Unpaired)":
                test_result = test_fn(x, y, nan_policy="omit")
            else:
                test_result = test_fn(x, y)
            stat = float(test_result.statistic)
            pvalue = float(test_result.pvalue)
        except ValueError:
            continue

        if not np.isfinite(pvalue):
            continue
        rows.append(
            [
                matrix.index[row_pos],
                float(log2fc_median[row_pos]),
                float(log2fc_mean[row_pos]),
                stat,
                pvalue,
            ]
        )

    diff = pd.DataFrame(
        rows,
        columns=["Feature", "Log2FC(median)", "Log2FC(mean)", stat_column, pvalue_column],
    )
    if diff.empty:
        return diff.set_index("Feature")

    diff["FDR"] = multipletests(diff[pvalue_column], method="fdr_bh")[1]
    diff["-Log10(FDR)"] = -np.log10(diff["FDR"].clip(lower=np.finfo(float).tiny))
    diff["Significance"] = [
        _decide_sig(row["FDR"], row["Log2FC(median)"], log2fc_cutoff, fdr_cutoff)
        for _index, row in diff.iterrows()
    ]
    return diff.set_index("Feature")


def _enrichr(
    genes: list[str],
    total_genes: list[str],
    output_dir: Path,
    gene_sets: list[str],
    organism: str,
    fdr_cutoff: float,
    enrichment_background_mode: str,
    prefer_local_gene_sets: bool,
    allow_remote_enrichr: bool,
) -> pd.DataFrame:
    columns = [
        "Gene_set",
        "Term",
        "Overlap",
        "P-value",
        "Adjusted P-value",
        "Old P-value",
        "Old Adjusted P-value",
        "Odds Ratio",
        "Combined Score",
        "Genes",
        "-Log10(Adj.P)",
    ]
    if not genes:
        return pd.DataFrame(columns=columns)

    if enrichment_background_mode not in ENRICHMENT_BACKGROUND_MODES:
        raise ValueError(
            "enrichment_background_mode must be one of: "
            f"{', '.join(sorted(ENRICHMENT_BACKGROUND_MODES))}."
        )

    resolved_gene_sets: list[str] = []
    missing_gene_sets: list[str] = []
    for gene_set in gene_sets:
        gene_set_path = Path(os.path.expanduser(gene_set))
        cached_gene_set_path = DEFAULT_GSEAPY_CACHE_DIR / f"Enrichr.{gene_set}.gmt"
        if gene_set_path.exists():
            resolved_gene_sets.append(str(gene_set_path))
        elif prefer_local_gene_sets and cached_gene_set_path.exists():
            resolved_gene_sets.append(str(cached_gene_set_path))
        else:
            missing_gene_sets.append(gene_set)

    if enrichment_background_mode in {"online", "notebook", "count"}:
        if missing_gene_sets and not allow_remote_enrichr:
            raise ValueError(
                f"enrichment_background_mode='{enrichment_background_mode}' with unresolved "
                "Enrichr library names would use online Enrichr. Pass local GMT paths, install "
                f"them in {DEFAULT_GSEAPY_CACHE_DIR}, or set allow_remote_enrichr=true."
            )
        effective_gene_sets = resolved_gene_sets + missing_gene_sets
        if not missing_gene_sets:
            return _local_count_enrichr(
                genes=genes,
                total_genes=total_genes,
                gene_set_paths=effective_gene_sets,
                fdr_cutoff=fdr_cutoff,
            )
        # Notebook-compatible mode:
        # enrichr(genes, total_genes, outdir, gene_sets=[...]) delegates to
        # gp.enrichr(..., background=len(total_genes), gene_sets=gene_sets, ...).
        enr = gp.enrichr(
            gene_list=list(genes),
            background=len(total_genes),
            gene_sets=effective_gene_sets,
            organism=organism,
            outdir=str(output_dir),
        )
        df = enr.res2d
        if df.empty:
            return pd.DataFrame(columns=columns)

        df = df[df["Adjusted P-value"] < fdr_cutoff].copy()
        if df.empty:
            df["-Log10(Adj.P)"] = []
            return df

        df["-Log10(Adj.P)"] = df["Adjusted P-value"].apply(
            lambda value: -math.log10(value)
        )
        return df.sort_values("-Log10(Adj.P)", ascending=False)

    if missing_gene_sets and not allow_remote_enrichr:
        missing = ", ".join(missing_gene_sets)
        raise ValueError(
            "Local gene set GMT file not found for "
            f"{missing}. Pass a local GMT path, install it in "
            f"{DEFAULT_GSEAPY_CACHE_DIR}, or set allow_remote_enrichr=true "
            "to permit Enrichr network submission."
        )

    if missing_gene_sets:
        resolved_gene_sets.extend(missing_gene_sets)

    enr = gp.enrichr(
        gene_list=list(genes),
        background=list(total_genes) if resolved_gene_sets else len(total_genes),
        gene_sets=resolved_gene_sets,
        organism=organism,
        outdir=str(output_dir),
    )
    df = enr.res2d
    if df.empty:
        return pd.DataFrame(columns=columns)

    df = df[df["Adjusted P-value"] < fdr_cutoff].copy()
    if df.empty:
        df["-Log10(Adj.P)"] = []
        return df

    df["-Log10(Adj.P)"] = df["Adjusted P-value"].apply(lambda value: -math.log10(value))
    return df.sort_values("-Log10(Adj.P)", ascending=False)


def _plot_volcano(
    diff: pd.DataFrame,
    title: str,
    xlabel: str,
    width: float,
    height: float,
    dpi: int,
    colors: dict[str, str],
    point_size: float,
    significant_point_size: float,
    fdr_cutoff: float,
    log2fc_cutoff: float,
    notebook_style: bool = False,
    show_title: bool = True,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    sns.set_style("white")

    up_df = diff[diff["Significance"] == "S-U"]
    down_df = diff[diff["Significance"] == "S-D"]
    sns.scatterplot(
        data=diff,
        x="Log2FC(median)",
        y="-Log10(FDR)",
        color=colors["background"],
        s=point_size,
        ax=ax,
    )
    sns.scatterplot(
        data=up_df,
        x="Log2FC(median)",
        y="-Log10(FDR)",
        color=colors["up"],
        s=significant_point_size,
        label=f"UP({up_df.shape[0]})",
        ax=ax,
    )
    sns.scatterplot(
        data=down_df,
        x="Log2FC(median)",
        y="-Log10(FDR)",
        color=colors["down"],
        s=significant_point_size,
        label=f"DOWN({down_df.shape[0]})",
        ax=ax,
    )

    y_upper = float(np.nanmax(diff["-Log10(FDR)"])) * 1.1 if not diff.empty else 1.0
    if notebook_style:
        plt.plot([0, 0], [0, y_upper], color="k", linewidth=0.5, linestyle="--")
        plt.plot([-log2fc_cutoff, -log2fc_cutoff], [0, y_upper], color="k", linewidth=0.5, linestyle="--")
        plt.plot([log2fc_cutoff, log2fc_cutoff], [0, y_upper], color="k", linewidth=0.5, linestyle="--")
        plt.plot([-5, 5], [-math.log10(fdr_cutoff), -math.log10(fdr_cutoff)], color="k", linewidth=0.5, linestyle="--")
    else:
        x_bound = float(
            1.1
            * max(
                abs(np.nanmin(diff["Log2FC(median)"])),
                abs(np.nanmax(diff["Log2FC(median)"])),
                log2fc_cutoff,
            )
        )
        ax.set_xlim(-x_bound, x_bound)
        ax.set_ylim(0, y_upper)
        ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
        ax.axvline(-log2fc_cutoff, color="black", linewidth=0.5, linestyle="--")
        ax.axvline(log2fc_cutoff, color="black", linewidth=0.5, linestyle="--")
        ax.axhline(-math.log10(fdr_cutoff), color="black", linewidth=0.5, linestyle="--")
    ax.legend()
    if show_title:
        ax.set_title(title)
    ax.set_xlabel(xlabel)
    if not notebook_style:
        fig.tight_layout()
    return fig


def _plot_enrichment_bars(
    up_df: pd.DataFrame,
    down_df: pd.DataFrame,
    title: str,
    min_x: float,
    max_x: float,
    width: float,
    height: float,
    dpi: int,
    colors: dict[str, str],
    notebook_style: bool = False,
    show_title: bool = True,
) -> tuple[plt.Figure, pd.DataFrame]:
    rows = []
    for _index, row in down_df.iterrows():
        rows.append([f"{row['Term']} ", -float(row["-Log10(Adj.P)"]), "Down"])
    for _index, row in up_df.iterrows():
        rows.append([f" {row['Term']}", float(row["-Log10(Adj.P)"]), "Up"])

    plot_df = pd.DataFrame(rows, columns=["Term", "-Log10(FDR)", "Class"])
    if plot_df.empty:
        fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
        ax.set_title(title)
        ax.text(0.5, 0.5, "No enriched pathways", ha="center", va="center")
        ax.axis("off")
        fig.tight_layout()
        return fig, plot_df

    plot_df = plot_df.sort_values("-Log10(FDR)", ascending=False)
    if not notebook_style:
        min_x = min(min_x, float(plot_df["-Log10(FDR)"].min()) - 1)
        max_x = max(max_x, float(plot_df["-Log10(FDR)"].max()) + 1)

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    sns.set_style("white")
    sns.barplot(
        data=plot_df,
        x="-Log10(FDR)",
        y="Term",
        hue="Class",
        palette={"Up": colors["up"], "Down": colors["down"]},
        ax=ax,
    )

    for n, (_index, row) in enumerate(plot_df.iterrows()):
        if row["-Log10(FDR)"] > 0:
            ax.text(-0.5, n + 0.5, row["Term"], ha="right", va="bottom")
        else:
            ax.text(0.5, n + 0.5, row["Term"], ha="left", va="bottom")

    ax.set_yticks([])
    ax.legend(bbox_to_anchor=(1, 0.85))
    ax.set_xlim(min_x, max_x)
    if show_title:
        ax.set_title(title)
    ax.set_ylabel("")
    ax.axvline(x=0, color="0.55", lw=0.8)
    if notebook_style:
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.axhline(y=up_df.shape[0] + down_df.shape[0], color="k", lw=0.5)
    else:
        ax.spines["bottom"].set_visible(True)
        ax.spines["bottom"].set_color("0.6")
        ax.spines["bottom"].set_linewidth(0.8)
        for spine_name in ["top", "left", "right"]:
            ax.spines[spine_name].set_visible(False)
        ax.tick_params(axis="x", length=3, color="0.6")
        ax.tick_params(axis="y", length=0)
        fig.tight_layout()
    return fig, plot_df


def _save_figure(
    fig: plt.Figure,
    png_path: Path,
    pdf_path: Path,
    tiff_path: Path,
    dpi: int,
    tiff_dpi: int,
) -> None:
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(tiff_path, dpi=tiff_dpi, bbox_inches="tight")
    plt.close(fig)


def _write_result_log(path: Path, values: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for key, value in values.items():
            handle.write(f"{key} = {value}\n")


def _write_n8n_script(path: Path, api_url: str, payload: dict[str, Any]) -> None:
    script = f"""const payload = {json.dumps(payload, indent=2)};

const response = await fetch("{api_url}", {{
  method: "POST",
  headers: {{ "Content-Type": "application/json" }},
  body: JSON.stringify(payload),
}});

if (!response.ok) {{
  throw new Error(`OmicsOne API failed: ${{response.status}} ${{await response.text()}}`);
}}

return await response.json();
"""
    path.write_text(script, encoding="utf-8")


def _write_html_report(
    path: Path,
    *,
    result_name: str,
    parameters: dict[str, Any],
    log_values: dict[str, Any],
    volcano_png: Path,
    enrichment_png: Path,
    up_enrichment: pd.DataFrame,
    down_enrichment: pd.DataFrame,
) -> None:
    def table_from_mapping(values: dict[str, Any]) -> str:
        rows = "\n".join(
            "<tr><th>{}</th><td>{}</td></tr>".format(escape(str(key)), escape(str(value)))
            for key, value in values.items()
        )
        return f"<table>{rows}</table>"

    def table_from_df(df: pd.DataFrame, max_rows: int = 10) -> str:
        if df.empty:
            return "<p>No enriched terms passed the configured FDR cutoff.</p>"
        cols = [
            col
            for col in ["Gene_set", "Term", "Overlap", "Adjusted P-value", "-Log10(Adj.P)", "Genes"]
            if col in df.columns
        ]
        return df.loc[:, cols].head(max_rows).to_html(index=False, escape=True, float_format="{:.4g}".format)

    report = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{escape(result_name)} volcano and enrichment report</title>
  <style>
    body {{ font-family: 'Liberation Sans', Arial, sans-serif; margin: 32px; color: #222; line-height: 1.45; }}
    h1, h2 {{ margin-top: 1.4em; }}
    table {{ border-collapse: collapse; width: 100%; margin: 12px 0 24px; font-size: 13px; }}
    th, td {{ border: 1px solid #ddd; padding: 6px 8px; vertical-align: top; }}
    th {{ text-align: left; background: #f2f2f2; width: 260px; }}
    img {{ max-width: 100%; border: 1px solid #ddd; margin: 8px 0 24px; }}
    .method {{ max-width: 980px; }}
  </style>
</head>
<body>
  <h1>{escape(result_name)} volcano and enrichment report</h1>
  <section class="method">
    <h2>Process and Method</h2>
    <p>Normal and tumor phosphosite abundance matrices were joined by feature ID. Differential abundance was tested with the configured two-group method after filtering features by global and per-group missingness and minimum non-missing sample counts. FDR values were calculated with Benjamini-Hochberg correction. Significant up and down features used the configured FDR and absolute log2 fold-change cutoffs.</p>
    <p>For enrichment, feature IDs were mapped to gene symbols through the configured FASTA file by taking the Ensembl gene prefix before the first pipe character. Up and down gene lists were de-duplicated, overlapping genes were removed from each direction, and Enrichr/GSEApy was run with the configured gene sets and background mode.</p>
  </section>
  <h2>Parameters</h2>
  {table_from_mapping(parameters)}
  <h2>Results Summary</h2>
  {table_from_mapping(log_values)}
  <h2>Volcano Plot</h2>
  <img src="{escape(volcano_png.name)}" alt="Volcano plot">
  <h2>Enrichment Plot</h2>
  <img src="{escape(enrichment_png.name)}" alt="Enrichment barchart">
  <h2>Top Up Enrichment Terms</h2>
  {table_from_df(up_enrichment)}
  <h2>Top Down Enrichment Terms</h2>
  {table_from_df(down_enrichment)}
</body>
</html>
"""
    path.write_text(report, encoding="utf-8")


def _split_config_list(value: str) -> list[str]:
    return [item.strip() for item in value.replace("\n", ",").split(",") if item.strip()]


def load_volcano_config(config_path: str | Path) -> dict[str, Any]:
    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Config file does not exist: {config_file}")

    parser = configparser.ConfigParser()
    parser.read(config_file, encoding="utf-8-sig")

    payload: dict[str, Any] = {}

    string_keys = {
        "analysis": ["cohort", "omics", "job_name", "method"],
        "input": ["normal_path", "tumor_path", "fasta_path"],
        "output": ["output_dir", "output_prefix", "api_url"],
        "enrichment": ["organism", "enrichment_background_mode"],
        "reference": ["up_enrichment_tsv", "down_enrichment_tsv"],
        "plot": [
            "title",
            "enrichment_title",
            "xlabel",
            "font_family",
            "background_color",
            "up_color",
            "down_color",
        ],
    }
    float_keys = {
        "diff": [
            "fdr_cutoff",
            "log2fc_cutoff",
            "max_miss_ratio_global",
            "max_miss_ratio_group",
        ],
        "enrichment": [
            "enrichment_fdr_cutoff",
            "enrichment_min_x",
            "enrichment_max_x",
        ],
        "plot": [
            "volcano_width",
            "volcano_height",
            "enrichment_width",
            "enrichment_height",
            "enrichment_width_ratio",
            "enrichment_height_ratio",
            "enrichment_size_scale",
            "font_size",
            "point_size",
            "significant_point_size",
        ],
    }
    int_keys = {
        "diff": ["min_sample_size"],
        "enrichment": ["enrichment_top_n"],
        "plot": ["dpi", "tiff_dpi"],
        "style": ["volcano_dpi", "enrichment_dpi"],
    }
    bool_keys = {
        "analysis": ["strip_feature_version"],
        "enrichment": ["prefer_local_gene_sets", "allow_remote_enrichr"],
        "output": ["write_n8n_script"],
        "plot": ["editable_pdf_text", "write_html_report"],
        "style": ["notebook_style_plots", "show_titles"],
    }
    list_keys = {"enrichment": ["gene_sets", "skip_pathways"]}

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
                payload[key] = _split_config_list(parser.get(section, key))

    return payload


def generate_volcano_enrichment_from_config(config_path: str | Path) -> VolcanoEnrichmentResult:
    return generate_volcano_enrichment(**load_volcano_config(config_path))


def generate_volcano_enrichment(
    *,
    normal_path: str,
    tumor_path: str,
    output_dir: str,
    fasta_path: str = DEFAULT_FASTA_PATH,
    cohort: str = "HNSCC",
    omics: str = "RNA",
    job_name: str | None = None,
    method: str = "Wilcoxon(Unpaired)",
    strip_feature_version: bool = True,
    fdr_cutoff: float = 0.01,
    log2fc_cutoff: float = 1.0,
    max_miss_ratio_global: float = 0.5,
    max_miss_ratio_group: float = 0.5,
    min_sample_size: int = 4,
    gene_sets: list[str] | None = None,
    enrichment_fdr_cutoff: float = 0.05,
    organism: str = "human",
    enrichment_background_mode: str = "gene_list",
    prefer_local_gene_sets: bool = True,
    allow_remote_enrichr: bool = False,
    up_enrichment_tsv: str | None = None,
    down_enrichment_tsv: str | None = None,
    skip_pathways: list[str] | None = None,
    enrichment_top_n: int = 10,
    title: str | None = None,
    enrichment_title: str | None = None,
    xlabel: str = "Log2FC(Tumor/NAT)",
    volcano_width: float = 4.0,
    volcano_height: float = 4.0,
    enrichment_width: float | None = None,
    enrichment_height: float | None = None,
    enrichment_width_ratio: float = DEFAULT_ENRICHMENT_WIDTH_RATIO,
    enrichment_height_ratio: float = DEFAULT_ENRICHMENT_HEIGHT_RATIO,
    enrichment_size_scale: float = DEFAULT_ENRICHMENT_SIZE_SCALE,
    enrichment_min_x: float = -25.0,
    enrichment_max_x: float = 40.0,
    dpi: int = 300,
    tiff_dpi: int = 600,
    volcano_dpi: int | None = None,
    enrichment_dpi: int | None = None,
    font_family: str = "Liberation Sans",
    font_size: float = 10.0,
    editable_pdf_text: bool = True,
    background_color: str = "#808080",
    up_color: str = "#FF0000",
    down_color: str = "#0000FF",
    point_size: float = 1.0,
    significant_point_size: float = 5.0,
    output_prefix: str | None = None,
    write_html_report: bool = True,
    write_n8n_script: bool = False,
    notebook_style_plots: bool = False,
    show_titles: bool = True,
    api_url: str = DEFAULT_API_URL,
) -> VolcanoEnrichmentResult:
    gene_sets = gene_sets or ["MSigDB_Hallmark_2020"]
    skip_pathways = DEFAULT_SKIP_PATHWAYS if skip_pathways is None else skip_pathways
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    normal = _read_group_matrix(Path(normal_path), ".N", strip_feature_version)
    tumor = _read_group_matrix(Path(tumor_path), ".T", strip_feature_version)
    data = tumor.join(normal)
    tumor_samples = [column for column in data.columns.values if column.endswith(".T")]
    normal_samples = [column for column in data.columns.values if column.endswith(".N")]

    diff = _compare_two_groups_fast(
        data,
        tumor_samples,
        normal_samples,
        method=method,
        fdr_cutoff=fdr_cutoff,
        log2fc_cutoff=log2fc_cutoff,
        max_miss_ratio_global=max_miss_ratio_global,
        max_miss_ratio_group=max_miss_ratio_group,
        min_sample_size=min_sample_size,
    )

    gene_map = _read_gene_map(fasta_path)
    total_gene_ids = [str(index) for index in data.index]
    up_gene_ids = [str(index) for index in diff[diff["Significance"] == "S-U"].index]
    down_gene_ids = [str(index) for index in diff[diff["Significance"] == "S-D"].index]
    total_genes = _map_gene_ids(total_gene_ids, gene_map)
    up_genes = set(_map_gene_ids(up_gene_ids, gene_map))
    down_genes = set(_map_gene_ids(down_gene_ids, gene_map))
    pure_up_genes = sorted(up_genes - down_genes)
    pure_down_genes = sorted(down_genes - up_genes)

    prefix = output_prefix or (job_name or f"{cohort}_{omics}").lower()
    diff_tsv = output_path / "diff.tsv"
    combined_matrix_tsv = output_path / "combined_matrix.tsv"
    up_genes_tsv = output_path / "pure_up_genes.tsv"
    down_genes_tsv = output_path / "pure_down_genes.tsv"
    total_genes_tsv = output_path / "total_genes.tsv"
    output_up_enrichment_tsv = output_path / "up_enrichr_df.tsv"
    output_down_enrichment_tsv = output_path / "down_enrichr_df.tsv"
    enrichment_plot_tsv = output_path / "enrichment_plot_table.tsv"
    volcano_png = output_path / f"{prefix}_volcano.png"
    volcano_pdf = output_path / f"{prefix}_volcano.pdf"
    volcano_tiff = output_path / f"{prefix}_volcano.tiff"
    enrichment_png = output_path / f"{prefix}_enrichment_barchart.png"
    enrichment_pdf = output_path / f"{prefix}_enrichment_barchart.pdf"
    enrichment_tiff = output_path / f"{prefix}_enrichment_barchart.tiff"
    report_html = output_path / f"{prefix}_report.html"
    result_log = output_path / "result.log"
    n8n_js = output_path / f"run_{prefix}_volcano_enrichment_n8n.js"

    data.to_csv(combined_matrix_tsv, sep="\t")
    diff.to_csv(diff_tsv, sep="\t")
    _save_gene_list(up_genes_tsv, pure_up_genes)
    _save_gene_list(down_genes_tsv, pure_down_genes)
    _save_gene_list(total_genes_tsv, total_genes)

    if up_enrichment_tsv and down_enrichment_tsv:
        up_enrichment = pd.read_csv(up_enrichment_tsv, sep="\t")
        down_enrichment = pd.read_csv(down_enrichment_tsv, sep="\t")
    else:
        up_enrichment = _enrichr(
            pure_up_genes,
            total_genes,
            output_path,
            gene_sets,
            organism,
            enrichment_fdr_cutoff,
            enrichment_background_mode,
            prefer_local_gene_sets,
            allow_remote_enrichr,
        )
        down_enrichment = _enrichr(
            pure_down_genes,
            total_genes,
            output_path,
            gene_sets,
            organism,
            enrichment_fdr_cutoff,
            enrichment_background_mode,
            prefer_local_gene_sets,
            allow_remote_enrichr,
        )
    up_enrichment.to_csv(output_up_enrichment_tsv, sep="\t", index=False)
    down_enrichment.to_csv(output_down_enrichment_tsv, sep="\t", index=False)

    filtered_up = (
        up_enrichment[~up_enrichment["Term"].isin(skip_pathways)]
        .sort_values("Adjusted P-value")
        .head(enrichment_top_n)
    )
    filtered_down = (
        down_enrichment[~down_enrichment["Term"].isin(skip_pathways)]
        .sort_values("Adjusted P-value")
        .head(enrichment_top_n)
    )
    effective_enrichment_width = (
        enrichment_width
        if enrichment_width is not None
        else enrichment_width_ratio * enrichment_size_scale
    )
    effective_enrichment_height = (
        enrichment_height
        if enrichment_height is not None
        else enrichment_height_ratio * enrichment_size_scale
    )

    colors = {"background": background_color, "up": up_color, "down": down_color}
    rc_updates: dict[str, Any] = {
        "font.family": "sans-serif",
        "font.sans-serif": [font_family],
        "font.size": font_size,
        "axes.titleweight": "regular",
        "axes.labelweight": "regular",
    }
    if editable_pdf_text:
        rc_updates.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})

    effective_volcano_dpi = volcano_dpi or dpi
    effective_enrichment_dpi = enrichment_dpi or dpi
    with matplotlib.rc_context(rc_updates):
        volcano_fig = _plot_volcano(
            diff,
            title or f"{cohort} {omics} differential expression analysis",
            xlabel,
            volcano_width,
            volcano_height,
            effective_volcano_dpi,
            colors,
            point_size,
            significant_point_size,
            fdr_cutoff,
            log2fc_cutoff,
            notebook_style_plots,
            show_titles,
        )
        _save_figure(
            volcano_fig,
            volcano_png,
            volcano_pdf,
            volcano_tiff,
            effective_volcano_dpi,
            tiff_dpi,
        )

        enrichment_fig, enrichment_plot_df = _plot_enrichment_bars(
            filtered_up,
            filtered_down,
            title=enrichment_title or f"{cohort} {omics} Pathways Enrichment Analysis",
            min_x=enrichment_min_x,
            max_x=enrichment_max_x,
            width=effective_enrichment_width,
            height=effective_enrichment_height,
            dpi=effective_enrichment_dpi,
            colors=colors,
            notebook_style=notebook_style_plots,
            show_title=show_titles,
        )
        _save_figure(
            enrichment_fig,
            enrichment_png,
            enrichment_pdf,
            enrichment_tiff,
            effective_enrichment_dpi,
            tiff_dpi,
        )
    enrichment_plot_df.to_csv(enrichment_plot_tsv, sep="\t", index=False)

    payload = {
        "normal_path": normal_path,
        "tumor_path": tumor_path,
        "output_dir": output_dir,
        "fasta_path": fasta_path,
        "cohort": cohort,
        "omics": omics,
        "job_name": job_name,
        "method": method,
        "strip_feature_version": strip_feature_version,
        "fdr_cutoff": fdr_cutoff,
        "log2fc_cutoff": log2fc_cutoff,
        "max_miss_ratio_global": max_miss_ratio_global,
        "max_miss_ratio_group": max_miss_ratio_group,
        "min_sample_size": min_sample_size,
        "gene_sets": gene_sets,
        "enrichment_fdr_cutoff": enrichment_fdr_cutoff,
        "organism": organism,
        "enrichment_background_mode": enrichment_background_mode,
        "prefer_local_gene_sets": prefer_local_gene_sets,
        "allow_remote_enrichr": allow_remote_enrichr,
        "up_enrichment_tsv": up_enrichment_tsv,
        "down_enrichment_tsv": down_enrichment_tsv,
        "skip_pathways": skip_pathways,
        "enrichment_top_n": enrichment_top_n,
        "title": title,
        "enrichment_title": enrichment_title,
        "xlabel": xlabel,
        "volcano_width": volcano_width,
        "volcano_height": volcano_height,
        "enrichment_width": enrichment_width,
        "enrichment_height": enrichment_height,
        "enrichment_width_ratio": enrichment_width_ratio,
        "enrichment_height_ratio": enrichment_height_ratio,
        "enrichment_size_scale": enrichment_size_scale,
        "enrichment_min_x": enrichment_min_x,
        "enrichment_max_x": enrichment_max_x,
        "dpi": dpi,
        "tiff_dpi": tiff_dpi,
        "volcano_dpi": volcano_dpi,
        "enrichment_dpi": enrichment_dpi,
        "font_family": font_family,
        "font_size": font_size,
        "editable_pdf_text": editable_pdf_text,
        "background_color": background_color,
        "up_color": up_color,
        "down_color": down_color,
        "point_size": point_size,
        "significant_point_size": significant_point_size,
        "output_prefix": output_prefix,
        "write_html_report": write_html_report,
        "write_n8n_script": write_n8n_script,
        "notebook_style_plots": notebook_style_plots,
        "show_titles": show_titles,
    }
    if write_n8n_script:
        _write_n8n_script(n8n_js, api_url, payload)

    log_values = {
        "cohort": cohort,
        "omics": omics,
        "normal_path": normal_path,
        "tumor_path": tumor_path,
        "fasta_path": fasta_path,
        "output_dir": output_dir,
        "method": method,
        "strip_feature_version": strip_feature_version,
        "fdr_cutoff": fdr_cutoff,
        "log2fc_cutoff": log2fc_cutoff,
        "gene_sets": ",".join(gene_sets),
        "enrichment_background_mode": enrichment_background_mode,
        "prefer_local_gene_sets": prefer_local_gene_sets,
        "allow_remote_enrichr": allow_remote_enrichr,
        "reference_up_enrichment_tsv": up_enrichment_tsv,
        "reference_down_enrichment_tsv": down_enrichment_tsv,
        "feature_count": data.shape[0],
        "sample_count": data.shape[1],
        "normal_sample_count": len(normal_samples),
        "tumor_sample_count": len(tumor_samples),
        "diff_feature_count": diff.shape[0],
        "up_count": len(up_gene_ids),
        "down_count": len(down_gene_ids),
        "pure_up_gene_count": len(pure_up_genes),
        "pure_down_gene_count": len(pure_down_genes),
        "total_gene_count": len(total_genes),
        "up_enrichment_count": up_enrichment.shape[0],
        "down_enrichment_count": down_enrichment.shape[0],
        "enrichment_width_ratio": enrichment_width_ratio,
        "enrichment_height_ratio": enrichment_height_ratio,
        "enrichment_size_scale": enrichment_size_scale,
        "dpi": dpi,
        "tiff_dpi": tiff_dpi,
        "effective_volcano_dpi": effective_volcano_dpi,
        "effective_enrichment_dpi": effective_enrichment_dpi,
        "effective_enrichment_width": effective_enrichment_width,
        "effective_enrichment_height": effective_enrichment_height,
        "volcano_png": volcano_png,
        "volcano_pdf": volcano_pdf,
        "volcano_tiff": volcano_tiff,
        "enrichment_png": enrichment_png,
        "enrichment_pdf": enrichment_pdf,
        "enrichment_tiff": enrichment_tiff,
        "report_html": report_html,
        "diff_tsv": diff_tsv,
        "up_enrichment_tsv": output_up_enrichment_tsv,
        "down_enrichment_tsv": output_down_enrichment_tsv,
        "n8n_js": n8n_js,
    }
    _write_result_log(result_log, log_values)
    if write_html_report:
        _write_html_report(
            report_html,
            result_name=prefix,
            parameters=payload,
            log_values=log_values,
            volcano_png=volcano_png,
            enrichment_png=enrichment_png,
            up_enrichment=up_enrichment,
            down_enrichment=down_enrichment,
        )

    return VolcanoEnrichmentResult(
        cohort=cohort,
        omics=omics,
        output_dir=output_path,
        diff_tsv=diff_tsv,
        combined_matrix_tsv=combined_matrix_tsv,
        up_genes_tsv=up_genes_tsv,
        down_genes_tsv=down_genes_tsv,
        total_genes_tsv=total_genes_tsv,
        up_enrichment_tsv=output_up_enrichment_tsv,
        down_enrichment_tsv=output_down_enrichment_tsv,
        enrichment_plot_tsv=enrichment_plot_tsv,
        volcano_png=volcano_png,
        volcano_pdf=volcano_pdf,
        volcano_tiff=volcano_tiff,
        enrichment_png=enrichment_png,
        enrichment_pdf=enrichment_pdf,
        enrichment_tiff=enrichment_tiff,
        report_html=report_html,
        result_log=result_log,
        n8n_js=n8n_js,
        feature_count=data.shape[0],
        diff_feature_count=diff.shape[0],
        up_count=len(up_gene_ids),
        down_count=len(down_gene_ids),
        pure_up_gene_count=len(pure_up_genes),
        pure_down_gene_count=len(pure_down_genes),
        total_gene_count=len(total_genes),
        up_enrichment_count=up_enrichment.shape[0],
        down_enrichment_count=down_enrichment.shape[0],
        method=method,
        fdr_cutoff=fdr_cutoff,
        log2fc_cutoff=log2fc_cutoff,
        gene_sets=gene_sets,
    )
