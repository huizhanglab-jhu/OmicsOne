from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from pyteomics import fasta
from scipy import stats


DEFAULT_FASTA_PATH = (
    r"F:\lab\HsinI\Head and Neck & Lung\fasta"
    r"\GENCODE.V42.basic.CHR.combined_contaminants.gpquest3.fasta"
)
DEFAULT_API_URL = "http://127.0.0.1:8001/api/v1/diff/boxplot/figures"
DEFAULT_EXAMPLE_GENES = ["MYL1", "MYL2", "MYL3", "TNNC1", "TNNC2", "CFD"]


@dataclass(frozen=True)
class BoxplotFigureRecord:
    gene: str
    gene_id: str
    pdf: Path
    pvalue: float
    significance: str
    normal_count: int
    tumor_count: int


@dataclass(frozen=True)
class BoxplotFiguresResult:
    cohort: str
    omics: str
    output_dir: Path
    boxplot_pdfs: list[Path]
    summary_tsv: Path
    result_log: Path
    n8n_js: Path
    generated_count: int
    missing_genes: list[str]
    records: list[BoxplotFigureRecord]


def _ensure_file(path: str | Path, label: str) -> Path:
    resolved = Path(path).expanduser()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} does not exist: {resolved}")
    return resolved


def _normalize_gene_id(value: object) -> str:
    return str(value).strip().split(".", maxsplit=1)[0]


def _safe_name(value: str) -> str:
    allowed = []
    for character in value.strip():
        if character.isalnum() or character in {"-", "_"}:
            allowed.append(character)
        else:
            allowed.append("_")
    return "".join(allowed).strip("_") or "boxplot"


def _read_gene_symbol_map(fasta_path: Path) -> dict[str, str]:
    gene_map: dict[str, str] = {}
    for description, _sequence in fasta.read(str(fasta_path)):
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


def _read_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", na_values=["NA", "NaN", ""])
    if df.empty:
        raise ValueError(f"Input matrix is empty: {path}")

    if "idx" in df.columns:
        df = df.set_index("idx")
    else:
        df = df.set_index(df.columns[0])

    df.index = [_normalize_gene_id(index) for index in df.index]
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.index.has_duplicates:
        df = df.groupby(level=0).median(numeric_only=True)
    return df


def _resolve_gene_ids(
    genes: list[str],
    available_gene_ids: set[str],
    gene_symbol_map: dict[str, str],
) -> tuple[list[tuple[str, str]], list[str]]:
    symbol_to_ids: dict[str, list[str]] = {}
    for gene_id, symbol in gene_symbol_map.items():
        if gene_id in available_gene_ids:
            symbol_to_ids.setdefault(symbol.upper(), []).append(gene_id)

    resolved = []
    missing = []
    for gene in genes:
        query = str(gene).strip()
        if not query:
            continue

        normalized_query = _normalize_gene_id(query)
        if normalized_query in available_gene_ids:
            symbol = gene_symbol_map.get(normalized_query, query)
            resolved.append((symbol, normalized_query))
            continue

        matching_ids = sorted(symbol_to_ids.get(query.upper(), []))
        if matching_ids:
            resolved.append((query, matching_ids[0]))
        else:
            missing.append(query)

    return resolved, missing


def _significance_symbol(pvalue: float) -> str:
    if pvalue < 0.001:
        return "***"
    if pvalue < 0.01:
        return "**"
    if pvalue < 0.05:
        return "*"
    return "ns"


def _plot_one_boxplot(
    data: pd.DataFrame,
    gene: str,
    ylabel: str,
    output_path: Path,
    width: float,
    height: float,
    dpi: int,
    palette: dict[str, str],
) -> None:
    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    sns.boxplot(
        data=data,
        x="Group",
        y="Value",
        hue="Group",
        order=["Tumor", "Normal"],
        hue_order=["Tumor", "Normal"],
        palette=palette,
        ax=ax,
        width=0.6,
        fliersize=3,
        linewidth=1,
    )
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()

    tumor_values = data.loc[data["Group"] == "Tumor", "Value"].dropna()
    normal_values = data.loc[data["Group"] == "Normal", "Value"].dropna()
    _stat, pvalue = stats.mannwhitneyu(
        tumor_values,
        normal_values,
        alternative="two-sided",
    )
    sig_symbol = _significance_symbol(float(pvalue))

    y_max = float(data["Value"].max())
    y_min = float(data["Value"].min())
    y_range = y_max - y_min
    if y_range == 0:
        y_range = max(abs(y_max) * 0.1, 1.0)

    bar_y = y_max + y_range * 0.10
    text_y = y_max + y_range * 0.13
    ax.plot([0, 0, 1, 1], [bar_y, bar_y + y_range * 0.05, bar_y + y_range * 0.05, bar_y], color="0.3", linewidth=1)
    ax.text(0.5, text_y, sig_symbol, ha="center", va="bottom", fontsize=8)
    ax.set_ylim(y_min - y_range * 0.08, y_max + y_range * 0.25)
    ax.set_title(gene)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.tick_params(axis="both", labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("0.6")
        spine.set_linewidth(1)
    fig.tight_layout()
    fig.savefig(output_path, format="pdf")
    plt.close(fig)


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


def generate_boxplot_figures(
    *,
    normal_path: str,
    tumor_path: str,
    output_dir: str,
    fasta_path: str = DEFAULT_FASTA_PATH,
    cohort: str = "LSCC",
    omics: str = "Protein",
    genes: list[str] | None = None,
    output_prefix: str | None = None,
    ylabel: str = "Log2 abundance",
    width: float = 2.25,
    height: float = 3.0,
    dpi: int = 300,
    font_family: str = "Liberation Sans",
    font_size: float = 9.0,
    editable_pdf_text: bool = True,
    tumor_color: str = "#1f77b4",
    normal_color: str = "#ff7f0e",
    write_n8n_script: bool = False,
    api_url: str = DEFAULT_API_URL,
) -> BoxplotFiguresResult:
    normal_file = _ensure_file(normal_path, "Normal matrix")
    tumor_file = _ensure_file(tumor_path, "Tumor matrix")
    fasta_file = _ensure_file(fasta_path, "FASTA file")
    output_path = Path(output_dir).expanduser()
    output_path.mkdir(parents=True, exist_ok=True)

    target_genes = genes or list(DEFAULT_EXAMPLE_GENES)
    normal = _read_matrix(normal_file)
    tumor = _read_matrix(tumor_file)
    available_gene_ids = set(normal.index) & set(tumor.index)
    if not available_gene_ids:
        raise ValueError("No overlapping gene IDs found between normal and tumor matrices.")

    gene_symbol_map = _read_gene_symbol_map(fasta_file)
    resolved_genes, missing_genes = _resolve_gene_ids(
        target_genes,
        available_gene_ids,
        gene_symbol_map,
    )
    if not resolved_genes:
        raise ValueError(
            "None of the requested genes were found in both matrices after FASTA mapping."
        )

    prefix = output_prefix or f"{cohort}_{omics}".lower()
    summary_tsv = output_path / "boxplot_summary.tsv"
    result_log = output_path / "result.log"
    n8n_js = output_path / f"run_{_safe_name(prefix)}_boxplots_n8n.js"
    records: list[BoxplotFigureRecord] = []

    rc_updates: dict[str, Any] = {
        "font.family": "sans-serif",
        "font.sans-serif": [font_family],
        "font.size": font_size,
        "axes.titleweight": "regular",
        "axes.labelweight": "regular",
    }
    if editable_pdf_text:
        rc_updates.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})

    with matplotlib.rc_context(rc_updates):
        for gene_symbol, gene_id in resolved_genes:
            tumor_values = tumor.loc[gene_id].dropna()
            normal_values = normal.loc[gene_id].dropna()
            if tumor_values.empty or normal_values.empty:
                missing_genes.append(gene_symbol)
                continue

            data = pd.DataFrame(
                {
                    "Sample": list(tumor_values.index) + list(normal_values.index),
                    "Value": list(tumor_values.values) + list(normal_values.values),
                    "Group": ["Tumor"] * len(tumor_values) + ["Normal"] * len(normal_values),
                }
            )
            _stat, pvalue = stats.mannwhitneyu(
                tumor_values,
                normal_values,
                alternative="two-sided",
            )
            pdf = output_path / f"{_safe_name(prefix)}_{_safe_name(gene_symbol)}_boxplot.pdf"
            _plot_one_boxplot(
                data,
                gene_symbol,
                ylabel,
                pdf,
                width,
                height,
                dpi,
                {"Tumor": tumor_color, "Normal": normal_color},
            )
            records.append(
                BoxplotFigureRecord(
                    gene=gene_symbol,
                    gene_id=gene_id,
                    pdf=pdf,
                    pvalue=float(pvalue),
                    significance=_significance_symbol(float(pvalue)),
                    normal_count=int(normal_values.shape[0]),
                    tumor_count=int(tumor_values.shape[0]),
                )
            )

    summary = pd.DataFrame(
        [
            {
                "gene": record.gene,
                "gene_id": record.gene_id,
                "pdf": str(record.pdf),
                "pvalue": record.pvalue,
                "significance": record.significance,
                "normal_count": record.normal_count,
                "tumor_count": record.tumor_count,
            }
            for record in records
        ]
    )
    summary.to_csv(summary_tsv, sep="\t", index=False)

    payload = {
        "normal_path": normal_path,
        "tumor_path": tumor_path,
        "output_dir": output_dir,
        "fasta_path": fasta_path,
        "cohort": cohort,
        "omics": omics,
        "genes": target_genes,
        "output_prefix": output_prefix,
        "ylabel": ylabel,
        "width": width,
        "height": height,
        "dpi": dpi,
        "font_family": font_family,
        "font_size": font_size,
        "editable_pdf_text": editable_pdf_text,
        "tumor_color": tumor_color,
        "normal_color": normal_color,
        "write_n8n_script": write_n8n_script,
    }
    if write_n8n_script:
        _write_n8n_script(n8n_js, api_url, payload)

    log_values = {
        "cohort": cohort,
        "omics": omics,
        "normal_path": normal_file,
        "tumor_path": tumor_file,
        "fasta_path": fasta_file,
        "output_dir": output_path,
        "requested_genes": ",".join(target_genes),
        "generated_count": len(records),
        "missing_genes": ",".join(missing_genes),
        "summary_tsv": summary_tsv,
        "n8n_js": n8n_js,
    }
    with result_log.open("w", encoding="utf-8") as handle:
        for key, value in log_values.items():
            handle.write(f"{key} = {value}\n")

    return BoxplotFiguresResult(
        cohort=cohort,
        omics=omics,
        output_dir=output_path,
        boxplot_pdfs=[record.pdf for record in records],
        summary_tsv=summary_tsv,
        result_log=result_log,
        n8n_js=n8n_js,
        generated_count=len(records),
        missing_genes=missing_genes,
        records=records,
    )
