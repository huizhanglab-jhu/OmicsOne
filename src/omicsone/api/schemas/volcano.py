from typing import Optional

from pydantic import BaseModel, Field

from omicsone.services.volcano_enrichment import (
    ENRICHMENT_BACKGROUND_MODES,
    DEFAULT_API_URL,
    DEFAULT_ENRICHMENT_HEIGHT_RATIO,
    DEFAULT_ENRICHMENT_SIZE_SCALE,
    DEFAULT_ENRICHMENT_WIDTH_RATIO,
    DEFAULT_FASTA_PATH,
    DEFAULT_PRESETS,
    DEFAULT_SKIP_PATHWAYS,
)


class VolcanoEnrichmentRequest(BaseModel):
    normal_path: str = Field(description="Path to the normal/NAT expression matrix TSV.")
    tumor_path: str = Field(description="Path to the tumor expression matrix TSV.")
    output_dir: str = Field(description="Directory for tables, plots, result.log, and optional n8n JS.")
    fasta_path: str = Field(
        default=DEFAULT_FASTA_PATH,
        description="FASTA path used to map Ensembl gene IDs to gene symbols.",
    )
    cohort: str = Field(default="HNSCC", description="Cohort label, for example HNSCC or LSCC.")
    omics: str = Field(default="RNA", description="Omics label, for example RNA or Protein.")
    job_name: Optional[str] = Field(default=None, description="Optional job name used in output naming.")
    method: str = Field(default="Wilcoxon(Unpaired)", description="Differential test method.")
    strip_feature_version: bool = Field(
        default=True,
        description=(
            "Strip text after the first dot in feature IDs. Set false for phosphosite IDs "
            "such as ENSG...|ENSP...|S123 so each site remains unique."
        ),
    )
    fdr_cutoff: float = Field(default=0.01, gt=0, lt=1, description="FDR cutoff for volcano significance.")
    log2fc_cutoff: float = Field(default=1.0, ge=0, description="Absolute log2 fold-change cutoff.")
    max_miss_ratio_global: float = Field(default=0.5, ge=0, le=1)
    max_miss_ratio_group: float = Field(default=0.5, ge=0, le=1)
    min_sample_size: int = Field(default=4, ge=1)
    gene_sets: list[str] = Field(default_factory=lambda: ["MSigDB_Hallmark_2020"])
    enrichment_fdr_cutoff: float = Field(default=0.05, gt=0, lt=1)
    organism: str = Field(default="human")
    enrichment_background_mode: str = Field(
        default="gene_list",
        description=(
            "Use 'gene_list' for local GMT enrichment with the actual background "
            "gene list, or 'count'/'notebook' to match notebooks that call "
            "enrichr(pure_*_genes, total_genes, job_dir, gene_sets=[...]), "
            "which uses background=len(total_genes). Local cached GMT files are used "
            "when available; unresolved Enrichr library names require allow_remote_enrichr=true. "
            f"Allowed values: {', '.join(sorted(ENRICHMENT_BACKGROUND_MODES))}. "
        ),
    )
    prefer_local_gene_sets: bool = Field(
        default=True,
        description="Use local GMT files from the gseapy cache or explicit GMT paths before Enrichr.",
    )
    allow_remote_enrichr: bool = Field(
        default=False,
        description="Permit remote Enrichr submission when a requested gene set is not available locally.",
    )
    up_enrichment_tsv: Optional[str] = Field(
        default=None,
        description="Optional TSV of precomputed up enrichment results to use for plotting/reporting.",
    )
    down_enrichment_tsv: Optional[str] = Field(
        default=None,
        description="Optional TSV of precomputed down enrichment results to use for plotting/reporting.",
    )
    skip_pathways: list[str] = Field(default_factory=lambda: list(DEFAULT_SKIP_PATHWAYS))
    enrichment_top_n: int = Field(default=10, ge=1)
    title: Optional[str] = Field(default=None, description="Volcano plot title.")
    enrichment_title: Optional[str] = Field(default=None, description="Enrichment barchart title.")
    xlabel: str = Field(default="Log2FC(Tumor/NAT)")
    volcano_width: float = Field(default=4.0, gt=0)
    volcano_height: float = Field(default=4.0, gt=0)
    enrichment_width: Optional[float] = Field(
        default=None,
        gt=0,
        description="Explicit enrichment barchart width in inches. Overrides ratio sizing when set.",
    )
    enrichment_height: Optional[float] = Field(
        default=None,
        gt=0,
        description="Explicit enrichment barchart height in inches. Overrides ratio sizing when set.",
    )
    enrichment_width_ratio: float = Field(
        default=DEFAULT_ENRICHMENT_WIDTH_RATIO,
        gt=0,
        description="Width component of enrichment barchart aspect ratio. Default ratio is 2:1.5.",
    )
    enrichment_height_ratio: float = Field(
        default=DEFAULT_ENRICHMENT_HEIGHT_RATIO,
        gt=0,
        description="Height component of enrichment barchart aspect ratio. Default ratio is 2:1.5.",
    )
    enrichment_size_scale: float = Field(
        default=DEFAULT_ENRICHMENT_SIZE_SCALE,
        gt=0,
        description="Scale multiplier applied to enrichment_width_ratio and enrichment_height_ratio.",
    )
    enrichment_min_x: float = Field(default=-25.0)
    enrichment_max_x: float = Field(default=40.0)
    dpi: int = Field(default=300, ge=72)
    tiff_dpi: int = Field(default=600, ge=72)
    volcano_dpi: Optional[int] = Field(default=None, ge=72)
    enrichment_dpi: Optional[int] = Field(default=None, ge=72)
    font_family: str = Field(default="Liberation Sans")
    font_size: float = Field(default=10.0, gt=0)
    editable_pdf_text: bool = Field(
        default=True,
        description="When true, PDF text is kept editable by using TrueType font embedding.",
    )
    background_color: str = Field(default="#808080")
    up_color: str = Field(default="#FF0000")
    down_color: str = Field(default="#0000FF")
    point_size: float = Field(default=1.0, gt=0)
    significant_point_size: float = Field(default=5.0, gt=0)
    output_prefix: Optional[str] = None
    write_html_report: bool = Field(default=True)
    notebook_style_plots: bool = Field(
        default=False,
        description="Use legacy omicsone_core plotting style for notebook-matched figures.",
    )
    show_titles: bool = Field(default=True, description="Show plot titles on generated figures.")
    write_n8n_script: bool = Field(default=False, description="Write a JavaScript helper for n8n/API automation.")
    api_url: str = Field(default=DEFAULT_API_URL, description="URL written into generated n8n JS when enabled.")


class VolcanoPresetRequest(BaseModel):
    preset: str = Field(
        default="HNSCC_RNA",
        description=f"One of: {', '.join(DEFAULT_PRESETS.keys())}",
    )
    overrides: dict = Field(
        default_factory=dict,
        description="Optional request field overrides applied after resolving the preset.",
    )


class VolcanoEnrichmentResponse(BaseModel):
    cohort: str
    omics: str
    output_dir: str
    diff_tsv: str
    combined_matrix_tsv: str
    up_genes_tsv: str
    down_genes_tsv: str
    total_genes_tsv: str
    up_enrichment_tsv: str
    down_enrichment_tsv: str
    enrichment_plot_tsv: str
    volcano_png: str
    volcano_pdf: str
    volcano_tiff: str
    enrichment_png: str
    enrichment_pdf: str
    enrichment_tiff: str
    report_html: str
    result_log: str
    n8n_js: str
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
