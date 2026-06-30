from typing import Optional

from pydantic import BaseModel, Field

from omicsone.services.boxplot_figures import (
    DEFAULT_EXAMPLE_GENES,
    DEFAULT_FASTA_PATH,
)


class BoxplotFiguresRequest(BaseModel):
    normal_path: str = Field(description="Path to the normal/NAT matrix TSV.")
    tumor_path: str = Field(description="Path to the tumor matrix TSV.")
    output_dir: str = Field(description="Directory where generated PDF boxplots are written.")
    fasta_path: str = Field(
        default=DEFAULT_FASTA_PATH,
        description="FASTA path used to map Ensembl gene IDs to gene symbols.",
    )
    cohort: str = Field(default="LSCC", description="Cohort label, for example HNSCC or LSCC.")
    omics: str = Field(default="Protein", description="Omics label, for example RNA or Protein.")
    genes: list[str] = Field(
        default_factory=lambda: list(DEFAULT_EXAMPLE_GENES),
        description="Gene symbols or Ensembl gene IDs to plot. One PDF is generated per gene.",
    )
    output_prefix: Optional[str] = Field(
        default=None,
        description="Optional prefix used in generated PDF filenames.",
    )
    ylabel: str = Field(default="Log2 abundance")
    width: float = Field(default=2.25, gt=0, description="Single boxplot width in inches.")
    height: float = Field(default=3.0, gt=0, description="Single boxplot height in inches.")
    dpi: int = Field(default=300, ge=72)
    font_family: str = Field(default="Liberation Sans")
    font_size: float = Field(default=9.0, gt=0)
    editable_pdf_text: bool = Field(
        default=True,
        description="When true, PDF text is kept editable by using TrueType font embedding.",
    )
    tumor_color: str = Field(default="#1f77b4")
    normal_color: str = Field(default="#ff7f0e")
    write_n8n_script: bool = Field(default=False, description="Write a JavaScript helper for n8n/API automation.")


class BoxplotFigureRecordResponse(BaseModel):
    gene: str
    gene_id: str
    pdf: str
    pvalue: float
    significance: str
    normal_count: int
    tumor_count: int


class BoxplotFiguresResponse(BaseModel):
    cohort: str
    omics: str
    output_dir: str
    boxplot_pdfs: list[str]
    summary_tsv: str
    result_log: str
    n8n_js: str
    generated_count: int
    missing_genes: list[str]
    records: list[BoxplotFigureRecordResponse]
