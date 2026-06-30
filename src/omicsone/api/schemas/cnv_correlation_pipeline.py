from typing import Any

from pydantic import BaseModel, Field, model_validator


class CnvCorrelationPipelineRequest(BaseModel):
    cohort: str = Field(
        default="omicsone",
        description="Cohort label used in output filenames, for example HNSCC or LSCC.",
    )
    cnv_path: str = Field(description="Raw CNV matrix with an idx column.")
    rna_path: str | None = Field(
        default=None,
        description="Optional raw RNA matrix with an idx column.",
    )
    protein_path: str | None = Field(
        default=None,
        description="Optional raw protein matrix with an idx column.",
    )
    gistic_path: str | None = Field(
        default=None,
        description="Optional GISTIC-level CNV matrix used for the GISTIC panel.",
    )
    fasta_file: str = Field(
        description="GENCODE FASTA file used for gene location annotation.",
    )
    chromosomes_file: str = Field(description="Chromosome length XLSX file.")
    cytoband_file: str = Field(description="UCSC-style cytoband TSV file.")
    output_dir: str = Field(
        description="Output folder for pairwise results, figures, and reports.",
    )
    min_valid_pairs: int = Field(default=2, ge=1)
    correlation_threshold: float = Field(default=0.5, ge=0, le=1)
    chunksize: int = Field(default=1_000_000, ge=1)
    dpi: int = Field(
        default=150,
        ge=72,
        description="PNG DPI. TIFF files are always saved at 600 DPI.",
    )
    use_three_way_common_genes: bool = Field(
        default=False,
        description=(
            "When CNV, RNA, and protein are provided, restrict each pairwise "
            "correlation to the shared CNV/RNA/protein gene and sample universe."
        ),
    )
    generate_clean_figures: bool = Field(
        default=True,
        description=(
            "Generate standardized clean CNV correlation figures with simple "
            "filenames after each pairwise figure set is created."
        ),
    )

    @model_validator(mode="after")
    def validate_targets(self):
        if self.rna_path is None and self.protein_path is None:
            raise ValueError("At least one of rna_path or protein_path must be provided")
        if self.use_three_way_common_genes and (
            self.rna_path is None or self.protein_path is None
        ):
            raise ValueError(
                "use_three_way_common_genes requires rna_path and protein_path"
            )
        return self


class CnvCorrelationPipelineStartResponse(BaseModel):
    job_id: str
    status: str
    status_url: str
    output_dir: str
    current_step: str


class CnvCorrelationPipelineStatusResponse(BaseModel):
    job_id: str
    status: str
    created_at: str
    started_at: str | None
    ended_at: str | None
    current_step: str
    output_dir: str
    html_report: str | None
    json_report: str | None
    error: str | None
    jobs: list[dict[str, Any]]
