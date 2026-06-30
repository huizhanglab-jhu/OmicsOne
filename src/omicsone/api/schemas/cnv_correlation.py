from typing import Literal

from pydantic import BaseModel, Field


TargetType = Literal["protein", "rna"]


class CnvCorrelationFiguresRequest(BaseModel):
    correlation_file: str = Field(
        description="Spearman output file with row_index_1, row_index_2, correlation columns.",
    )
    cnv_for_corr_file: str = Field(
        description="CNV *_for_corr.txt file with idx column and sample headers.",
    )
    target_for_corr_file: str = Field(
        description="Protein or RNA *_for_corr.txt file with idx column and sample headers.",
    )
    fasta_file: str = Field(
        description="GENCODE FASTA file used to map Ensembl genes to genome locations.",
    )
    chromosomes_file: str = Field(
        description="Chromosome length XLSX file.",
    )
    cytoband_file: str = Field(
        description="UCSC-style cytoband TSV file.",
    )
    output_dir: str = Field(
        description="Directory where figures and intermediate annotation tables are written.",
    )
    gistic_file: str | None = Field(
        default=None,
        description="Optional GISTIC-level CNV file. When provided, GISTIC panel files are generated.",
    )
    target_type: TargetType = Field(
        default="protein",
        description="Target omics type paired with CNV.",
    )
    correlation_threshold: float = Field(
        default=0.5,
        ge=0,
        le=1,
        description="Absolute Spearman correlation threshold used for plotting.",
    )
    output_prefix: str | None = Field(
        default=None,
        description="Optional prefix used in output filenames.",
    )
    chunksize: int = Field(
        default=1_000_000,
        ge=1,
        description="Rows per chunk when reading the Spearman correlation file.",
    )
    dpi: int = Field(
        default=150,
        ge=72,
        description="DPI for PNG rendering. TIFF output is always saved at 600 DPI.",
    )


class CnvCorrelationFigureSetResponse(BaseModel):
    output_dir: str
    annotated_correlations_file: str
    cnv_distribution_counts_file: str
    gistic_counts_file: str | None
    png_files: dict[str, str | None]
    pdf_files: dict[str, str | None]
    tiff_files: dict[str, str | None]
    filtered_correlation_count: int
    annotated_correlation_count: int
