from typing import Optional

from pydantic import BaseModel, Field

from omicsone.services.mutation_figures import (
    DEFAULT_COHORT,
    DEFAULT_MAF_PATH,
    DEFAULT_MUTATION_BINARY_PATH,
    DEFAULT_MUTATION_EXCEL_PATH,
    DEFAULT_MUTATION_META_PATH,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_V2_COHORT,
    DEFAULT_V2_MAF_PATH,
    DEFAULT_V2_OUTPUT_DIR,
)


class HnscMutationFiguresRequest(BaseModel):
    mutation_excel_path: str = Field(
        default=DEFAULT_MUTATION_EXCEL_PATH,
        description="Path to HNSCC_somatic_mutation_gene_level_revised.xlsx.",
    )
    maf_path: str = Field(
        default=DEFAULT_MAF_PATH,
        description="Path to the HNSCC somatic mutation MAF file.",
    )
    output_dir: str = Field(
        default=DEFAULT_OUTPUT_DIR,
        description="Directory where generated PDF figures are written.",
    )
    mutation_threshold: float = Field(
        default=0.15,
        ge=0,
        le=1,
        description="Minimum NUM_MUT threshold used to keep genes for plotting.",
    )
    cohort: str = Field(
        default=DEFAULT_COHORT,
        description="Cohort name used for default output filenames, for example hnscc or lscc.",
    )
    heatmap_filename: Optional[str] = Field(
        default=None,
        description="Optional explicit heatmap PDF filename. When omitted, derived from cohort.",
    )
    mutation_type_filename: Optional[str] = Field(
        default=None,
        description="Optional explicit mutation-type distribution PDF filename. When omitted, derived from cohort.",
    )
    output_prefix: Optional[str] = Field(
        default=None,
        description="Optional prefix added to generated PDF filenames.",
    )


class HnscMutationFiguresResponse(BaseModel):
    heatmap_pdf: str
    mutation_type_pdf: str
    result_log: str
    gene_summary_tsv: str
    binary_matrix_tsv: str
    mutation_type_matrix_tsv: str
    encoded_matrix_tsv: str
    sample_annotations_tsv: str
    sample_annotation_colors_tsv: str
    mutation_color_table_tsv: str
    filtered_maf_tsv: str
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


class MutationHeatmapV2Request(BaseModel):
    mutation_binary_path: str = Field(
        default=DEFAULT_MUTATION_BINARY_PATH,
        description="Path to HNSCC_somatic_mutation_gene_level_binary.txt.",
    )
    meta_path: str = Field(
        default=DEFAULT_MUTATION_META_PATH,
        description="Path to HNSCC_meta.txt with Stage and Histologic_Grade columns.",
    )
    maf_path: str = Field(
        default=DEFAULT_V2_MAF_PATH,
        description="Path to the somatic mutation MAF file.",
    )
    output_dir: str = Field(
        default=DEFAULT_V2_OUTPUT_DIR,
        description="Directory where generated V2 PDF figures are written.",
    )
    mutation_threshold: float = Field(
        default=0.15,
        ge=0,
        le=1,
        description="Minimum mutation ratio used to keep genes for plotting.",
    )
    cohort: str = Field(
        default=DEFAULT_V2_COHORT,
        description="Cohort name used for default output filenames.",
    )
    species: str = Field(
        default="human",
        description="Gene-symbol map species. Use human or mouse.",
    )
    gene_symbol_map_path: Optional[str] = Field(
        default=None,
        description=(
            "Optional TSV with columns gene_id and gene_symbol. "
            "When omitted, species=human uses the generated human map and "
            "species=mouse uses the generated mouse map."
        ),
    )
    heatmap_filename: Optional[str] = Field(
        default=None,
        description="Optional explicit heatmap PDF filename. When omitted, derived from cohort.",
    )
    mutation_type_filename: Optional[str] = Field(
        default=None,
        description="Optional explicit mutation-type distribution PDF filename. When omitted, derived from cohort.",
    )
    output_prefix: Optional[str] = Field(
        default=None,
        description="Optional prefix added to generated PDF filenames.",
    )


class MutationHeatmapV2Response(HnscMutationFiguresResponse):
    mutation_binary_path: str
    meta_path: str
    gene_symbol_map_path: str
    species: str
