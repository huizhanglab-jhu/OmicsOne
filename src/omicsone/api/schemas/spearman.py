from typing import Literal

from pydantic import BaseModel, Field, model_validator


SpearmanOmicsDataType = Literal["cnv", "rna", "protein"]


class SpearmanOmicsInput(BaseModel):
    path: str = Field(
        description="Path to an omics matrix file with an idx column and sample columns.",
    )
    data_type: SpearmanOmicsDataType = Field(
        description="Omics data type. Supported values: cnv, rna, protein.",
    )


class SpearmanComputeFileRequest(BaseModel):
    input_file1: str = Field(
        description="Path to the first whitespace-delimited numeric matrix file.",
    )
    input_file2: str = Field(
        description="Path to the second whitespace-delimited numeric matrix file.",
    )
    output_file: str = Field(
        description="Path where row-pair Spearman correlations will be written.",
    )
    min_valid_pairs: int = Field(
        default=2,
        ge=1,
        description="Minimum paired non-NaN values required to compute a correlation.",
    )


class SpearmanComputeFileResponse(BaseModel):
    input_file1: str
    input_file2: str
    output_file: str
    min_valid_pairs: int
    backend: str
    result_rows: int


class SpearmanPairedOmicsRequest(BaseModel):
    input1: SpearmanOmicsInput
    input2: SpearmanOmicsInput
    input3: SpearmanOmicsInput | None = Field(
        default=None,
        description=(
            "Optional third omics matrix. When provided, common genes and samples "
            "are intersected across all three inputs and all pairwise Spearman "
            "correlations are generated."
        ),
    )
    output_dir: str = Field(
        description=(
            "Output folder where *_for_corr.txt files and the Spearman correlation "
            "output file will be written."
        ),
    )
    output_prefix: str | None = Field(
        default=None,
        description=(
            "Optional prefix added to generated Spearman correlation filenames."
        ),
    )
    min_valid_pairs: int = Field(
        default=2,
        ge=1,
        description="Minimum paired non-NaN values required to compute a correlation.",
    )

    @model_validator(mode="after")
    def validate_unique_data_types(self):
        data_types = [self.input1.data_type, self.input2.data_type]
        if self.input3 is not None:
            data_types.append(self.input3.data_type)
        if len(set(data_types)) != len(data_types):
            raise ValueError("Each omics input must have a unique data_type")
        return self


class SpearmanPairResponse(BaseModel):
    data_type1: str
    data_type2: str
    matrix_file1: str
    matrix_file2: str
    correlation_file: str
    result_rows: int


class SpearmanPairedOmicsResponse(BaseModel):
    input_files: dict[str, str]
    output_dir: str
    matrix_files: dict[str, str]
    common_gene_count: int
    common_sample_count: int
    matrix_shapes: dict[str, list[int]]
    min_valid_pairs: int
    backend: str
    pairs: list[SpearmanPairResponse]
