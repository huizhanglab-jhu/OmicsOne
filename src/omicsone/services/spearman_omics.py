from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from omicsone.services.spearman_preprocessing import write_rust_spearman_input
from omicsone.utils import spearmanr


OmicsDataType = Literal["cnv", "rna", "protein"]


@dataclass(frozen=True)
class OmicsInput:
    path: Path
    data_type: OmicsDataType


@dataclass(frozen=True)
class SpearmanPairResult:
    data_type1: str
    data_type2: str
    matrix_file1: Path
    matrix_file2: Path
    correlation_file: Path
    result_rows: int


@dataclass(frozen=True)
class PairedOmicsSpearmanResult:
    inputs: tuple[OmicsInput, ...]
    output_dir: Path
    matrix_files: dict[str, Path]
    pairs: list[SpearmanPairResult]
    common_gene_count: int
    common_sample_count: int
    matrix_shapes: dict[str, tuple[int, int]]
    min_valid_pairs: int
    backend: str


def compute_paired_omics_spearman(
    inputs: list[OmicsInput],
    output_dir: Path,
    min_valid_pairs: int = 2,
    output_prefix: str | None = None,
    universe_inputs: list[OmicsInput] | None = None,
) -> PairedOmicsSpearmanResult:
    if len(inputs) not in (2, 3):
        raise ValueError("Provide either two or three omics inputs")
    if min_valid_pairs < 1:
        raise ValueError("min_valid_pairs must be at least 1")

    data_types = [item.data_type for item in inputs]
    if len(set(data_types)) != len(data_types):
        raise ValueError("Each omics input must have a unique data_type")

    output_dir.mkdir(parents=True, exist_ok=True)

    matrices = {
        item.data_type: _read_and_preprocess_matrix(item.path, item.data_type)
        for item in inputs
    }

    universe_matrices = matrices
    if universe_inputs is not None:
        universe_data_types = [item.data_type for item in universe_inputs]
        if len(set(universe_data_types)) != len(universe_data_types):
            raise ValueError("Each universe input must have a unique data_type")
        universe_matrices = {}
        for item in universe_inputs:
            if item.data_type in matrices:
                universe_matrices[item.data_type] = matrices[item.data_type]
            else:
                universe_matrices[item.data_type] = _read_and_preprocess_matrix(
                    item.path,
                    item.data_type,
                )

    common_genes = sorted(
        set.intersection(*(set(df.index) for df in universe_matrices.values()))
    )
    common_samples = sorted(
        set.intersection(*(set(df.columns) for df in universe_matrices.values()))
    )

    if not common_genes:
        raise ValueError("No common genes found across the provided input matrices")
    if not common_samples:
        raise ValueError("No common samples found across the provided input matrices")

    matrix_files = {}
    rust_matrix_files = {}
    matrix_shapes = {}
    rust_input_dir = output_dir / "_rust_spearman_inputs"

    for data_type, df in matrices.items():
        matrix = df.loc[common_genes, common_samples]
        matrix = matrix[~matrix.index.duplicated(keep="first")]
        matrix_shapes[data_type] = matrix.shape

        matrix_file = output_dir / _for_corr_filename(data_type)
        matrix.to_csv(matrix_file, header=True, index=True, index_label="idx", sep="\t")
        matrix_files[data_type] = matrix_file

        rust_matrix_file = rust_input_dir / _for_rust_spearman_filename(data_type)
        rust_matrix_files[data_type] = write_rust_spearman_input(
            matrix_file,
            rust_matrix_file,
            has_header=True,
            has_index=True,
            sep="\t",
        )

    pairs = []
    for data_type1, data_type2 in combinations(data_types, 2):
        correlation_file = output_dir / _correlation_filename(
            data_type1,
            data_type2,
            output_prefix=output_prefix,
        )
        result_rows = spearmanr.compute_file(
            rust_matrix_files[data_type1],
            rust_matrix_files[data_type2],
            correlation_file,
            min_valid_pairs=min_valid_pairs,
        )
        pairs.append(
            SpearmanPairResult(
                data_type1=data_type1,
                data_type2=data_type2,
                matrix_file1=matrix_files[data_type1],
                matrix_file2=matrix_files[data_type2],
                correlation_file=correlation_file,
                result_rows=result_rows,
            )
        )

    return PairedOmicsSpearmanResult(
        inputs=tuple(inputs),
        output_dir=output_dir,
        matrix_files=matrix_files,
        pairs=pairs,
        common_gene_count=len(common_genes),
        common_sample_count=len(common_samples),
        matrix_shapes=matrix_shapes,
        min_valid_pairs=min_valid_pairs,
        backend=spearmanr.backend(),
    )


def _read_and_preprocess_matrix(path: Path, data_type: OmicsDataType) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Input file does not exist: {path}")

    df = pd.read_csv(path, sep="\t")
    if "idx" not in df.columns:
        raise ValueError(f"Input file must contain an 'idx' column: {path}")

    df = df.set_index("idx")
    df.index = [str(value).split(".")[0] for value in df.index]
    df = df.apply(pd.to_numeric, errors="coerce")

    if data_type == "rna":
        df = df.replace(0, np.nan).dropna()
    elif data_type == "protein":
        df = df.dropna()
    elif data_type == "cnv":
        pass
    else:
        raise ValueError(f"Unsupported omics data type: {data_type}")

    df = df.loc[~df.index.duplicated(keep="first"), ~df.columns.duplicated(keep="first")]
    return df


def _for_corr_filename(data_type: OmicsDataType) -> str:
    if data_type == "protein":
        return "pro_for_corr.txt"
    return f"{data_type}_for_corr.txt"


def _for_rust_spearman_filename(data_type: OmicsDataType) -> str:
    if data_type == "protein":
        return "pro_for_rust_spearman.txt"
    return f"{data_type}_for_rust_spearman.txt"


def _correlation_filename(
    data_type1: str,
    data_type2: str,
    output_prefix: str | None = None,
) -> str:
    prefix = f"{output_prefix}_" if output_prefix else ""
    return f"{prefix}{data_type1}_vs_{data_type2}_spearman_correlations.txt"
