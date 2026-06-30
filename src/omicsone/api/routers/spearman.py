from pathlib import Path

from fastapi import APIRouter, HTTPException

from omicsone.api.schemas.spearman import (
    SpearmanComputeFileRequest,
    SpearmanComputeFileResponse,
    SpearmanPairedOmicsRequest,
    SpearmanPairedOmicsResponse,
)
from omicsone.services.spearman_omics import OmicsInput, compute_paired_omics_spearman
from omicsone.utils import spearmanr


router = APIRouter()


@router.post("/compute-file", response_model=SpearmanComputeFileResponse)
def compute_spearman_file(request: SpearmanComputeFileRequest):
    input_file1 = Path(request.input_file1).expanduser()
    input_file2 = Path(request.input_file2).expanduser()
    output_file = Path(request.output_file).expanduser()

    if not input_file1.is_file():
        raise HTTPException(
            status_code=404,
            detail=f"Input file 1 does not exist: {input_file1}",
        )
    if not input_file2.is_file():
        raise HTTPException(
            status_code=404,
            detail=f"Input file 2 does not exist: {input_file2}",
        )
    if output_file.parent and not output_file.parent.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Output directory does not exist: {output_file.parent}",
        )

    try:
        result_rows = spearmanr.compute_file(
            input_file1,
            input_file2,
            output_file,
            min_valid_pairs=request.min_valid_pairs,
        )
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except OSError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    return {
        "input_file1": str(input_file1),
        "input_file2": str(input_file2),
        "output_file": str(output_file),
        "min_valid_pairs": request.min_valid_pairs,
        "backend": spearmanr.backend(),
        "result_rows": result_rows,
    }


@router.post("/paired-omics", response_model=SpearmanPairedOmicsResponse)
def compute_paired_omics_spearman_file(request: SpearmanPairedOmicsRequest):
    try:
        inputs = [
            OmicsInput(
                path=Path(request.input1.path).expanduser(),
                data_type=request.input1.data_type,
            ),
            OmicsInput(
                path=Path(request.input2.path).expanduser(),
                data_type=request.input2.data_type,
            ),
        ]
        if request.input3 is not None:
            inputs.append(
                OmicsInput(
                    path=Path(request.input3.path).expanduser(),
                    data_type=request.input3.data_type,
                )
            )

        result = compute_paired_omics_spearman(
            inputs=inputs,
            output_dir=Path(request.output_dir).expanduser(),
            min_valid_pairs=request.min_valid_pairs,
            output_prefix=request.output_prefix,
        )
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except OSError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    return {
        "input_files": {
            item.data_type: str(item.path)
            for item in result.inputs
        },
        "output_dir": str(result.output_dir),
        "matrix_files": {
            data_type: str(path)
            for data_type, path in result.matrix_files.items()
        },
        "common_gene_count": result.common_gene_count,
        "common_sample_count": result.common_sample_count,
        "matrix_shapes": {
            data_type: list(shape)
            for data_type, shape in result.matrix_shapes.items()
        },
        "min_valid_pairs": result.min_valid_pairs,
        "backend": result.backend,
        "pairs": [
            {
                "data_type1": pair.data_type1,
                "data_type2": pair.data_type2,
                "matrix_file1": str(pair.matrix_file1),
                "matrix_file2": str(pair.matrix_file2),
                "correlation_file": str(pair.correlation_file),
                "result_rows": pair.result_rows,
            }
            for pair in result.pairs
        ],
    }
