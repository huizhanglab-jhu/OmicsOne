from pathlib import Path

from fastapi import APIRouter, BackgroundTasks, HTTPException

from omicsone.api.schemas.cnv_correlation import (
    CnvCorrelationFiguresRequest,
    CnvCorrelationFigureSetResponse,
)
from omicsone.api.schemas.cnv_correlation_pipeline import (
    CnvCorrelationPipelineRequest,
    CnvCorrelationPipelineStartResponse,
    CnvCorrelationPipelineStatusResponse,
)
from omicsone.services.cnv_correlation_pipeline import (
    create_pipeline_job,
    read_pipeline_job,
    run_pipeline_job,
)
from omicsone.plots.cnv_correlation import (
    CnvCorrelationFigureResult,
    generate_cnv_correlation_figures,
)


router = APIRouter()


@router.post("/pipeline", response_model=CnvCorrelationPipelineStartResponse)
def start_cnv_correlation_pipeline(
    request: CnvCorrelationPipelineRequest,
    background_tasks: BackgroundTasks,
):
    try:
        record = create_pipeline_job(request.model_dump())
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except OSError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    background_tasks.add_task(run_pipeline_job, record["job_id"])
    return {
        "job_id": record["job_id"],
        "status": record["status"],
        "status_url": f"/api/v1/cnv-correlation/pipeline/{record['job_id']}",
        "output_dir": record["output_dir"],
        "current_step": record["current_step"],
    }


@router.get(
    "/pipeline/{job_id}",
    response_model=CnvCorrelationPipelineStatusResponse,
)
def get_cnv_correlation_pipeline_status(job_id: str):
    try:
        record = read_pipeline_job(job_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    return {
        "job_id": record["job_id"],
        "status": record["status"],
        "created_at": record["created_at"],
        "started_at": record["started_at"],
        "ended_at": record["ended_at"],
        "current_step": record["current_step"],
        "output_dir": record["output_dir"],
        "html_report": record["html_report"],
        "json_report": record["json_report"],
        "error": record["error"],
        "jobs": record["jobs"],
    }


@router.post("/figures", response_model=CnvCorrelationFigureSetResponse)
def create_cnv_correlation_figures(request: CnvCorrelationFiguresRequest):
    try:
        result = generate_cnv_correlation_figures(
            correlation_file=Path(request.correlation_file).expanduser(),
            cnv_for_corr_file=Path(request.cnv_for_corr_file).expanduser(),
            target_for_corr_file=Path(request.target_for_corr_file).expanduser(),
            fasta_file=Path(request.fasta_file).expanduser(),
            chromosomes_file=Path(request.chromosomes_file).expanduser(),
            cytoband_file=Path(request.cytoband_file).expanduser(),
            output_dir=Path(request.output_dir).expanduser(),
            gistic_file=(
                Path(request.gistic_file).expanduser()
                if request.gistic_file is not None
                else None
            ),
            target_type=request.target_type,
            correlation_threshold=request.correlation_threshold,
            output_prefix=request.output_prefix,
            chunksize=request.chunksize,
            dpi=request.dpi,
        )
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except OSError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    return _response_from_result(result)


def _response_from_result(result: CnvCorrelationFigureResult) -> dict:
    png_files = {
        "corr_heatmap": result.corr_heatmap_file,
        "cnv_distribution": result.cnv_distribution_file,
        "chromosome": result.chromosome_file,
        "gistic": result.gistic_file,
        "combined": result.combined_file,
    }

    return {
        "output_dir": str(result.output_dir),
        "annotated_correlations_file": str(result.annotated_correlations_file),
        "cnv_distribution_counts_file": str(result.cnv_distribution_counts_file),
        "gistic_counts_file": (
            str(result.gistic_counts_file)
            if result.gistic_counts_file is not None
            else None
        ),
        "png_files": _stringify_paths(png_files),
        "pdf_files": _stringify_paths(
            {
                name: path.with_suffix(".pdf") if path is not None else None
                for name, path in png_files.items()
            }
        ),
        "tiff_files": _stringify_paths(
            {
                name: path.with_suffix(".tiff") if path is not None else None
                for name, path in png_files.items()
            }
        ),
        "filtered_correlation_count": result.filtered_correlation_count,
        "annotated_correlation_count": result.annotated_correlation_count,
    }


def _stringify_paths(paths: dict[str, Path | None]) -> dict[str, str | None]:
    return {
        name: str(path) if path is not None else None
        for name, path in paths.items()
    }
