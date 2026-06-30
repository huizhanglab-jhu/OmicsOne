from fastapi import APIRouter, HTTPException

from omicsone.api.schemas.boxplots import (
    BoxplotFiguresRequest,
    BoxplotFiguresResponse,
)
from omicsone.services.boxplot_figures import generate_boxplot_figures


router = APIRouter()


@router.post("/boxplot/figures", response_model=BoxplotFiguresResponse)
def create_boxplot_figures(request: BoxplotFiguresRequest):
    try:
        result = generate_boxplot_figures(**request.model_dump())
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    return {
        "cohort": result.cohort,
        "omics": result.omics,
        "output_dir": str(result.output_dir),
        "boxplot_pdfs": [str(path) for path in result.boxplot_pdfs],
        "summary_tsv": str(result.summary_tsv),
        "result_log": str(result.result_log),
        "n8n_js": str(result.n8n_js),
        "generated_count": result.generated_count,
        "missing_genes": result.missing_genes,
        "records": [
            {
                "gene": record.gene,
                "gene_id": record.gene_id,
                "pdf": str(record.pdf),
                "pvalue": record.pvalue,
                "significance": record.significance,
                "normal_count": record.normal_count,
                "tumor_count": record.tumor_count,
            }
            for record in result.records
        ],
    }
