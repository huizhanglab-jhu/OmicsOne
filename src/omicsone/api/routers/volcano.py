from fastapi import APIRouter, HTTPException

from omicsone.api.schemas.volcano import (
    VolcanoEnrichmentRequest,
    VolcanoEnrichmentResponse,
    VolcanoPresetRequest,
)
from omicsone.services.volcano_enrichment import (
    generate_volcano_enrichment,
    resolve_volcano_preset,
)


router = APIRouter()


def _response_from_result(result) -> dict:
    return {
        "cohort": result.cohort,
        "omics": result.omics,
        "output_dir": str(result.output_dir),
        "diff_tsv": str(result.diff_tsv),
        "combined_matrix_tsv": str(result.combined_matrix_tsv),
        "up_genes_tsv": str(result.up_genes_tsv),
        "down_genes_tsv": str(result.down_genes_tsv),
        "total_genes_tsv": str(result.total_genes_tsv),
        "up_enrichment_tsv": str(result.up_enrichment_tsv),
        "down_enrichment_tsv": str(result.down_enrichment_tsv),
        "enrichment_plot_tsv": str(result.enrichment_plot_tsv),
        "volcano_png": str(result.volcano_png),
        "volcano_pdf": str(result.volcano_pdf),
        "volcano_tiff": str(result.volcano_tiff),
        "enrichment_png": str(result.enrichment_png),
        "enrichment_pdf": str(result.enrichment_pdf),
        "enrichment_tiff": str(result.enrichment_tiff),
        "report_html": str(result.report_html),
        "result_log": str(result.result_log),
        "n8n_js": str(result.n8n_js),
        "feature_count": result.feature_count,
        "diff_feature_count": result.diff_feature_count,
        "up_count": result.up_count,
        "down_count": result.down_count,
        "pure_up_gene_count": result.pure_up_gene_count,
        "pure_down_gene_count": result.pure_down_gene_count,
        "total_gene_count": result.total_gene_count,
        "up_enrichment_count": result.up_enrichment_count,
        "down_enrichment_count": result.down_enrichment_count,
        "method": result.method,
        "fdr_cutoff": result.fdr_cutoff,
        "log2fc_cutoff": result.log2fc_cutoff,
        "gene_sets": result.gene_sets,
    }


@router.post("/volcano/enrichment", response_model=VolcanoEnrichmentResponse)
def create_volcano_enrichment(request: VolcanoEnrichmentRequest):
    try:
        result = generate_volcano_enrichment(**request.model_dump())
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    return _response_from_result(result)


@router.post("/volcano/enrichment/preset", response_model=VolcanoEnrichmentResponse)
def create_volcano_enrichment_from_preset(request: VolcanoPresetRequest):
    try:
        payload = resolve_volcano_preset(request.preset)
        payload.update(request.overrides)
        validated = VolcanoEnrichmentRequest(**payload)
        result = generate_volcano_enrichment(**validated.model_dump())
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    return _response_from_result(result)
