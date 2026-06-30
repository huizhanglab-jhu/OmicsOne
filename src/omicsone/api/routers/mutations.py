from fastapi import APIRouter, HTTPException

from omicsone.api.schemas.mutations import (
    HnscMutationFiguresRequest,
    HnscMutationFiguresResponse,
)
from omicsone.services.mutation_figures import generate_hnsc_mutation_figures


router = APIRouter()


@router.post("/heatmap/figures", response_model=HnscMutationFiguresResponse)
def create_hnsc_mutation_figures(request: HnscMutationFiguresRequest):
    try:
        result = generate_hnsc_mutation_figures(
            mutation_excel_path=request.mutation_excel_path,
            maf_path=request.maf_path,
            output_dir=request.output_dir,
            mutation_threshold=request.mutation_threshold,
            cohort=request.cohort,
            heatmap_filename=request.heatmap_filename,
            mutation_type_filename=request.mutation_type_filename,
            output_prefix=request.output_prefix,
        )
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    return {
        "heatmap_pdf": str(result.heatmap_pdf),
        "mutation_type_pdf": str(result.mutation_type_pdf),
        "result_log": str(result.result_log),
        "gene_summary_tsv": str(result.gene_summary_tsv),
        "binary_matrix_tsv": str(result.binary_matrix_tsv),
        "mutation_type_matrix_tsv": str(result.mutation_type_matrix_tsv),
        "encoded_matrix_tsv": str(result.encoded_matrix_tsv),
        "sample_annotations_tsv": str(result.sample_annotations_tsv),
        "sample_annotation_colors_tsv": str(result.sample_annotation_colors_tsv),
        "mutation_color_table_tsv": str(result.mutation_color_table_tsv),
        "filtered_maf_tsv": str(result.filtered_maf_tsv),
        "gene_count": result.gene_count,
        "sample_count": result.sample_count,
        "cohort": result.cohort,
        "total_maf_rows": result.total_maf_rows,
        "filtered_maf_rows": result.filtered_maf_rows,
        "found_mutations": result.found_mutations,
        "heatmap_width_inch": result.heatmap_width_inch,
        "heatmap_height_inch": result.heatmap_height_inch,
        "heatmap_aspect": result.heatmap_aspect,
        "mutation_type_width_inch": result.mutation_type_width_inch,
        "mutation_type_height_inch": result.mutation_type_height_inch,
        "mutation_type_aspect": result.mutation_type_aspect,
        "sample_gene_ratio": result.sample_gene_ratio,
    }


@router.post("/hnsc/figures", response_model=HnscMutationFiguresResponse, include_in_schema=False)
def create_hnsc_mutation_figures_legacy(request: HnscMutationFiguresRequest):
    return create_hnsc_mutation_figures(request)
