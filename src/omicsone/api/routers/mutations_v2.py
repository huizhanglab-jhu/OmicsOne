from fastapi import APIRouter, HTTPException

from omicsone.api.schemas.mutations import (
    MutationHeatmapV2Request,
    MutationHeatmapV2Response,
)
from omicsone.services.mutation_figures import generate_mutation_figures_from_binary


router = APIRouter()


@router.post("/heatmap/figures", response_model=MutationHeatmapV2Response)
def create_mutation_heatmap_figures_v2(request: MutationHeatmapV2Request):
    try:
        result = generate_mutation_figures_from_binary(
            mutation_binary_path=request.mutation_binary_path,
            meta_path=request.meta_path,
            maf_path=request.maf_path,
            output_dir=request.output_dir,
            mutation_threshold=request.mutation_threshold,
            cohort=request.cohort,
            species=request.species,
            gene_symbol_map_path=request.gene_symbol_map_path,
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
        "mutation_binary_path": str(result.mutation_binary_path),
        "meta_path": str(result.meta_path),
        "gene_symbol_map_path": str(result.gene_symbol_map_path),
        "species": result.species,
    }
