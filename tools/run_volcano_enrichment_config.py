from __future__ import annotations

import argparse

from omicsone.services.volcano_enrichment import generate_volcano_enrichment_from_config


def main() -> None:
    parser = argparse.ArgumentParser(description="Run volcano and enrichment analysis from config.ini.")
    parser.add_argument("--config", required=True, help="Path to config.ini.")
    args = parser.parse_args()

    result = generate_volcano_enrichment_from_config(args.config)
    print(
        {
            "output_dir": str(result.output_dir),
            "diff_feature_count": result.diff_feature_count,
            "up_count": result.up_count,
            "down_count": result.down_count,
            "pure_up_gene_count": result.pure_up_gene_count,
            "pure_down_gene_count": result.pure_down_gene_count,
            "up_enrichment_count": result.up_enrichment_count,
            "down_enrichment_count": result.down_enrichment_count,
            "volcano_png": str(result.volcano_png),
            "volcano_pdf": str(result.volcano_pdf),
            "volcano_tiff": str(result.volcano_tiff),
            "enrichment_png": str(result.enrichment_png),
            "enrichment_pdf": str(result.enrichment_pdf),
            "enrichment_tiff": str(result.enrichment_tiff),
            "report_html": str(result.report_html),
            "result_log": str(result.result_log),
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
