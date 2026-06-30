from __future__ import annotations

import argparse

from omicsone.services.pathway_scatter import run_pathway_scatter_analysis_from_config


def main() -> None:
    parser = argparse.ArgumentParser(description="Run protein/phosphosite pathway scatter plots.")
    parser.add_argument("--config", required=True, help="Path to pathway scatter config.ini.")
    args = parser.parse_args()

    results = run_pathway_scatter_analysis_from_config(args.config)
    for result in results:
        print(
            {
                "pathway": result.pathway,
                "output_dir": str(result.output_dir),
                "point_count": result.point_count,
                "highlight_count": result.highlight_count,
                "trend_slope": result.trend_slope,
                "trend_intercept": result.trend_intercept,
                "trend_r": result.trend_r,
                "png": str(result.png),
                "pdf": str(result.pdf),
                "tiff": str(result.tiff),
                "missing_highlights_tsv": str(result.missing_highlights_tsv),
            },
            flush=True,
        )


if __name__ == "__main__":
    main()
