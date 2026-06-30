from __future__ import annotations

import argparse

from omicsone.replay.pathway_scatter import run_phosphosite_protein_pathway_pipeline


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate phosphosite-protein pathway tables and scatter plots."
    )
    parser.add_argument("--config", required=True, help="Path to pipeline config.ini.")
    args = parser.parse_args()

    run_phosphosite_protein_pathway_pipeline(args.config)


if __name__ == "__main__":
    main()

