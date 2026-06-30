from __future__ import annotations

import argparse
from typing import Sequence

from omicsone.cli.common import add_config_run_arguments, print_json_enabled


def run_differential(args: argparse.Namespace) -> int:
    from omicsone.replay.differential import run_differential_analysis

    run_differential_analysis(
        args.config,
        output_dir=args.output_dir,
        print_json=print_json_enabled(args),
    )
    return 0


def run_phospho_differential(args: argparse.Namespace) -> int:
    from omicsone.replay.differential import run_phospho_differential_analysis

    run_phospho_differential_analysis(
        args.config,
        output_dir=args.output_dir,
        print_json=print_json_enabled(args),
    )
    return 0


def run_mutations(args: argparse.Namespace) -> int:
    from omicsone.replay.mutations import run_mutation_figures

    run_mutation_figures(
        args.config,
        output_dir=args.output_dir,
        print_json=print_json_enabled(args),
    )
    return 0


def post_mutations_api(args: argparse.Namespace) -> int:
    from omicsone.replay.mutations import post_mutation_figures_api

    post_mutation_figures_api(
        args.config,
        output_dir=args.output_dir,
        api_url=args.api_url,
        print_json=print_json_enabled(args),
    )
    return 0


def run_boxplots(args: argparse.Namespace) -> int:
    from omicsone.replay.boxplots import run_boxplot_figures

    run_boxplot_figures(
        args.config,
        output_dir=args.output_dir,
        print_json=print_json_enabled(args),
    )
    return 0


def run_pathway_scatter(args: argparse.Namespace) -> int:
    from omicsone.replay.pathway_scatter import run_pathway_scatter_plots

    run_pathway_scatter_plots(
        args.config,
        output_dir=args.output_dir,
        print_json=print_json_enabled(args),
    )
    return 0


def run_phosphosite_protein_pathway(args: argparse.Namespace) -> int:
    from omicsone.replay.pathway_scatter import run_phosphosite_protein_pathway_pipeline

    run_phosphosite_protein_pathway_pipeline(
        args.config,
        output_dir=args.output_dir,
        print_json=print_json_enabled(args),
    )
    return 0


def run_cnv_correlation(args: argparse.Namespace) -> int:
    from omicsone.replay.cnv_correlation_pipeline import run_cnv_correlation_pipeline

    run_cnv_correlation_pipeline(
        args.config,
        output_dir=args.output_dir,
        print_json=print_json_enabled(args),
    )
    return 0


def run_app(_args: argparse.Namespace) -> int:
    from omicsone_streamlit.__main__ import main as streamlit_main

    result = streamlit_main()
    return int(result or 0)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="omicsone",
        description="Run OmicsOne local workflows and utilities.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    differential = subparsers.add_parser("differential", help="Differential analysis workflows.")
    differential_subparsers = differential.add_subparsers(dest="differential_command", required=True)
    differential_run = differential_subparsers.add_parser("run", help="Run differential volcano/enrichment analysis.")
    add_config_run_arguments(differential_run)
    differential_run.set_defaults(func=run_differential)
    differential_phospho = differential_subparsers.add_parser(
        "phospho",
        help="Run phosphosite-level differential volcano/enrichment analysis.",
    )
    add_config_run_arguments(differential_phospho)
    differential_phospho.set_defaults(func=run_phospho_differential)

    mutations = subparsers.add_parser("mutations", help="Mutation figure workflows.")
    mutations_subparsers = mutations.add_subparsers(dest="mutations_command", required=True)
    mutations_run = mutations_subparsers.add_parser("run", help="Run mutation figure generation locally.")
    add_config_run_arguments(mutations_run)
    mutations_run.set_defaults(func=run_mutations)
    mutations_post_api = mutations_subparsers.add_parser("post-api", help="POST mutation figure payload to OmicsOne API.")
    add_config_run_arguments(mutations_post_api)
    mutations_post_api.add_argument("--api-url", help="Override the configured API endpoint URL.")
    mutations_post_api.set_defaults(func=post_mutations_api)

    boxplots = subparsers.add_parser("boxplots", help="Boxplot figure workflows.")
    boxplots_subparsers = boxplots.add_subparsers(dest="boxplots_command", required=True)
    boxplots_run = boxplots_subparsers.add_parser("run", help="Run normal-vs-tumor boxplot generation.")
    add_config_run_arguments(boxplots_run)
    boxplots_run.set_defaults(func=run_boxplots)

    pathway = subparsers.add_parser("pathway-scatter", help="Pathway scatter workflows.")
    pathway_subparsers = pathway.add_subparsers(dest="pathway_scatter_command", required=True)
    pathway_run = pathway_subparsers.add_parser("run", help="Run protein-vs-phosphosite pathway scatter plots.")
    add_config_run_arguments(pathway_run)
    pathway_run.set_defaults(func=run_pathway_scatter)
    pathway_pipeline = pathway_subparsers.add_parser(
        "phosphosite-protein",
        help="Generate phosphosite-protein pathway tables and scatter plots.",
    )
    add_config_run_arguments(pathway_pipeline)
    pathway_pipeline.set_defaults(func=run_phosphosite_protein_pathway)

    cnv = subparsers.add_parser("cnv-correlation", help="CNV correlation workflows.")
    cnv_subparsers = cnv.add_subparsers(dest="cnv_correlation_command", required=True)
    cnv_run = cnv_subparsers.add_parser("run", help="Run the CNV correlation pipeline.")
    add_config_run_arguments(cnv_run)
    cnv_run.set_defaults(func=run_cnv_correlation)

    app = subparsers.add_parser("app", help="Launch the Streamlit UI app.")
    app.set_defaults(func=run_app)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


