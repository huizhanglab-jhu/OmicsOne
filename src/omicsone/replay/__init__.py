"""Replay adapters for OmicsOne analysis services."""

from __future__ import annotations

from importlib import import_module
from typing import Any


_EXPORTS = {
    "DEFAULT_OUTPUT_CONFIG_NAME": "omicsone.replay.differential",
    "default_boxplot_config_path": "omicsone.replay.boxplots",
    "default_cnv_correlation_pipeline_config_path": "omicsone.replay.cnv_correlation_pipeline",
    "default_differential_config_path": "omicsone.replay.differential",
    "default_mutation_config_path": "omicsone.replay.mutations",
    "default_pathway_scatter_config_path": "omicsone.replay.pathway_scatter",
    "ensure_output_boxplot_config": "omicsone.replay.boxplots",
    "ensure_output_cnv_correlation_pipeline_config": "omicsone.replay.cnv_correlation_pipeline",
    "ensure_output_differential_config": "omicsone.replay.differential",
    "ensure_output_mutation_config": "omicsone.replay.mutations",
    "ensure_output_pathway_scatter_config": "omicsone.replay.pathway_scatter",
    "load_boxplot_config": "omicsone.replay.boxplots",
    "load_cnv_correlation_pipeline_config": "omicsone.replay.cnv_correlation_pipeline",
    "load_mutation_config": "omicsone.replay.mutations",
    "load_pathway_scatter_replay_config": "omicsone.replay.pathway_scatter",
    "post_mutation_figures_api": "omicsone.replay.mutations",
    "run_boxplot_figures": "omicsone.replay.boxplots",
    "run_cnv_correlation_pipeline": "omicsone.replay.cnv_correlation_pipeline",
    "run_cnv_correlation_clean_figures": "omicsone.replay.cnv_correlation_clean_figures",
    "run_differential_analysis": "omicsone.replay.differential",
    "run_mutation_figures": "omicsone.replay.mutations",
    "run_pathway_scatter_plots": "omicsone.replay.pathway_scatter",
    "run_phospho_differential_analysis": "omicsone.replay.differential",
    "run_phosphosite_protein_pathway_pipeline": "omicsone.replay.pathway_scatter",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str) -> Any:
    module_name = _EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(module_name)
    value = getattr(module, name)
    globals()[name] = value
    return value

