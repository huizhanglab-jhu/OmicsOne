from __future__ import annotations

import sys
import types

from omicsone.cli import commands


def _install_module(monkeypatch, name: str, function_name: str, calls: list[tuple[str, tuple, dict]]):
    module = types.ModuleType(name)

    def fake_function(*args, **kwargs):
        calls.append((function_name, args, kwargs))
        return {"ok": True}

    setattr(module, function_name, fake_function)
    monkeypatch.setitem(sys.modules, name, module)


def test_differential_run_dispatches_to_replay_adapter(monkeypatch):
    calls = []
    _install_module(monkeypatch, "omicsone.replay.differential", "run_differential_analysis", calls)

    assert commands.main(["differential", "run", "--config", "input.ini", "--output-dir", "out", "--quiet"]) == 0

    assert calls == [
        (
            "run_differential_analysis",
            ("input.ini",),
            {"output_dir": "out", "print_json": False},
        )
    ]


def test_phospho_differential_dispatches_to_replay_adapter(monkeypatch):
    calls = []
    _install_module(monkeypatch, "omicsone.replay.differential", "run_phospho_differential_analysis", calls)

    assert commands.main(["differential", "phospho", "--config", "input.ini"]) == 0

    assert calls == [
        (
            "run_phospho_differential_analysis",
            ("input.ini",),
            {"output_dir": None, "print_json": True},
        )
    ]


def test_mutations_run_dispatches_to_replay_adapter(monkeypatch):
    calls = []
    _install_module(monkeypatch, "omicsone.replay.mutations", "run_mutation_figures", calls)

    assert commands.main(["mutations", "run", "--config", "input.ini", "--quiet"]) == 0

    assert calls == [
        (
            "run_mutation_figures",
            ("input.ini",),
            {"output_dir": None, "print_json": False},
        )
    ]


def test_mutations_post_api_dispatches_to_replay_adapter(monkeypatch):
    calls = []
    _install_module(monkeypatch, "omicsone.replay.mutations", "post_mutation_figures_api", calls)

    assert commands.main(
        [
            "mutations",
            "post-api",
            "--config",
            "input.ini",
            "--api-url",
            "http://127.0.0.1:8001/api/v1/mutations/heatmap/figures",
        ]
    ) == 0

    assert calls == [
        (
            "post_mutation_figures_api",
            ("input.ini",),
            {
                "output_dir": None,
                "api_url": "http://127.0.0.1:8001/api/v1/mutations/heatmap/figures",
                "print_json": True,
            },
        )
    ]


def test_boxplots_run_dispatches_to_replay_adapter(monkeypatch):
    calls = []
    _install_module(monkeypatch, "omicsone.replay.boxplots", "run_boxplot_figures", calls)

    assert commands.main(["boxplots", "run", "--config", "input.ini"]) == 0

    assert calls == [
        (
            "run_boxplot_figures",
            ("input.ini",),
            {"output_dir": None, "print_json": True},
        )
    ]


def test_pathway_scatter_run_dispatches_to_replay_adapter(monkeypatch):
    calls = []
    _install_module(monkeypatch, "omicsone.replay.pathway_scatter", "run_pathway_scatter_plots", calls)

    assert commands.main(["pathway-scatter", "run", "--config", "input.ini"]) == 0

    assert calls == [
        (
            "run_pathway_scatter_plots",
            ("input.ini",),
            {"output_dir": None, "print_json": True},
        )
    ]


def test_pathway_scatter_pipeline_dispatches_to_replay_adapter(monkeypatch):
    calls = []
    _install_module(monkeypatch, "omicsone.replay.pathway_scatter", "run_phosphosite_protein_pathway_pipeline", calls)

    assert commands.main(["pathway-scatter", "phosphosite-protein", "--config", "input.ini"]) == 0

    assert calls == [
        (
            "run_phosphosite_protein_pathway_pipeline",
            ("input.ini",),
            {"output_dir": None, "print_json": True},
        )
    ]


def test_cnv_correlation_run_dispatches_to_replay_adapter(monkeypatch):
    calls = []
    _install_module(monkeypatch, "omicsone.replay.cnv_correlation_pipeline", "run_cnv_correlation_pipeline", calls)

    assert commands.main(["cnv-correlation", "run", "--config", "input.ini", "--output-dir", "out"]) == 0

    assert calls == [
        (
            "run_cnv_correlation_pipeline",
            ("input.ini",),
            {"output_dir": "out", "print_json": True},
        )
    ]


def test_app_dispatches_to_streamlit_launcher(monkeypatch):
    module = types.ModuleType("omicsone_streamlit.__main__")
    calls = []

    def fake_main():
        calls.append("app")
        return 0

    module.main = fake_main
    monkeypatch.setitem(sys.modules, "omicsone_streamlit.__main__", module)

    assert commands.main(["app"]) == 0
    assert calls == ["app"]


