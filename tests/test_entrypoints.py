from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def _run_module(*args: str) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC) + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, "-m", *args],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def test_omicsone_module_help():
    result = _run_module("omicsone", "--help")

    assert result.returncode == 0
    assert "Run OmicsOne local workflows" in result.stdout
    assert "cnv-correlation" in result.stdout


def test_pyproject_scripts_define_four_entry_layers():
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")

    assert 'omicsone = "omicsone.cli.__main__:main"' in pyproject
    assert 'omicsone-api = "omicsone.api.__main__:main"' in pyproject
    assert 'omicsone-app = "omicsone_streamlit.__main__:main"' in pyproject
    assert 'omicsone-replay-cnv-correlation = "omicsone.replay.cnv_correlation_pipeline:main"' in pyproject


def test_omicsone_differential_help():
    result = _run_module("omicsone", "differential", "run", "--help")

    assert result.returncode == 0
    assert "--config" in result.stdout
    assert "--quiet" in result.stdout


def test_omicsone_cnv_correlation_help():
    result = _run_module("omicsone", "cnv-correlation", "run", "--help")

    assert result.returncode == 0
    assert "--config" in result.stdout
    assert "--output-dir" in result.stdout


def test_api_main_is_importable():
    from omicsone.api.__main__ import main

    assert callable(main)


def test_streamlit_main_is_importable():
    from omicsone_streamlit.__main__ import main

    assert callable(main)


def test_compatibility_imports():
    from omicsone.replay import run_cnv_correlation_pipeline
    from omicsone_streamlit.plots.cnv_correlation import generate_cnv_correlation_figures
    from omicsone_streamlit.utils.fasta import get_gene_map

    assert callable(run_cnv_correlation_pipeline)
    assert callable(generate_cnv_correlation_figures)
    assert callable(get_gene_map)

