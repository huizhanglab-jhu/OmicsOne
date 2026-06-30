from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def _python_files(path: Path) -> list[Path]:
    return sorted(file for file in path.rglob("*.py") if file.is_file())


def _imported_modules(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            modules.add(node.module)
    return modules


def test_services_do_not_depend_on_entry_layers():
    banned_exact = {"fastapi", "streamlit"}
    banned_prefixes = ("omicsone_streamlit", "omicsone.cli")
    violations = []

    for path in _python_files(SRC / "omicsone" / "services"):
        for module in _imported_modules(path):
            if module in banned_exact or module.startswith(banned_prefixes):
                violations.append(f"{path.relative_to(ROOT)} imports {module}")

    assert violations == []


def test_api_routers_do_not_import_streamlit_modules():
    violations = []

    for path in _python_files(SRC / "omicsone" / "api" / "routers"):
        for module in _imported_modules(path):
            if module.startswith("omicsone_streamlit"):
                violations.append(f"{path.relative_to(ROOT)} imports {module}")

    assert violations == []


def test_cli_does_not_import_services_from_services_layer():
    violations = []

    for path in _python_files(SRC / "omicsone" / "services"):
        for module in _imported_modules(path):
            if module.startswith("omicsone.cli"):
                violations.append(f"{path.relative_to(ROOT)} imports {module}")

    assert violations == []
