"""Spearman correlation utilities.

The public API is intentionally small:

    from omicsone.utils import spearmanr
    spearmanr.compute_file("cnv.txt", "rna.txt", "out.txt", min_valid_pairs=4)

When the Rust extension is installed, calls are dispatched to it. A pure Python
fallback keeps the module usable in source checkouts before the extension is
built.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Optional, Sequence, Union

try:
    from . import _spearmanr as _rust
except ImportError:  # pragma: no cover - exercised when Rust extension is absent
    _rust = None


Number = Optional[Union[float, int]]


def backend() -> str:
    """Return the active implementation backend."""

    return "rust" if _rust is not None else "python"


def spearman(
    x: Sequence[Number],
    y: Sequence[Number],
    min_valid_pairs: int = 2,
) -> float:
    """Compute Spearman correlation for two vectors.

    ``None`` and ``NaN`` are treated as missing values and removed pairwise.
    ``NaN`` is returned when fewer than ``min_valid_pairs`` valid pairs remain
    or when either ranked vector has zero variance.
    """

    if _rust is not None:
        return float(_rust.spearman(x, y, min_valid_pairs))

    return _spearman_python(x, y, min_valid_pairs)


def compute_file(
    input_file1: str | Path,
    input_file2: str | Path,
    output_file: str | Path,
    min_valid_pairs: int = 2,
) -> int:
    """Compute all row-pair Spearman correlations between two matrix files.

    Input files are whitespace-delimited numeric matrices. Tokens equal to
    ``NaN`` are treated as missing values. The output format matches the C#
    tool: ``row_index_1<TAB>row_index_2<TAB>correlation``.
    """

    if _rust is not None:
        return int(
            _rust.compute_file(
                str(input_file1),
                str(input_file2),
                str(output_file),
                min_valid_pairs,
            )
        )

    matrix1 = _read_matrix(input_file1)
    matrix2 = _read_matrix(input_file2)

    if matrix1 and matrix2 and len(matrix1[0]) != len(matrix2[0]):
        raise ValueError("input matrices must have the same number of columns")

    rows: list[str] = []
    for i, row1 in enumerate(matrix1):
        for j, row2 in enumerate(matrix2):
            corr = _spearman_python(row1, row2, min_valid_pairs)
            rows.append(f"{i}\t{j}\t{_format_float(corr)}")

    Path(output_file).write_text("\n".join(rows) + ("\n" if rows else ""))
    return len(rows)


def _read_matrix(path: str | Path) -> list[list[float]]:
    matrix: list[list[float]] = []
    expected_cols: int | None = None

    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            values = [_parse_token(token) for token in line.split()]
            if expected_cols is None:
                expected_cols = len(values)
            elif len(values) != expected_cols:
                raise ValueError(
                    f"{path} has {len(values)} columns on line {line_number}; "
                    f"expected {expected_cols}"
                )
            matrix.append(values)

    return matrix


def _parse_token(token: str) -> float:
    if token.lower() == "nan":
        return math.nan
    try:
        return float(token)
    except ValueError:
        return math.nan


def _spearman_python(
    x: Sequence[Number],
    y: Sequence[Number],
    min_valid_pairs: int,
) -> float:
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if min_valid_pairs < 1:
        raise ValueError("min_valid_pairs must be at least 1")

    valid_x: list[float] = []
    valid_y: list[float] = []

    for a, b in zip(x, y):
        fa = _as_float(a)
        fb = _as_float(b)
        if fa is None or fb is None:
            continue
        valid_x.append(fa)
        valid_y.append(fb)

    if len(valid_x) < min_valid_pairs:
        return math.nan

    return _pearson(_rank(valid_x), _rank(valid_y))


def _as_float(value: Number) -> float | None:
    if value is None:
        return None
    result = float(value)
    if math.isnan(result):
        return None
    return result


def _rank(values: Sequence[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    i = 0

    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1

        rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[indexed[k][0]] = rank
        i = j

    return ranks


def _pearson(x: Sequence[float], y: Sequence[float]) -> float:
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)

    cov = 0.0
    var_x = 0.0
    var_y = 0.0
    for a, b in zip(x, y):
        da = a - mean_x
        db = b - mean_y
        cov += da * db
        var_x += da * da
        var_y += db * db

    denom = math.sqrt(var_x * var_y)
    if denom == 0.0:
        return math.nan
    return cov / denom


def _format_float(value: float) -> str:
    if math.isnan(value):
        return "NaN"
    return f"{value:.17g}"
