from __future__ import annotations

from pathlib import Path


def write_rust_spearman_input(
    input_file: str | Path,
    output_file: str | Path,
    *,
    has_header: bool = True,
    has_index: bool = True,
    sep: str = "\t",
) -> Path:
    """Write a numeric-only matrix file suitable for Rust Spearman.

    The canonical OmicsOne intermediate matrix keeps sample headers and a gene
    index. The Rust implementation expects only numeric matrix values, so this
    adapter strips the optional first row and first column without loading the
    full matrix into memory.
    """

    input_path = Path(input_file)
    output_path = Path(output_file)

    if not input_path.is_file():
        raise FileNotFoundError(f"Input matrix does not exist: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", encoding="utf-8", newline="") as source:
        with output_path.open("w", encoding="utf-8", newline="\n") as target:
            if has_header:
                next(source, None)

            for line in source:
                stripped = line.rstrip("\r\n")
                if not stripped:
                    target.write("\n")
                    continue

                values = stripped.split(sep)
                if has_index:
                    values = values[1:]

                target.write(sep.join(values))
                target.write("\n")

    return output_path
