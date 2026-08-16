from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from omicsone.omicsx.pipeline import run_omicsx_config


def run_omicsx_analysis(
    input_ini: str | Path,
    *,
    output_dir: str | Path | None = None,
    print_json: bool = True,
) -> dict[str, Any]:
    """Replay the reusable OmicsX package from an INI file."""
    return run_omicsx_config(input_ini, output_dir=output_dir, print_json=print_json)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run pairwise OmicsX analysis from a replay config.")
    parser.add_argument("--config", required=True, help="Path to the OmicsX config.ini.")
    parser.add_argument("--output-dir", help="Override output_dir from the config.")
    parser.add_argument("--quiet", action="store_true", help="Do not print JSON output.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    run_omicsx_analysis(args.config, output_dir=args.output_dir, print_json=not args.quiet)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
