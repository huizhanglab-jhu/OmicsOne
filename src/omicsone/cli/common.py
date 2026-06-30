from __future__ import annotations

import argparse


def add_config_run_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", required=True, help="Path to the input config.ini.")
    parser.add_argument("--output-dir", help="Override output_dir from the config.")
    parser.add_argument("--quiet", action="store_true", help="Do not print JSON output.")


def print_json_enabled(args: argparse.Namespace) -> bool:
    return not bool(getattr(args, "quiet", False))

