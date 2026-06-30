from pathlib import Path
def read_version():
    return Path(__file__).parents[1].joinpath("VERSION").read_text(encoding="utf-8").strip()
