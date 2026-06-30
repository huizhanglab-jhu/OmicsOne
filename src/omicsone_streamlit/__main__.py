def main():
    import subprocess
    import sys
    from pathlib import Path

    app_path = Path(__file__).parent / "app.py"
    if not app_path.exists():
        raise SystemExit("Cannot find app.py inside package.")

    cmd = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(app_path),
        "--server.port",
        "8501",
        "--server.headless",
        "true",
        "--browser.gatherUsageStats",
        "false",
    ]
    cmd.extend(sys.argv[1:])
    raise SystemExit(subprocess.call(cmd))


if __name__ == "__main__":
    main()
