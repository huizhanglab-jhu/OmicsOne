---
name: omicsone
description: >-
  Orient agents in the OmicsOne repository and choose the correct OmicsOne
  runtime surface. Use when a task involves using, debugging, documenting, or
  extending the omicsone Python package, Streamlit app, FastAPI service, replay
  adapters, analysis services, or bundled command-line tools.
---

# OmicsOne

## Overview

Use this skill first when working in the OmicsOne repository. Choose the
smallest runtime surface that matches the user request.

## Repository Map

- `src/omicsone/services/`: reusable analysis engines. Prefer this layer for
  new programmatic workflows.
- `src/omicsone/api/`: FastAPI app, routers, and Pydantic request/response
  schemas.
- `src/omicsone/replay/`: config-driven adapters for local replay and batch execution.
- `src/omicsone/utils/`: reusable low-level utilities, including the Spearman
  Python wrapper for the Rust backend.
- `src/omicsone/resources/`: packaged static resources.
- `src/omicsone_streamlit/`: Streamlit app, UI pages, plotting code, and older
  UI-linked utilities.
- `packages/rust_spearmanr/`: Rust/PyO3 Spearman correlation backend.
- `tools/`: developer/helper scripts that are not all installed as package
  commands.

## Runtime Surfaces

- Use Python services when the user wants direct computation or code changes.
- Use FastAPI routers when the user asks about HTTP endpoints or service
  integration.
- Use replay adapters when the user asks for config.ini workflows or batch runs.
- Use Streamlit code only for UI behavior.
- Use Rust Spearman only through `omicsone.utils.spearmanr` unless editing the
  native backend itself.

## Commands

Run the Streamlit app:

```powershell
omicsone
```

Run the API:

```powershell
omicsone-api
```

Run installed Replay CNV correlation pipeline:

```powershell
omicsone-replay-cnv-correlation --help
```

Run helper scripts from the repository root:

```powershell
python tools/run_volcano_enrichment_config.py --config path\to\config.ini
python tools/run_pathway_scatter_config.py --config path\to\config.ini
```

## Core Rules

- Prefer `src/omicsone/services` for reusable computation.
- Keep UI-only code in `src/omicsone_streamlit`.
- Do not duplicate statistical logic in API or replay wrappers; call services.
- Treat `docs/` as human documentation and `skills/` as agent instructions.
- Before changing public behavior, update the matching docs and skills.
- Avoid committing generated artifacts: `build/`, `dist/`, `__pycache__/`,
  local logs, and temporary output folders.

## Related Skills

- `omicsone-differential-analysis`: volcano, differential testing, enrichment.
- `omicsone-cnv-correlation`: CNV/RNA/protein correlation pipeline and figures.
- `omicsone-mutation-figures`: mutation heatmaps and mutation-type plots.
- `omicsone-boxplots`: normal-vs-tumor gene boxplots.
- `omicsone-pathway-scatter`: protein/phosphosite pathway scatter plots.

