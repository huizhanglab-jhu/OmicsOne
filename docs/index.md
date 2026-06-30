# OmicsOne Documentation

OmicsOne is a Python package and Streamlit application for multi-omics analysis
workflows. The package exposes reusable Python services, a FastAPI application,
replay-friendly config adapters, and helper scripts for batch execution.

## Main Runtime Surfaces

- Python package: `omicsone`
- Streamlit command: `omicsone`
- FastAPI command: `omicsone-api`
- Replay CNV pipeline command: `omicsone-replay-cnv-correlation`

## Source Layout

```text
src/omicsone/
  api/        FastAPI app, routers, request/response schemas
  services/   reusable analysis engines
  replay/      config.ini and replay-friendly adapters
  resources/  packaged reference resources
  utils/      shared low-level utilities

src/omicsone_streamlit/
  app.py      Streamlit application entry point
  mypages/    Streamlit pages
  plots/      plotting helpers used by the UI and services
  utils/      older Streamlit-linked analysis utilities

packages/rust_spearmanr/
  Rust/PyO3 Spearman correlation backend

tools/
  repository helper scripts
```

## Recommended Use

- Use `src/omicsone/services` for direct Python workflows and reusable analysis.
- Use `src/omicsone/api` for HTTP service integration.
- Use `src/omicsone/replay` for config-file workflows and local workflow callers.
- Use `src/omicsone_streamlit` for interactive UI changes.
- Use `skills/` for AI-agent instructions and `docs/` for human documentation.

## Documentation Map

- [Installation](installation.md)
- [API Overview](api/overview.md)
- [FastAPI Endpoints](api/fastapi-endpoints.md)
- [Python API](api/python-api.md)
- [CLI And Tools](tools/cli-tools.md)
- [Config Files](tools/config-files.md)
- [Differential Analysis](workflows/differential-analysis.md)
- [CNV Correlation](workflows/cnv-correlation.md)
- [Mutation Figures](workflows/mutation-figures.md)
- [Boxplots](workflows/boxplots.md)
- [Pathway Scatter](workflows/pathway-scatter.md)
- [Agent Skills](skills.md)

