---
name: omicsone-boxplots
description: >-
  Run, debug, document, or extend OmicsOne normal-vs-tumor boxplot generation,
  including gene ID mapping through FASTA, per-gene PDF output, summary TSVs,
  FastAPI boxplot endpoints, and replay/config-driven boxplot workflows.
---

# OmicsOne Boxplots

## Overview

Primary files:

- `src/omicsone/services/boxplot_figures.py`
- `src/omicsone/api/routers/boxplots.py`
- `src/omicsone/api/schemas/boxplots.py`
- `src/omicsone/replay/boxplots.py`
- `src/omicsone/replay/default_boxplot_config.ini`

## API

- `POST /api/v1/diff/boxplot/figures`

## Quick Start

Python:

```python
from omicsone.services.boxplot_figures import generate_boxplot_figures

result = generate_boxplot_figures(
    normal_path="normal.tsv",
    tumor_path="tumor.tsv",
    output_dir="out/boxplots",
    genes=["TP53", "EGFR"],
)
```

replay/config:

```python
from omicsone.replay.boxplots import run_boxplot_figures

run_boxplot_figures("input.ini")
```

## Inputs

- normal/NAT TSV matrix.
- tumor TSV matrix.
- FASTA path for Ensembl-to-symbol mapping.
- gene list as symbols or Ensembl IDs.
- output directory.

## Outputs

- one PDF per generated gene boxplot.
- `summary.tsv`.
- `result.log`.
- n8n replay JavaScript.
- missing gene list in the API/replay response.

## Core Rules

- Use the service layer for new workflow code.
- Use `genes_path` for large gene lists in config workflows.
- Report missing genes rather than silently ignoring them.
- Keep plot dimensions small and publication-oriented unless the user asks for
  a different layout.

