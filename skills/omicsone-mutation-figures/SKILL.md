---
name: omicsone-mutation-figures
description: >-
  Run, debug, document, or extend OmicsOne mutation figure workflows, including
  mutation heatmaps, mutation-type distribution PDFs, MAF filtering, binary
  mutation matrix workflows, v1/v2 mutation FastAPI endpoints, and replay
  mutation config execution.
---

# OmicsOne Mutation Figures

## Overview

Primary files:

- `src/omicsone/services/mutation_figures.py`
- `src/omicsone/api/routers/mutations.py`
- `src/omicsone/api/routers/mutations_v2.py`
- `src/omicsone/api/schemas/mutations.py`
- `src/omicsone/replay/mutations.py`
- `src/omicsone/replay/default_mutation_config.ini`

## API

V1:

- `POST /api/v1/mutations/heatmap/figures`

V2:

- `POST /api/v2/mutations/heatmap/figures`

The legacy `/api/v1/mutations/hnsc/figures` route exists but is hidden from the
OpenAPI schema.

## Inputs

V1 uses:

- mutation Excel file.
- MAF file.
- output directory.
- mutation threshold.

V2 uses:

- binary mutation matrix.
- metadata table.
- MAF file.
- output directory.
- optional gene-symbol map path or species-specific packaged map.

## Replay Workflow

Use:

```python
from omicsone.replay.mutations import run_mutation_figures

run_mutation_figures("input.ini")
```

To call a running API from config:

```python
from omicsone.replay.mutations import post_mutation_figures_api

post_mutation_figures_api("input.ini")
```

## Outputs

- heatmap PDF.
- mutation-type distribution PDF.
- `result.log`.
- summary counts including gene count, sample count, total MAF rows, filtered
  MAF rows, and found mutation types.

## Core Rules

- Use V2 for binary mutation matrix plus metadata workflows.
- Keep MAF filtering behavior in `mutation_figures.py`.
- Do not add cohort-specific defaults to routers; keep defaults in schemas or
  services.
- Preserve editable vector PDF output when modifying plots.

