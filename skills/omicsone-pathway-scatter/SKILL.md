---
name: omicsone-pathway-scatter
description: >-
  Run, debug, document, or extend OmicsOne protein-vs-phosphosite pathway
  scatter workflows, including pathway GMT loading, highlight phosphosite
  tables, replay config parsing, generated pathway scatter plots, and direct
  service/config execution.
---

# OmicsOne Pathway Scatter

## Overview

Primary files:

- `src/omicsone/services/pathway_scatter.py`
- `src/omicsone/replay/pathway_scatter.py`
- `src/omicsone/replay/default_pathway_scatter_config.ini`
- `tools/run_pathway_scatter_config.py`

## Quick Start

Direct config helper:

```powershell
python tools/run_pathway_scatter_config.py --config path\to\config.ini
```

replay/config:

```python
from omicsone.replay.pathway_scatter import run_pathway_scatter_plots

run_pathway_scatter_plots("input.ini")
```

Service:

```python
from omicsone.services.pathway_scatter import run_pathway_scatter_analysis

results = run_pathway_scatter_analysis(
    protein_diff_path="protein_diff.tsv",
    phospho_diff_path="phospho_diff.tsv",
    output_dir="out/pathway_scatter",
    pathways={"MY_PATHWAY": "HALLMARK_MYC_TARGETS_V1"},
    highlight_paths={"MY_PATHWAY": "highlights.tsv"},
)
```

## Inputs

- protein differential result table.
- phosphosite differential result table.
- FASTA path and optional GMT path.
- pathway definitions from config sections or a pathway/highlight table.
- highlight files for pathway-specific phosphosites.

## Outputs

Per pathway:

- points TSV.
- highlight TSV.
- missing highlights TSV.
- highlight text file.
- PNG/PDF/TIFF scatter plots.

Batch workflows also write `pathway_scatter_summary.tsv`.

## Core Rules

- Use `run_pathway_scatter_plots` for replay-friendly input.ini workflows.
- Use `run_pathway_scatter_analysis_from_config` only for direct service-style
  configs.
- Ensure pathway folder names are unique.
- Check highlight file paths before assuming missing phosphosites are biological
  rather than path/config errors.

