---
name: omicsone-differential-analysis
description: >-
  Run, debug, document, or extend OmicsOne differential analysis workflows,
  including two-group Wilcoxon/t-test comparisons, volcano plots, gene mapping,
  Enrichr/GSEApy pathway enrichment, config.ini/replay execution, FastAPI
  endpoints, and Streamlit differential-analysis UI behavior.
---

# OmicsOne Differential Analysis

## Overview

Use the existing OmicsOne differential-analysis engine. Do not reimplement the
statistics unless the user explicitly asks to change the engine.

Primary files:

- `src/omicsone/services/volcano_enrichment.py`
- `src/omicsone/replay/differential.py`
- `src/omicsone/replay/default_differential_config.ini`
- `src/omicsone/api/routers/volcano.py`
- `src/omicsone/api/schemas/volcano.py`
- `tools/run_volcano_enrichment_config.py`

Legacy Streamlit UI files:

- `src/omicsone_streamlit/utils/omicsone_diff.py`
- `src/omicsone_streamlit/utils/omicsone_volcano.py`
- `src/omicsone_streamlit/mypages/analysis_diff.py`

## Inputs

Expect two tab-separated matrices:

- `normal_path`: normal/NAT matrix.
- `tumor_path`: tumor matrix.
- Rows are features and columns are samples.
- The first column is used as the feature index, or `idx` when present.
- Values are expected to already be log2-like, because Log2FC is computed by
  subtraction.

For RNA/protein gene-level matrices, usually keep `strip_feature_version=true`.
For phosphosite IDs such as `ENSG|ENSP|S123`, set `strip_feature_version=false`.

## Quick Start

Config-file workflow:

```powershell
python tools/run_volcano_enrichment_config.py --config path\to\config.ini
```

Python workflow:

```python
from omicsone.services.volcano_enrichment import generate_volcano_enrichment

result = generate_volcano_enrichment(
    normal_path="normal.tsv",
    tumor_path="tumor.tsv",
    output_dir="out/diff",
)
```

replay-style workflow:

```python
from omicsone.replay.differential import run_differential_analysis

run_differential_analysis("input.ini")
```

Phosphosite workflow:

```python
from omicsone.replay.differential import run_phospho_differential_analysis

run_phospho_differential_analysis("input.ini")
```

## API

- `POST /api/v1/diff/volcano/enrichment`
- `POST /api/v1/diff/volcano/enrichment/preset`

The preset endpoint resolves presets from `DEFAULT_PRESETS` in
`volcano_enrichment.py`, then applies request overrides.

## Outputs

A successful run writes:

- `diff.tsv`
- `combined_matrix.tsv`
- `pure_up_genes.tsv`
- `pure_down_genes.tsv`
- `total_genes.tsv`
- `up_enrichr_df.tsv`
- `down_enrichr_df.tsv`
- `enrichment_plot_table.tsv`
- volcano PNG/PDF/TIFF
- enrichment PNG/PDF/TIFF
- HTML report
- `result.log`
- n8n replay JavaScript

## Core Rules

- Keep tumor-vs-normal direction clear: Log2FC is tumor minus normal.
- Paired methods require equal sample counts and matched sample order.
- `S-U` and `S-D` pass both FDR and absolute Log2FC thresholds.
- `U` and `D` pass FDR but not the strong absolute Log2FC threshold.
- Do not enable remote Enrichr unless the user accepts sending gene lists to a
  remote service.
- Prefer local GMT files or cached GSEApy libraries for reproducible enrichment.
- If results are unexpectedly empty, inspect missingness thresholds,
  `min_sample_size`, sample overlap, and feature ID stripping before changing
  statistical code.

## References

Read `references/config-and-outputs.md` before editing default configs,
documenting outputs, or diagnosing config-driven runs.

