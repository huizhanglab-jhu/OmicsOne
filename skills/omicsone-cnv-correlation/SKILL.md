---
name: omicsone-cnv-correlation
description: >-
  Run, debug, document, or extend OmicsOne CNV correlation workflows, including
  Rust-backed Spearman row-pair correlations, paired CNV/RNA/protein matrix
  preprocessing, CNV correlation figures, GISTIC panels, FastAPI endpoints,
  and replay/config-driven CNV correlation pipeline execution.
---

# OmicsOne CNV Correlation

## Overview

Use the service and utility layers rather than rebuilding correlation logic.

Primary files:

- `src/omicsone/services/cnv_correlation_pipeline.py`
- `src/omicsone/services/spearman_omics.py`
- `src/omicsone/services/spearman_preprocessing.py`
- `src/omicsone/utils/spearmanr.py`
- `packages/rust_spearmanr/src/lib.rs`
- `src/omicsone_streamlit/plots/cnv_correlation.py`
- `src/omicsone/api/routers/cnv_correlation.py`
- `src/omicsone/api/routers/spearman.py`
- `src/omicsone/replay/cnv_correlation_pipeline.py`

## Quick Start

Use the installed replay command when available:

```powershell
omicsone-replay-cnv-correlation --help
```

Use the Python service for direct computation:

```python
from pathlib import Path
from omicsone.services.spearman_omics import OmicsInput, compute_paired_omics_spearman

result = compute_paired_omics_spearman(
    inputs=[
        OmicsInput(path=Path("cnv.tsv"), data_type="cnv"),
        OmicsInput(path=Path("protein.tsv"), data_type="protein"),
    ],
    output_dir=Path("out/spearman"),
    min_valid_pairs=4,
)
```

## API

Spearman:

- `POST /api/v1/spearman/compute-file`
- `POST /api/v1/spearman/paired-omics`

CNV correlation:

- `POST /api/v1/cnv-correlation/pipeline`
- `GET /api/v1/cnv-correlation/pipeline/{job_id}`
- `POST /api/v1/cnv-correlation/figures`

## Rust Backend

The Rust source is `packages/rust_spearmanr/src/lib.rs`. It exposes:

- `spearman(x, y, min_valid_pairs=2)`
- `compute_file(input_file1, input_file2, output_file, min_valid_pairs=2)`

Use it from Python through `omicsone.utils.spearmanr`. The wrapper falls back to
Python if the Rust extension is unavailable.

## Inputs

Raw omics matrix inputs should contain an `idx` column plus sample columns.
Preprocessed Rust Spearman inputs are whitespace-delimited numeric-only matrix
files where missing values are represented as `NaN`.

CNV figure inputs include:

- Spearman correlation output.
- CNV `*_for_corr.txt` file.
- RNA or protein `*_for_corr.txt` file.
- GENCODE FASTA file.
- chromosome length XLSX file.
- cytoband TSV file.
- optional GISTIC file.

## Core Rules

- Do not call the Rust extension directly from workflow code; use
  `omicsone.utils.spearmanr`.
- Keep CNV/RNA/protein sample and gene intersections explicit.
- Use `min_valid_pairs` to control missing-data tolerance.
- Large Spearman outputs should stay on disk; read them in chunks for plotting.
- For figure changes, preserve PNG, editable PDF, and TIFF outputs.
- For API pipeline runs, use the status endpoint to inspect async progress.

