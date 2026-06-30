# CNV Correlation

The CNV correlation workflow computes pairwise Spearman correlations between
CNV and RNA/protein matrices, then generates genome-positioned correlation
figures.

## Implementation

- Pipeline service: `src/omicsone/services/cnv_correlation_pipeline.py`
- Spearman service: `src/omicsone/services/spearman_omics.py`
- Spearman preprocessing: `src/omicsone/services/spearman_preprocessing.py`
- Spearman wrapper: `src/omicsone/utils/spearmanr.py`
- Rust backend: `packages/rust_spearmanr/src/lib.rs`
- Figure generation: `src/omicsone/plots/cnv_correlation.py`
- API router: `src/omicsone/api/routers/cnv_correlation.py`
- Spearman API router: `src/omicsone/api/routers/spearman.py`
- replay adapter: `src/omicsone/replay/cnv_correlation_pipeline.py`

## Rust Spearman

The Rust extension exposes:

- `spearman`
- `compute_file`

Use it through:

```python
from omicsone.utils import spearmanr
```

The wrapper reports which backend is active:

```python
spearmanr.backend()
```

## API Pipeline

Start:

```text
POST /api/v1/cnv-correlation/pipeline
```

Check status:

```text
GET /api/v1/cnv-correlation/pipeline/{job_id}
```

At least one of `rna_path` or `protein_path` is required.

## Figure Generation

Use:

```text
POST /api/v1/cnv-correlation/figures
```

Required inputs:

- correlation file.
- CNV `*_for_corr.txt` file.
- target RNA/protein `*_for_corr.txt` file.
- FASTA file.
- chromosome length XLSX file.
- cytoband TSV file.
- output directory.

Optional:

- GISTIC file.
- target type.
- correlation threshold.
- chunksize.
- DPI.

## Outputs

The pipeline writes intermediate files, Spearman outputs, figure files, and
JSON/HTML reports. Figure generation writes:

- annotated correlations table.
- CNV distribution counts.
- optional GISTIC counts.
- panel PNG/PDF/TIFF files.
- combined figure PNG/PDF/TIFF.

