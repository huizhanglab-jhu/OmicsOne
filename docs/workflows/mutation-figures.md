# Mutation Figures

Mutation workflows generate mutation heatmaps and mutation-type distribution
figures.

## Implementation

- Service: `src/omicsone/services/mutation_figures.py`
- API v1 router: `src/omicsone/api/routers/mutations.py`
- API v2 router: `src/omicsone/api/routers/mutations_v2.py`
- Schema: `src/omicsone/api/schemas/mutations.py`
- replay adapter: `src/omicsone/replay/mutations.py`
- Config template: `src/omicsone/replay/default_mutation_config.ini`

## API

V1 endpoint:

```text
POST /api/v1/mutations/heatmap/figures
```

V2 endpoint:

```text
POST /api/v2/mutations/heatmap/figures
```

## V1 Inputs

- mutation Excel path.
- MAF path.
- output directory.
- mutation threshold.
- cohort.
- optional output filenames or prefix.

## V2 Inputs

- binary mutation matrix.
- metadata table.
- MAF path.
- output directory.
- mutation threshold.
- cohort.
- species.
- optional gene symbol map path.
- optional output filenames or prefix.

## Python Example

```python
from omicsone.services.mutation_figures import generate_hnsc_mutation_figures

result = generate_hnsc_mutation_figures(
    mutation_excel_path="mutation.xlsx",
    maf_path="variants.maf",
    output_dir="out/mutations",
)
```

## Replay Example

```python
from omicsone.replay.mutations import run_mutation_figures

run_mutation_figures("input.ini")
```

## Outputs

- heatmap PDF.
- mutation-type distribution PDF.
- `result.log`.
- gene count.
- sample count.
- total and filtered MAF row counts.
- found mutation types.

