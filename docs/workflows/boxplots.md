# Boxplots

The boxplot workflow creates normal-vs-tumor PDF boxplots for selected genes.

## Implementation

- Service: `src/omicsone/services/boxplot_figures.py`
- API router: `src/omicsone/api/routers/boxplots.py`
- Schema: `src/omicsone/api/schemas/boxplots.py`
- replay adapter: `src/omicsone/replay/boxplots.py`
- Config template: `src/omicsone/replay/default_boxplot_config.ini`

## API

```text
POST /api/v1/diff/boxplot/figures
```

## Inputs

- normal/NAT matrix TSV.
- tumor matrix TSV.
- output directory.
- FASTA path for gene ID mapping.
- gene symbols or Ensembl IDs.
- plot styling parameters.

## Python Example

```python
from omicsone.services.boxplot_figures import generate_boxplot_figures

result = generate_boxplot_figures(
    normal_path="normal.tsv",
    tumor_path="tumor.tsv",
    output_dir="out/boxplots",
    genes=["TP53", "EGFR"],
)
```

## Replay Example

```python
from omicsone.replay.boxplots import run_boxplot_figures

run_boxplot_figures("input.ini")
```

## Outputs

- one boxplot PDF per generated gene.
- `summary.tsv`.
- `result.log`.
- optional n8n replay JavaScript when `write_n8n_script = true`.
- generated count.
- missing gene list.

