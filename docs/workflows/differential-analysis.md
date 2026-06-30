# Differential Analysis

Differential analysis compares two matrices, computes feature-level statistical
tests, writes volcano plots, maps gene IDs, and optionally runs pathway
enrichment.

## Implementation

- Service: `src/omicsone/services/volcano_enrichment.py`
- replay adapter: `src/omicsone/replay/differential.py`
- Config template: `src/omicsone/replay/default_differential_config.ini`
- FastAPI router: `src/omicsone/api/routers/volcano.py`
- Schema: `src/omicsone/api/schemas/volcano.py`
- Helper script: `tools/run_volcano_enrichment_config.py`

## Input Matrices

Inputs are TSV matrices:

- rows are features.
- columns are samples.
- first column is used as the index unless an `idx` column exists.
- values should already be log2-like.

## Methods

Supported methods:

- `Wilcoxon(Unpaired)`
- `Wilcoxon(Paired)`
- `T-test(Unpaired)`
- `T-test(Paired)`

Paired methods require the two groups to have equal sample counts in matched
order.

## Python Example

```python
from omicsone.services.volcano_enrichment import generate_volcano_enrichment

result = generate_volcano_enrichment(
    normal_path="normal.tsv",
    tumor_path="tumor.tsv",
    output_dir="out/diff",
    method="Wilcoxon(Unpaired)",
    fdr_cutoff=0.01,
    log2fc_cutoff=1.0,
)
```

## Config Example

```powershell
python tools/run_volcano_enrichment_config.py --config path\to\config.ini
```

## API Example

```json
{
  "normal_path": "normal.tsv",
  "tumor_path": "tumor.tsv",
  "output_dir": "out/diff",
  "method": "Wilcoxon(Unpaired)",
  "fdr_cutoff": 0.01,
  "log2fc_cutoff": 1.0
}
```

Post to:

```text
POST /api/v1/diff/volcano/enrichment
```

## Outputs

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
- optional n8n replay JavaScript when `write_n8n_script = true`

## Interpretation

Log2FC is tumor minus normal.

Significance labels:

- `S-U`: FDR-significant and Log2FC above the positive cutoff.
- `S-D`: FDR-significant and Log2FC below the negative cutoff.
- `U`: FDR-significant and positive, but not beyond the strong Log2FC cutoff.
- `D`: FDR-significant and negative, but not beyond the strong Log2FC cutoff.

