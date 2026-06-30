# Pathway Scatter

Pathway scatter plots compare protein and phosphosite differential results for
selected pathways and highlighted sites.

## Implementation

- Service: `src/omicsone/services/pathway_scatter.py`
- replay adapter: `src/omicsone/replay/pathway_scatter.py`
- Config template: `src/omicsone/replay/default_pathway_scatter_config.ini`
- Helper script: `tools/run_pathway_scatter_config.py`

## Inputs

- protein differential result table.
- phosphosite differential result table.
- output directory.
- pathway definitions.
- highlight phosphosite files.
- optional FASTA path.
- optional GMT path.

Pathways can be supplied through config sections or by a pathway/highlight TSV.
Recognized table columns include:

- `pathway_term`
- `highlight_path`
- optional folder column such as `pathway_folder`, `folder_name`, `folder`, or
  `name`.

## Python Example

```python
from omicsone.services.pathway_scatter import run_pathway_scatter_analysis

results = run_pathway_scatter_analysis(
    protein_diff_path="protein_diff.tsv",
    phospho_diff_path="phospho_diff.tsv",
    output_dir="out/pathway_scatter",
    pathways={"my_pathway": "HALLMARK_MYC_TARGETS_V1"},
    highlight_paths={"my_pathway": "highlights.tsv"},
)
```

## Config Example

```powershell
python tools/run_pathway_scatter_config.py --config path\to\config.ini
```

## Replay Example

```python
from omicsone.replay.pathway_scatter import run_pathway_scatter_plots

run_pathway_scatter_plots("input.ini")
```

## Phosphosite-Protein Table Pipeline

Use this workflow when the pathway scatter plots should be generated from
phosphosite-protein tables first, with a configurable gene selection rule.

```powershell
python tools/run_phosphosite_protein_pathway_pipeline.py --config path\to\config.ini
```

Set `[selection] gene_selection_mode` in the config:

- `filtered`: E2F/G2M use `pure_up_genes.tsv`; myogenesis uses
  `pure_down_genes.tsv`, unless overridden in `[directions]`.
- `all_pathway`: use all genes in the GMT pathway that are present in the
  phospho differential table.

Required extra inputs for `filtered` mode:

- `phospho_run_dir`, or
- explicit `pure_up_genes_path` and `pure_down_genes_path`.

## Outputs

Per pathway:

- points TSV.
- highlight TSV.
- missing highlights TSV.
- highlight text file.
- PNG/PDF/TIFF scatter plot.

Batch workflows write:

- `pathway_scatter_summary.tsv`.
- `phosphosite_protein_table_summary.tsv` for the phosphosite-protein table
  pipeline.
- `phosphosite_protein_pathway_tables.xlsx` for the phosphosite-protein table
  pipeline.

