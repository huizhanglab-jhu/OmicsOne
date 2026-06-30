# Python API

Prefer the service layer for direct Python workflows.

The stable Python entrypoint is `omicsone.services`. Services contain the
business logic and do not require CLI, FastAPI, or Streamlit entry layers.

## Differential Analysis

```python
from omicsone.services.volcano_enrichment import generate_volcano_enrichment

result = generate_volcano_enrichment(
    normal_path="normal.tsv",
    tumor_path="tumor.tsv",
    output_dir="out/diff",
)
```

Config-driven:

```python
from omicsone.replay.differential import run_differential_analysis

run_differential_analysis("input.ini")
```

Phosphosite config-driven:

```python
from omicsone.replay.differential import run_phospho_differential_analysis

run_phospho_differential_analysis("input.ini")
```

## Boxplots

```python
from omicsone.services.boxplot_figures import generate_boxplot_figures

result = generate_boxplot_figures(
    normal_path="normal.tsv",
    tumor_path="tumor.tsv",
    output_dir="out/boxplots",
    genes=["TP53", "EGFR"],
)
```

Config-driven:

```python
from omicsone.replay.boxplots import run_boxplot_figures

run_boxplot_figures("input.ini")
```

## Spearman Correlation

```python
from omicsone.utils import spearmanr

r = spearmanr.spearman([1, 2, 3], [3, 2, 1], min_valid_pairs=2)
rows = spearmanr.compute_file("a.txt", "b.txt", "out.txt", min_valid_pairs=4)
backend = spearmanr.backend()
```

## Paired Omics Spearman

```python
from pathlib import Path
from omicsone.services.spearman_omics import OmicsInput, compute_paired_omics_spearman

result = compute_paired_omics_spearman(
    inputs=[
        OmicsInput(path=Path("cnv.tsv"), data_type="cnv"),
        OmicsInput(path=Path("rna.tsv"), data_type="rna"),
        OmicsInput(path=Path("protein.tsv"), data_type="protein"),
    ],
    output_dir=Path("out/spearman"),
    min_valid_pairs=4,
)
```

## Mutation Figures

```python
from omicsone.services.mutation_figures import generate_hnsc_mutation_figures

result = generate_hnsc_mutation_figures(
    mutation_excel_path="mutation.xlsx",
    maf_path="variants.maf",
    output_dir="out/mutations",
)
```

V2 binary-matrix workflow:

```python
from omicsone.services.mutation_figures import generate_mutation_figures_from_binary

result = generate_mutation_figures_from_binary(
    mutation_binary_path="mutation_binary.tsv",
    meta_path="meta.tsv",
    maf_path="variants.maf",
    output_dir="out/mutations_v2",
)
```

## Pathway Scatter

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

