# Config Files

OmicsOne uses INI files for replay-friendly and batch workflows. Packaged default
templates live in:

```text
src/omicsone/replay/
```

## Templates

| Workflow | Template |
| --- | --- |
| Differential volcano/enrichment | `default_differential_config.ini` |
| Boxplots | `default_boxplot_config.ini` |
| Mutation figures | `default_mutation_config.ini` |
| CNV correlation pipeline | `default_cnv_correlation_pipeline_config.ini` |
| Pathway scatter | `default_pathway_scatter_config.ini` |

## Adapter Pattern

Most replay adapters follow the same pattern:

1. Read `input.ini`.
2. Resolve the workflow output directory.
3. Copy the packaged default config into the output directory if missing.
4. Load the editable output config.
5. Merge default config, input settings, and explicit overrides.
6. Run the service.
7. Print or return a JSON-compatible result dictionary.

## Common Input Keys

Common path keys include:

- `normal_path`
- `tumor_path`
- `fasta_path`
- `output_dir`
- `out_dir`
- `cohort`
- `omics`

Workflow-specific keys are documented in the workflow pages.

## Release Note

The default INI templates are included in package data through
`pyproject.toml`:

```toml
[tool.setuptools.package-data]
"omicsone" = [
    "resources/**/*.tsv",
    "replay/*.ini",
]
```

