# OmicsOne

OmicsOne is a local multi-omics analysis toolkit for reproducible command-line
workflows. It focuses on turning proteomics, phosphoproteomics, RNA, CNV, and
mutation input tables into analysis outputs and publication-oriented figures.

Current status:

- CLI: available
- Streamlit app: coming soon
- FastAPI service: coming soon

## Available CLI Workflows

The command-line interface currently supports:

- Differential analysis for protein and phosphosite matrices
- Mutation figure generation
- Normal-vs-tumor protein boxplots
- Protein-vs-phosphosite pathway scatter plots
- CNV correlation workflows
- CNV correlation clean figure replay
- Rust/PyO3-backed Spearman correlation through `omicsone.utils.spearmanr`

## Install

Install from the repository root:

```powershell
cd <path-to-omicsone-streamlit>
python -m pip install -e .
```

The package metadata and pinned runtime dependencies are defined in
`pyproject.toml`.

After installation, these commands are available:

```powershell
omicsone
omicsone-replay-cnv-correlation
omicsone-replay-cnv-correlation-clean-figures
```

The Streamlit and FastAPI entrypoints may exist in the package, but they should
be treated as coming soon until their interfaces are finalized.

## CLI Help

Show top-level help:

```powershell
omicsone --help
```

Show workflow help:

```powershell
omicsone differential --help
omicsone mutations --help
omicsone boxplots --help
omicsone pathway-scatter --help
omicsone cnv-correlation --help
```

Show command-specific help:

```powershell
omicsone differential run --help
omicsone differential phospho --help
omicsone mutations run --help
omicsone boxplots run --help
omicsone pathway-scatter phosphosite-protein --help
omicsone cnv-correlation run --help
omicsone cnv-correlation clean-figures --help
```

Most replay commands support these common arguments:

```text
--config      Required. Path to an input/config INI file.
--output-dir  Optional. Override output_dir or out_dir from the config.
--quiet       Optional. Suppress JSON output on stdout.
```

## Replay Workflow

A replay workflow reruns an analysis from explicit input paths and settings.
The usual flow is:

1. Prepare an INI config with input file paths, output directories, and analysis
   settings.
2. Run the matching `omicsone <workflow> <action> --config <config.ini>`
   command.
3. OmicsOne writes outputs into the configured output directory.
4. Unless `--quiet` is used, the CLI prints a JSON summary containing output
   paths and analysis statistics.

Common commands:

```powershell
omicsone differential run --config configs\HNSCC_Protein.ini
omicsone differential phospho --config configs\HNSCC_phospho.ini
omicsone mutations run --config configs\HNSCC_Mutations.ini
omicsone boxplots run --config configs\HNSCC_Protein_boxplots.ini
omicsone pathway-scatter phosphosite-protein --config configs\HNSCC_pathways.ini
omicsone cnv-correlation run --config configs\HNCC_CNV.ini --output-dir runs\HNCC_CNV
omicsone cnv-correlation clean-figures --config configs\HNCC_CNV_clean_figures.ini
```

Compatibility entrypoints are also available for CNV correlation:

```powershell
omicsone-replay-cnv-correlation --config configs\HNCC_CNV.ini
omicsone-replay-cnv-correlation-clean-figures --config configs\HNCC_CNV_clean_figures.ini
```

## Config Examples

Minimal differential config:

```ini
[input]
normal_path = C:\data\normal.tsv
tumor_path = C:\data\tumor.tsv
fasta_path = C:\data\protein.fasta
cohort = HNSCC
omics = Protein

[output]
output_dir = C:\runs\HNSCC_Protein
```

Minimal CNV correlation config:

```ini
[task]
cohort = HNSCC

[paths]
cnv_path = C:\data\cnv.tsv
rna_path = C:\data\rna.tsv
protein_path = C:\data\protein.tsv
gistic_path = C:\data\gistic.tsv
fasta_file = C:\data\protein.fasta
chromosomes_file = C:\data\chromosomes.xlsx
cytoband_file = C:\data\cytoBand.txt
output_dir = C:\runs\HNCC_CNV

[settings]
min_valid_pairs = 4
correlation_threshold = 0.5
chunksize = 50000
dpi = 600
use_three_way_common_genes = true
generate_clean_figures = true
```

Different workflows require different fields. If required fields are missing,
the CLI reports which settings are missing.

## Docker Replay

The Docker image fixes the Python environment to Python 3.11 and builds the
Rust/PyO3 Spearman extension during image build.

### Build The Image

Run this from the OmicsOne repository:

```powershell
cd <path-to-omicsone-streamlit>
docker build -t omicsone:py311-cnv .
```

The build should print:

```text
spearman backend: rust
```

You can also verify the backend manually:

```powershell
docker run --rm --entrypoint python omicsone:py311-cnv -c "from omicsone.utils import spearmanr; print(spearmanr.backend())"
```

Expected output:

```text
rust
```

### Replay Package Layout

The Docker image expects a run folder mounted to `/runs/current`. The minimum
replay package is:

```text
<run-folder>/
|-- data/
|-- configs.docker/
`-- replay_all.sh
```

`data/` alone is not enough. `configs.docker/` stores paths, thresholds,
enrichment settings, output directories, and CNV correlation settings.
`replay_all.sh` preserves the correct run order.

Generated output folders can be deleted and regenerated:

```text
HNCC_CNV/
LSCC_CNV/
HNSCC_Protein/
LSCC_Protein/
HNSCC_phospho/
LSCC_phospho/
HNSCC_Mutations/
LSCC_Mutations/
HNSCC_Protein_boxplots/
LSCC_Protein_boxplots/
HNSCC_pathways/
LSCC_pathways/
```

### Replay All Analyses

Mount any compatible run folder to `/runs/current`:

```powershell
docker run --rm `
  -e RUN_ROOT=/runs/current `
  -v "<run-folder>:/runs/current" `
  omicsone:py311-cnv
```

The default image command runs `/app/replay_all.sh`. A full replay runs:

1. HNSCC protein differential analysis
2. LSCC protein differential analysis
3. HNSCC phospho differential analysis
4. LSCC phospho differential analysis
5. HNSCC mutation figures
6. LSCC mutation figures
7. HNSCC protein boxplots
8. LSCC protein boxplots
9. HNSCC pathway scatter plots
10. LSCC pathway scatter plots
11. HNSCC CNV correlation pipeline
12. LSCC CNV correlation pipeline

To run the script from the mounted folder instead:

```powershell
docker run --rm `
  -e RUN_ROOT=/runs/current `
  -v "<run-folder>:/runs/current" `
  --entrypoint bash `
  omicsone:py311-cnv `
  /runs/current/replay_all.sh
```

### Run One Analysis Only

CNV correlation only:

```powershell
docker run --rm `
  -e RUN_ROOT=/runs/current `
  -v "<run-folder>:/runs/current" `
  omicsone:py311-cnv `
  omicsone cnv-correlation run `
    --config /runs/current/configs.docker/HNCC_CNV.ini `
    --output-dir /runs/current/HNCC_CNV `
    --quiet
```

Phospho differential only:

```powershell
docker run --rm `
  -e RUN_ROOT=/runs/current `
  -v "<run-folder>:/runs/current" `
  omicsone:py311-cnv `
  omicsone differential phospho `
    --config /runs/current/configs.docker/LSCC_phospho.ini `
    --quiet
```

## Enrichment Modes

For stable offline replay, configs can point to copied reference enrichment
TSVs:

```ini
[reference]
up_enrichment_tsv = /runs/current/data/references/<cohort>/<up_enrichment>.tsv
down_enrichment_tsv = /runs/current/data/references/<cohort>/<down_enrichment>.tsv
```

For local GMT enrichment without network access, clear the reference TSVs and
use a local GMT file:

```ini
[enrichment]
gene_sets = /runs/current/data/references/Enrichr.MSigDB_Hallmark_2020.gmt
enrichment_background_mode = gene_list
prefer_local_gene_sets = true
allow_remote_enrichr = false

[reference]
up_enrichment_tsv =
down_enrichment_tsv =
```

For remote Enrichr, use:

```ini
[enrichment]
enrichment_background_mode = online
prefer_local_gene_sets = false
allow_remote_enrichr = true

[reference]
up_enrichment_tsv =
down_enrichment_tsv =
```

Remote Enrichr submits gene lists to an external service, so use it only when
that is acceptable for the data.

## Fonts

The default plotting font is `Liberation Sans`. It is an open-source,
Arial-compatible font and is installed in the Docker image through the Debian
`fonts-liberation` package.

To change the plotting font, edit the relevant config file:

```ini
[plot]
font_family = Liberation Sans
```

You can replace `Liberation Sans` with any font installed in the runtime
environment, for example:

```ini
[plot]
font_family = DejaVu Sans
```

For Docker replay, put the setting in the relevant file under
`configs.docker/`, such as:

```text
configs.docker/HNSCC_Protein.ini
configs.docker/LSCC_phospho.ini
configs.docker/HNCC_CNV.ini
```

If a private run folder contains fonts under `data/fonts`, the Docker entry
script registers them at container startup:

```text
<run-folder>/data/fonts/
```

Do not publish Microsoft Arial font files in a public Docker image or public
data package. If exact Arial rendering is required, users should provide their
own properly licensed Arial files privately and set:

```ini
[plot]
font_family = Arial
```

## Validation Files

A reproducible run folder can include validation outputs such as:

```text
validation_summary.csv
validation_summary_normalized.csv
validation_diff_tolerance_summary.csv
```

## Streamlit

Coming soon.

## FastAPI

Coming soon.

## Development

Run tests:

```powershell
python -m pytest
```

Run only CLI tests:

```powershell
python -m pytest tests\test_cli_dispatch.py tests\test_entrypoints.py
```
