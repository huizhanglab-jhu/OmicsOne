# CLI And Helper Tools

## Installed Commands

The package exposes console scripts from `pyproject.toml`.

### `omicsone`

Runs the local OmicsOne CLI for terminal, replay, n8n, and batch workflows.

```powershell
omicsone --help
omicsone differential run --config path\to\input.ini
omicsone cnv-correlation run --config path\to\input.ini
```

Common workflow commands:

```powershell
omicsone differential run --config path\to\input.ini
omicsone differential phospho --config path\to\input.ini
omicsone mutations run --config path\to\input.ini
omicsone mutations post-api --config path\to\input.ini --api-url http://127.0.0.1:8001/api/v1/mutations/heatmap/figures
omicsone boxplots run --config path\to\input.ini
omicsone pathway-scatter run --config path\to\input.ini
omicsone pathway-scatter phosphosite-protein --config path\to\input.ini
omicsone cnv-correlation run --config path\to\input.ini
```

### `omicsone-app`

Runs the Streamlit app. This replaces the old `omicsone` app-launch behavior.

```powershell
omicsone-app
```

The app can also be launched from the CLI:

```powershell
omicsone app
```

### `omicsone-api`

Runs the FastAPI service.

```powershell
omicsone-api
```

### `omicsone-replay-cnv-correlation`

Runs the replay-friendly CNV correlation pipeline.

```powershell
omicsone-replay-cnv-correlation --help
```

## Repository Helper Scripts

These scripts live under `tools/` and are usually run from the repository root.

### `tools/run_volcano_enrichment_config.py`

Runs differential volcano/enrichment analysis from a config file.

```powershell
python tools/run_volcano_enrichment_config.py --config path\to\config.ini
```

### `tools/run_pathway_scatter_config.py`

Runs protein/phosphosite pathway scatter plots from a config file.

```powershell
python tools/run_pathway_scatter_config.py --config path\to\config.ini
```

### `tools/run_volcano_presets.py`

Runs predefined volcano/enrichment presets from
`src/omicsone/services/volcano_enrichment.py`.

### n8n JavaScript Payloads

Several `tools/*_online_enrichr_n8n.js` files are replay examples for posting
volcano/enrichment payloads to the local API endpoint:

```text
http://127.0.0.1:8001/api/v1/diff/volcano/enrichment
```

Generated differential and boxplot runs write n8n replay scripts only when
`write_n8n_script = true` is set in the config or request payload.

