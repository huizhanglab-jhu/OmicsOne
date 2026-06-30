# API Overview

OmicsOne exposes a FastAPI app from:

```text
src/omicsone/api/app.py
```

Start it with:

```powershell
omicsone-api
```

The API is organized by routers:

| Area | Prefix | Router |
| --- | --- | --- |
| Health | `/health` | `src/omicsone/api/routers/health.py` |
| FASTA | `/api/v1/fasta` | `src/omicsone/api/routers/fasta.py` |
| Mutations v1 | `/api/v1/mutations` | `src/omicsone/api/routers/mutations.py` |
| Mutations v2 | `/api/v2/mutations` | `src/omicsone/api/routers/mutations_v2.py` |
| Spearman | `/api/v1/spearman` | `src/omicsone/api/routers/spearman.py` |
| CNV correlation | `/api/v1/cnv-correlation` | `src/omicsone/api/routers/cnv_correlation.py` |
| Differential volcano | `/api/v1/diff` | `src/omicsone/api/routers/volcano.py` |
| Differential boxplots | `/api/v1/diff` | `src/omicsone/api/routers/boxplots.py` |

Request and response models live in:

```text
src/omicsone/api/schemas/
```

The API generally returns artifact paths rather than large binary payloads.
Generated tables and figures are written to the requested `output_dir`.
