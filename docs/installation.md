# Installation

## Editable Development Install

From the repository root:

```powershell
C:\Users\yhu39\AppData\Local\anaconda3\envs\omicsone\python.exe -m pip install -e .
```

The distribution name is `omicsone`.

## Commands

Run the Streamlit app:

```powershell
omicsone
```

Run the FastAPI service:

```powershell
omicsone-api
```

Run the installed Replay CNV correlation command:

```powershell
omicsone-replay-cnv-correlation --help
```

## Rust Spearman Backend

OmicsOne includes a Rust/PyO3 Spearman correlation extension. The source is in:

```text
packages/rust_spearmanr/src/lib.rs
```

The extension is built into:

```python
omicsone.utils._spearmanr
```

The Python wrapper is:

```python
from omicsone.utils import spearmanr
```

For source installs, users need a working Rust/Cargo toolchain and the platform
compiler required by PyO3. For releases to non-developer users, prefer
publishing prebuilt wheels.

