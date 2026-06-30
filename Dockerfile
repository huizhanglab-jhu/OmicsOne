FROM python:3.11-slim AS builder

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV CARGO_HOME=/opt/cargo
ENV RUSTUP_HOME=/opt/rustup
ENV PATH="/opt/cargo/bin:${PATH}"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

RUN curl https://sh.rustup.rs -sSf | sh -s -- -y --profile minimal

WORKDIR /app

COPY pyproject.toml setup.py README.md LICENSE VERSION ./
COPY src ./src
COPY packages/rust_spearmanr ./packages/rust_spearmanr

RUN python -m pip install --upgrade pip setuptools wheel setuptools-rust \
    && python -m pip wheel . -w /wheels


FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV RUN_ROOT=/runs/current

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        fonts-liberation \
        fontconfig \
        libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /wheels /wheels
COPY docker/replay_all.sh /app/replay_all.sh

RUN python -m pip install --upgrade pip \
    && python -m pip install --no-cache-dir --no-index --find-links=/wheels omicsone \
    && rm -rf /wheels \
    && chmod +x /app/replay_all.sh \
    && python - <<'PY'
from omicsone.utils import spearmanr
print("spearman backend:", spearmanr.backend())
if spearmanr.backend() != "rust":
    raise SystemExit("Rust Spearman backend was not compiled into the image")
PY

EXPOSE 8000

CMD ["/app/replay_all.sh"]
