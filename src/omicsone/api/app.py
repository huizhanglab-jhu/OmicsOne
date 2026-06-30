from fastapi import FastAPI

from omicsone.api.routers import (
    boxplots,
    cnv_correlation,
    fasta,
    health,
    mutations,
    mutations_v2,
    spearman,
    volcano,
)


app = FastAPI(
    title="OmicsOne API",
    version="0.1.0",
    description="Programmatic API for OmicsOne analysis services.",
)

app.include_router(health.router, prefix="/health", tags=["health"])
app.include_router(fasta.router, prefix="/api/v1/fasta", tags=["fasta"])
app.include_router(mutations.router, prefix="/api/v1/mutations", tags=["mutations"])
app.include_router(mutations_v2.router, prefix="/api/v2/mutations", tags=["mutations-v2"])
app.include_router(spearman.router, prefix="/api/v1/spearman", tags=["spearman"])
app.include_router(
    cnv_correlation.router,
    prefix="/api/v1/cnv-correlation",
    tags=["cnv-correlation"],
)
app.include_router(volcano.router, prefix="/api/v1/diff", tags=["diff-volcano"])
app.include_router(boxplots.router, prefix="/api/v1/diff", tags=["diff-boxplots"])
