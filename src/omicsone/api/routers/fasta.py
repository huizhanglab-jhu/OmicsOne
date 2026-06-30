import configparser
from pathlib import Path
from typing import Dict, Optional

from fastapi import APIRouter, HTTPException

from omicsone.api.schemas.fasta import FastaGeneMapRequest, FastaGeneMapResponse
from omicsone.services.fasta import get_gene_map


router = APIRouter()


def get_configured_fasta_path() -> str:
    settings_path = Path(__file__).resolve().parents[3] / "omicsone_streamlit" / "config" / "settings.ini"
    settings = configparser.ConfigParser()
    settings.read(settings_path, encoding="utf-8")
    try:
        return settings["paths"]["fasta_path"]
    except KeyError as exc:
        raise HTTPException(
            status_code=500,
            detail="Missing paths.fasta_path in settings.ini",
        ) from exc


def resolve_fasta_path(fasta_path: Optional[str]) -> Path:
    path = Path(fasta_path or get_configured_fasta_path()).expanduser()
    if not path.is_file():
        raise HTTPException(
            status_code=404,
            detail=f"FASTA file does not exist: {path}",
        )
    return path


@router.get("/default-path")
def get_default_fasta_path():
    return {"fasta_path": get_configured_fasta_path()}


@router.post("/gene-map", response_model=FastaGeneMapResponse)
def build_gene_map(request: FastaGeneMapRequest):
    fasta_path = resolve_fasta_path(request.fasta_path)
    gene_map = get_gene_map(str(fasta_path))

    if request.limit is not None:
        limited_items = list(gene_map.items())[: request.limit]
        result: Dict[str, dict] = dict(limited_items)
    else:
        result = gene_map

    return {
        "fasta_path": str(fasta_path),
        "total": len(gene_map),
        "gene_map": result,
    }
