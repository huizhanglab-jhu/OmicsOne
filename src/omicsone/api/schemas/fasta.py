from typing import Dict, Optional

from pydantic import BaseModel, Field


class FastaGeneMapRequest(BaseModel):
    fasta_path: Optional[str] = Field(
        default=None,
        description="Optional FASTA path. When omitted, the configured default FASTA path is used.",
    )
    limit: Optional[int] = Field(
        default=100,
        ge=1,
        description="Maximum number of gene-map records to return. Set to null to return all records.",
    )


class FastaGeneMapResponse(BaseModel):
    fasta_path: str
    total: int
    gene_map: Dict[str, dict]

