from __future__ import annotations

from pyteomics import fasta


def get_gene_map(fasta_path: str) -> dict:
    chr_gene_map = {}
    for description, _sequence in fasta.read(fasta_path):
        items = description.split("|")
        try:
            gene_id = [item for item in items if item[:4] == "ENSG"][0].split(".")[0]
            gene = [item for item in items if item[:3] == "GN="][0].split("=")[1]
            chrom = [item for item in items if item[:3] == "chr"][0].split("-")[0].split(":")
            chrom[0] = chrom[0][3:]
            chrom[1] = int(chrom[1])
            chr_gene_map[gene_id] = {
                "chr": str(chrom[0]),
                "offset": int(chrom[1]),
                "gene": gene,
            }
        except (IndexError, ValueError):
            continue

    return chr_gene_map
