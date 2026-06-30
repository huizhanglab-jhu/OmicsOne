from omicsone.services.volcano_enrichment import (
    generate_volcano_enrichment,
    resolve_volcano_preset,
)


def main() -> None:
    for preset in ["HNSCC_RNA", "HNSCC_Protein", "LSCC_RNA", "LSCC_Protein"]:
        payload = resolve_volcano_preset(preset)
        print(f"RUN {preset} -> {payload['output_dir']}", flush=True)
        result = generate_volcano_enrichment(**payload)
        print(
            {
                "preset": preset,
                "output_dir": str(result.output_dir),
                "diff_feature_count": result.diff_feature_count,
                "up_count": result.up_count,
                "down_count": result.down_count,
                "up_enrichment_count": result.up_enrichment_count,
                "down_enrichment_count": result.down_enrichment_count,
                "result_log": str(result.result_log),
            },
            flush=True,
        )


if __name__ == "__main__":
    main()
