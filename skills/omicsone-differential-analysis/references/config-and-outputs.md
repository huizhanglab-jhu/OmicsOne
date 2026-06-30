# Differential Config And Outputs

## Config Sections

`src/omicsone/replay/default_differential_config.ini` is the packaged template.
The parser accepts settings from sections such as `[analysis]`, `[diff]`,
`[enrichment]`, `[plot]`, and `[style]`.

Important analysis keys:

- `method`: `Wilcoxon(Unpaired)`, `Wilcoxon(Paired)`, `T-test(Unpaired)`, or
  `T-test(Paired)`.
- `strip_feature_version`: true for gene-level Ensembl IDs; false for
  phosphosite-style IDs.

Important differential keys:

- `fdr_cutoff`
- `log2fc_cutoff`
- `max_miss_ratio_global`
- `max_miss_ratio_group`
- `min_sample_size`

Important enrichment keys:

- `gene_sets`: Enrichr library names or local GMT paths.
- `enrichment_background_mode`: `gene_list`, `online`, `notebook`, or `count`.
- `prefer_local_gene_sets`
- `allow_remote_enrichr`
- `skip_pathways`
- `enrichment_top_n`

Important plotting keys:

- `xlabel`
- `volcano_width`, `volcano_height`
- `enrichment_width`, `enrichment_height`
- `dpi`, `tiff_dpi`
- `font_family`, `font_size`
- `background_color`, `up_color`, `down_color`

## Output Meaning

- `diff.tsv`: feature-level test statistics, p-values, FDR, Log2FC, and
  significance labels.
- `combined_matrix.tsv`: tumor and normal matrices after joining and suffixing
  columns.
- `pure_up_genes.tsv`: mapped gene symbols significant only in the up direction.
- `pure_down_genes.tsv`: mapped gene symbols significant only in the down
  direction.
- `total_genes.tsv`: mapped background gene symbols.
- `up_enrichr_df.tsv`, `down_enrichr_df.tsv`: enrichment tables.
- `enrichment_plot_table.tsv`: filtered terms actually plotted.
- `result.log`: key parameters, counts, and artifact paths.
- `report.html`: compact report with parameters, plots, and enrichment tables.

## Troubleshooting

- Empty `diff.tsv`: check file paths, matrix index column, missing values, and
  sample counts.
- Empty enrichment tables: check FASTA gene mapping, feature ID format, gene set
  availability, and `allow_remote_enrichr`.
- Wrong Log2FC direction: verify which file was passed as tumor and which as
  normal.
- Phosphosite IDs collapsed: set `strip_feature_version=false`.

