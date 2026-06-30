# FastAPI Endpoints

## Health

### `GET /health`

Returns:

```json
{"status": "ok"}
```

## FASTA

### `GET /api/v1/fasta/default-path`

Returns the configured default FASTA path from the packaged settings file.

### `POST /api/v1/fasta/gene-map`

Builds a gene map from a FASTA file.

Request fields:

- `fasta_path`: optional FASTA path. If omitted, the configured default is used.
- `limit`: optional maximum number of records to return. Set `null` to return
  the full map.

## Differential Volcano And Enrichment

### `POST /api/v1/diff/volcano/enrichment`

Runs two-group differential analysis, volcano plotting, gene mapping, and
pathway enrichment.

Required request fields:

- `normal_path`
- `tumor_path`
- `output_dir`

Common optional fields:

- `fasta_path`
- `cohort`
- `omics`
- `method`
- `strip_feature_version`
- `fdr_cutoff`
- `log2fc_cutoff`
- `gene_sets`
- `enrichment_background_mode`
- `allow_remote_enrichr`

### `POST /api/v1/diff/volcano/enrichment/preset`

Runs differential analysis from a named preset in
`src/omicsone/services/volcano_enrichment.py`.

Request fields:

- `preset`
- `overrides`

## Differential Boxplots

### `POST /api/v1/diff/boxplot/figures`

Generates one normal-vs-tumor PDF boxplot per requested gene.

Required request fields:

- `normal_path`
- `tumor_path`
- `output_dir`

Common optional fields:

- `fasta_path`
- `cohort`
- `omics`
- `genes`
- `ylabel`
- `width`
- `height`

## Spearman

### `POST /api/v1/spearman/compute-file`

Computes all row-pair Spearman correlations between two numeric matrix files.

Request fields:

- `input_file1`
- `input_file2`
- `output_file`
- `min_valid_pairs`

### `POST /api/v1/spearman/paired-omics`

Preprocesses two or three omics matrices, intersects common samples/features,
and writes pairwise Spearman outputs.

Request fields:

- `input1`
- `input2`
- `input3`
- `output_dir`
- `output_prefix`
- `min_valid_pairs`

Each input has:

- `path`
- `data_type`: `cnv`, `rna`, or `protein`

## CNV Correlation

### `POST /api/v1/cnv-correlation/pipeline`

Starts the asynchronous CNV correlation pipeline.

Request fields:

- `cohort`
- `cnv_path`
- `rna_path`
- `protein_path`
- `gistic_path`
- `fasta_file`
- `chromosomes_file`
- `cytoband_file`
- `output_dir`
- `min_valid_pairs`
- `correlation_threshold`
- `chunksize`
- `dpi`

At least one of `rna_path` or `protein_path` is required.

### `GET /api/v1/cnv-correlation/pipeline/{job_id}`

Returns async pipeline status and generated report paths.

### `POST /api/v1/cnv-correlation/figures`

Generates CNV correlation figure panels from precomputed Spearman outputs.

Required request fields:

- `correlation_file`
- `cnv_for_corr_file`
- `target_for_corr_file`
- `fasta_file`
- `chromosomes_file`
- `cytoband_file`
- `output_dir`

## Mutations

### `POST /api/v1/mutations/heatmap/figures`

Generates mutation heatmap and mutation-type PDF figures from the v1 mutation
inputs.

### `POST /api/v2/mutations/heatmap/figures`

Generates mutation heatmap and mutation-type PDF figures from binary mutation
matrix plus metadata inputs.
