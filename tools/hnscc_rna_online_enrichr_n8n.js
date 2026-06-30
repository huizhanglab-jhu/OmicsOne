const payload = {
  "normal_path": "F:\\lab\\HsinI\\Head and Neck & Lung\\HNSCC\\HNSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt",
  "tumor_path": "F:\\lab\\HsinI\\Head and Neck & Lung\\HNSCC\\HNSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
  "output_dir": "E:\\lab\\HSinI\\runs\\20260503_HNSCC_RNA_volcano",
  "fasta_path": "F:\\lab\\HsinI\\Head and Neck & Lung\\fasta\\GENCODE.V42.basic.CHR.combined_contaminants.gpquest3.fasta",
  "cohort": "HNSCC",
  "omics": "RNA",
  "job_name": "HNSCC_RNA",
  "method": "Wilcoxon(Unpaired)",
  "fdr_cutoff": 0.01,
  "log2fc_cutoff": 1.0,
  "max_miss_ratio_global": 0.5,
  "max_miss_ratio_group": 0.5,
  "min_sample_size": 4,
  "gene_sets": ["MSigDB_Hallmark_2020"],
  "enrichment_fdr_cutoff": 0.05,
  "organism": "human",
  "enrichment_background_mode": "online",
  "prefer_local_gene_sets": false,
  "allow_remote_enrichr": true,
  "skip_pathways": [
    "Phagosome",
    "Human papillomavirus infection",
    "Pertussis",
    "Malaria",
    "Arrhythmogenic right ventricular cardiomyopathy",
    "Staphylococcus aureus infection",
    "Regulation of actin cytoskeleton"
  ],
  "enrichment_top_n": 10,
  "title": "HNSCC RNA differential expression analysis",
  "enrichment_title": "HNSCC RNA MSigDB hallmark Pathways Enrichment Analysis",
  "xlabel": "Log2FC(Tumor/NAT)",
  "volcano_width": 4.0,
  "volcano_height": 4.0,
  "enrichment_width": null,
  "enrichment_height": null,
  "enrichment_width_ratio": 2.0,
  "enrichment_height_ratio": 1.5,
  "enrichment_size_scale": 4.0,
  "enrichment_min_x": -25.0,
  "enrichment_max_x": 60.0,
  "dpi": 300,
  "font_family": "Arial",
  "font_size": 10.0,
  "editable_pdf_text": true,
  "background_color": "#808080",
  "up_color": "#FF0000",
  "down_color": "#0000FF",
  "point_size": 1.0,
  "significant_point_size": 5.0,
  "output_prefix": null
};

const response = await this.helpers.httpRequest({
  method: "POST",
  url: "http://127.0.0.1:8001/api/v1/diff/volcano/enrichment",
  headers: {
    "Content-Type": "application/json"
  },
  body: payload,
  json: true
});

return response;
