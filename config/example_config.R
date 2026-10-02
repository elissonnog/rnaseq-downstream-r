analysis_config <- list(
  project_title = "Synthetic bulk RNA-seq mechanics check",
  data_label = "SYNTHETIC - mechanics check only",
  counts_path = "example/synthetic/counts.csv",
  metadata_path = "example/synthetic/metadata.csv",
  contrasts_path = "example/synthetic/contrasts.csv",
  output_dir = "output/synthetic",
  design = "~ condition",
  min_total_count = 10L,
  fdr = 0.05,
  pca_top_genes = 30L,
  heatmap_top_genes = 20L,
  gene_labels = list(
    enabled = FALSE,
    gene_id_type = "SYMBOL",
    organism_db = "org.Hs.eg.db",
    column = "SYMBOL"
  ),
  enrichment = list(
    enabled = FALSE,
    species = "Homo sapiens",
    gene_id_type = "SYMBOL",
    organism_db = "org.Hs.eg.db"
  )
)
