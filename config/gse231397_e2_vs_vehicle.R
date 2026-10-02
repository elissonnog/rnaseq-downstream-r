analysis_config <- list(
  project_title = "GSE231397: estradiol response in MCF-7 cells",
  data_label = "GSE231397: 1 nM E2 vs vehicle (24 h)",
  counts_path = "example/gse231397/data/counts.csv",
  metadata_path = "example/gse231397/metadata.csv",
  contrasts_path = "example/gse231397/contrasts.csv",
  output_dir = "output/gse231397_e2_vs_vehicle",
  design = "~ condition",
  min_total_count = 10L,
  fdr = 0.05,
  pca_top_genes = 500L,
  heatmap_top_genes = 30L,
  enrichment = list(
    enabled = FALSE,
    species = "Homo sapiens",
    gene_id_type = "ENTREZID",
    organism_db = "org.Hs.eg.db"
  )
)
