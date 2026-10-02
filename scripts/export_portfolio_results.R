#!/usr/bin/env Rscript

source(file.path("R", "helpers.R"))
config <- load_analysis_config(file.path("config", "gse231397_e2_vs_vehicle.R"))
source_dir <- file.path(config$output_dir, "tables")
destination <- file.path("results", "gse231397")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(source_dir, "contrast_summary.csv")
result_path <- file.path(source_dir, "E2_vs_vehicle_deseq2.csv")
if (!file.exists(summary_path) || !file.exists(result_path)) {
  stop("Run the GSE231397 analysis before exporting portfolio results.")
}
summary <- read.csv(summary_path, check.names = FALSE, stringsAsFactors = FALSE)
result <- read.csv(result_path, check.names = FALSE, stringsAsFactors = FALSE)
write.csv(summary, file.path(destination, "contrast_summary.csv"), row.names = FALSE)

targets <- data.frame(
  gene = c("GREB1", "PGR", "TFF1"), gene_id = c("9687", "5241", "7031"),
  expected_direction = "up", stringsAsFactors = FALSE
)
matched <- result[match(targets$gene_id, as.character(result$gene_id)), ]
if (anyNA(matched$gene_id)) stop("One or more declared published targets are absent from the result ledger.")
targets$log2FoldChange <- matched$log2FoldChange
targets$padj <- matched$padj
targets$test_status <- matched$test_status
targets$observed_direction <- matched$direction
targets$direction_concordant <- matched$test_status == "significant" & matched$direction == "up"
write.csv(targets, file.path(destination, "published_target_concordance.csv"), row.names = FALSE)
cat("Regenerated GSE231397 portfolio summaries.\n")
