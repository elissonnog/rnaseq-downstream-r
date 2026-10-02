#!/usr/bin/env Rscript

source(file.path("R", "helpers.R"))

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

assert_error <- function(expression, pattern) {
  message <- tryCatch({
    force(expression)
    NA_character_
  }, error = function(error) conditionMessage(error))
  if (is.na(message) || !grepl(pattern, message, fixed = TRUE)) {
    stop("Expected error containing `", pattern, "`; received: ", message, call. = FALSE)
  }
}

counts <- read_count_matrix(file.path("example", "synthetic", "counts.csv"))
metadata <- read_sample_metadata(file.path("example", "synthetic", "metadata.csv"))
contrasts <- read_contrasts(file.path("example", "synthetic", "contrasts.csv"))
metadata <- validate_analysis_inputs(counts, metadata, contrasts, "~ condition")

assert_true(identical(dim(counts), c(300L, 6L)), "Synthetic count dimensions changed unexpectedly.")
assert_true(identical(colnames(counts), rownames(metadata)), "Sample ordering was not aligned.")
assert_true(identical(safe_file_stem("treatment vs control"), "treatment_vs_control"), "Filename sanitization failed.")

reordered <- metadata[rev(rownames(metadata)), , drop = FALSE]
realigned <- validate_analysis_inputs(counts, reordered, contrasts, "~ condition")
assert_true(identical(colnames(counts), rownames(realigned)), "Metadata was not safely reordered.")

bad_contrasts <- contrasts
bad_contrasts$numerator <- "absent_level"
assert_error(
  validate_analysis_inputs(counts, metadata, bad_contrasts, "~ condition"),
  "Contrast levels absent"
)

assert_error(
  validate_analysis_inputs(counts, metadata, contrasts, "~ missing_variable"),
  "Design variables absent"
)

assert_error(
  validate_analysis_inputs(counts, metadata, contrasts, "~ batch"),
  "Contrast factor absent from design formula"
)

assert_error(
  validate_analysis_inputs(counts, metadata, contrasts, "~ condition * batch"),
  "Only additive categorical design terms"
)

numeric_metadata <- metadata
numeric_metadata$continuous <- seq_len(nrow(numeric_metadata))
assert_error(
  validate_analysis_inputs(counts, numeric_metadata, contrasts, "~ condition + continuous"),
  "Design variables must be categorical"
)

collision_file <- tempfile(fileext = ".csv")
writeLines(c(
  "label,factor,numerator,denominator",
  "A/B,condition,treatment,control",
  "A B,condition,treatment,control"
), collision_file)
assert_error(read_contrasts(collision_file), "colliding output filenames")

case_collision_file <- tempfile(fileext = ".csv")
writeLines(c(
  "label,factor,numerator,denominator",
  "Treatment,condition,treatment,control",
  "treatment,condition,treatment,control"
), case_collision_file)
assert_error(read_contrasts(case_collision_file), "colliding output filenames")

self_contrast_file <- tempfile(fileext = ".csv")
writeLines(c(
  "label,factor,numerator,denominator",
  "control_vs_control,condition,control,control"
), self_contrast_file)
assert_error(read_contrasts(self_contrast_file), "numerator and denominator must differ")

duplicate_header_file <- tempfile(fileext = ".csv")
writeLines(c("gene_id,Sample01,Sample01", "Gene001,1,2"), duplicate_header_file)
assert_error(read_count_matrix(duplicate_header_file), "column headers must be unique")

blank_header_file <- tempfile(fileext = ".csv")
writeLines(c("gene_id,,Sample02", "Gene001,1,2"), blank_header_file)
assert_error(read_count_matrix(blank_header_file), "column headers must be non-blank")

duplicate_metadata_header_file <- tempfile(fileext = ".csv")
writeLines(c("sample_id,condition,condition", "Sample01,control,A"), duplicate_metadata_header_file)
assert_error(read_sample_metadata(duplicate_metadata_header_file), "Metadata column headers must be unique")

blank_metadata_header_file <- tempfile(fileext = ".csv")
writeLines(c("sample_id,,batch", "Sample01,control,A"), blank_metadata_header_file)
assert_error(read_sample_metadata(blank_metadata_header_file), "Metadata column headers must be non-blank")

duplicate_contrast_header_file <- tempfile(fileext = ".csv")
writeLines(c(
  "label,factor,numerator,numerator",
  "comparison,condition,treatment,control"
), duplicate_contrast_header_file)
assert_error(read_contrasts(duplicate_contrast_header_file), "Contrast table column headers must be unique")

blank_contrast_header_file <- tempfile(fileext = ".csv")
writeLines(c(
  "label,factor,,denominator",
  "comparison,condition,treatment,control"
), blank_contrast_header_file)
assert_error(read_contrasts(blank_contrast_header_file), "Contrast table column headers must be non-blank")

config <- load_analysis_config(file.path("config", "example_config.R"))
assert_true(identical(config$design, "~ condition"), "Configuration did not load correctly.")
assert_true(grepl("SYNTHETIC", config$data_label, fixed = TRUE), "Synthetic status label is missing.")
assert_true(!isTRUE(config$enrichment$enabled), "Synthetic enrichment must remain disabled.")

ranking_fixture <- data.frame(
  gene_id = c("gene_a", "gene_b", "gene_c", "gene_d"),
  padj = c(0.01, 0.02, 0.01, 0.20),
  log2FoldChange = c(1, 5, -2, 8),
  stringsAsFactors = FALSE
)
assert_true(
  identical(select_significant_genes(ranking_fixture, 0.05, 2L), c("gene_c", "gene_a")),
  "Significant-gene heatmap ranking is not deterministic."
)
assert_true(
  identical(select_significant_genes(ranking_fixture, 0.001, 10L), character()),
  "Heatmap selection must not include genes outside the declared adjusted-p-value threshold."
)

cat("All dependency-free helper and fixture tests passed.\n")
