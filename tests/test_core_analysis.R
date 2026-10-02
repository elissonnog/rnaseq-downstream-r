#!/usr/bin/env Rscript

source(file.path("R", "helpers.R"))
assert_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)

config <- load_analysis_config(file.path("config", "example_config.R"))
counts <- read_count_matrix(config$counts_path)
metadata <- read_sample_metadata(config$metadata_path)
contrasts <- read_contrasts(config$contrasts_path)
metadata <- validate_analysis_inputs(counts, metadata, contrasts, config$design)
dds <- fit_deseq_model(prepare_deseq_dataset(counts, metadata, config$design, config$min_total_count))
result <- run_deseq_contrast(dds, contrasts[1L, , drop = FALSE], config$fdr,
                             input_counts = counts, min_total_count = config$min_total_count)

assert_true(nrow(result) == nrow(counts), "The full input gene ledger was not preserved.")
assert_true(all(result$test_status %in% c(
  "all_zero", "low_count_filtered", "unavailable", "independent_filtered",
  "tested_nonsignificant", "significant"
)), "Unexpected result status.")
assert_true(all(is.na(result$significant[result$test_status %in% c(
  "all_zero", "low_count_filtered", "unavailable", "independent_filtered"
)])), "Untested or independently filtered genes must have significant = NA.")
assert_true(all(!result$significant[result$test_status == "tested_nonsignificant"]),
            "Tested nonsignificant genes must have significant = FALSE.")
assert_true(sum(result$test_status == "significant") >= 30L,
            "Synthetic planted signal recovery changed unexpectedly.")
summary <- summarize_contrast_result(result, contrasts$label[[1L]])
assert_true(summary$input_rows == nrow(counts), "Contrast summary input count is wrong.")
assert_true(summary$retained_rows == nrow(dds), "Contrast summary retained count is wrong.")
cat("RNA-seq DESeq2 integration test passed.\n")
