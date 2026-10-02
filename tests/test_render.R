#!/usr/bin/env Rscript

source(file.path("R", "helpers.R"))
config <- load_analysis_config(file.path("config", "example_config.R"))
config$counts_path <- normalizePath(config$counts_path)
config$metadata_path <- normalizePath(config$metadata_path)
config$contrasts_path <- normalizePath(config$contrasts_path)
config$output_dir <- tempfile("rnaseq-render-")
config_path <- tempfile(fileext = ".R")
writeLines(c("analysis_config <-", capture.output(dput(config))), config_path)
dir.create(config$output_dir, recursive = TRUE)

rmarkdown::render(
  input = file.path("report", "rnaseq_downstream.Rmd"),
  output_file = "rnaseq_downstream.html", output_dir = config$output_dir,
  intermediates_dir = config$output_dir, params = list(config_path = config_path),
  knit_root_dir = getwd(), envir = new.env(parent = globalenv()), quiet = TRUE
)
required <- c(
  file.path(config$output_dir, "rnaseq_downstream.html"),
  file.path(config$output_dir, "tables", "contrast_summary.csv"),
  file.path(config$output_dir, "tables", "warnings.csv"),
  file.path(config$output_dir, "tables", "sample_correlation.csv")
)
if (!all(file.exists(required))) stop("End-to-end render did not create all required artifacts.")
cat("RNA-seq synthetic end-to-end render test passed.\n")
