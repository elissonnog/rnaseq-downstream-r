#!/usr/bin/env Rscript

config_path <- normalizePath(file.path("config", "gse231397_e2_vs_vehicle.R"))
rmarkdown::render(
  input = file.path("report", "rnaseq_downstream.Rmd"),
  output_file = "gse231397_e2_vs_vehicle_downstream.html",
  output_dir = normalizePath("report"),
  params = list(config_path = config_path, portfolio_snapshot = TRUE),
  knit_root_dir = getwd(), envir = new.env(parent = globalenv()), quiet = FALSE
)
