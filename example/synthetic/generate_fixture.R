#!/usr/bin/env Rscript

# Reproducible synthetic counts for interface and mechanics testing only.
# These values do not represent a real organism, experiment, or biological result.
set.seed(20261002)

gene_count <- 300L
sample_ids <- sprintf("Sample%02d", seq_len(6L))
conditions <- rep(c("control", "treatment"), each = 3L)
base_means <- exp(stats::rnorm(gene_count, mean = log(80), sd = 1.0))
fold_change <- rep(1, gene_count)
fold_change[1:20] <- 4
fold_change[21:40] <- 0.25
library_scale <- c(0.90, 1.05, 1.12, 0.95, 1.08, 1.02)

counts <- vapply(seq_along(sample_ids), function(index) {
  condition_effect <- if (conditions[index] == "treatment") fold_change else 1
  stats::rnbinom(
    gene_count,
    mu = base_means * condition_effect * library_scale[index],
    size = 5
  )
}, numeric(gene_count))

fixture <- data.frame(
  gene_id = sprintf("Gene%03d", seq_len(gene_count)),
  counts,
  check.names = FALSE
)
names(fixture)[-1L] <- sample_ids

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
fixture_dir <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
write.csv(fixture, file.path(fixture_dir, "counts.csv"), row.names = FALSE, quote = FALSE)

