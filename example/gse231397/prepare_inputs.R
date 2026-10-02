#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: prepare_inputs.R <counts.tsv.gz> <family.soft.gz> <metadata.csv> <output.csv>")
}

counts_path <- args[[1L]]
soft_path <- args[[2L]]
metadata_path <- args[[3L]]
output_path <- args[[4L]]

vehicle_ids <- c("GSM7270183", "GSM7270188", "GSM7270193")
estradiol_ids <- c("GSM7270186", "GSM7270190", "GSM7270195")
selected_ids <- c(vehicle_ids, estradiol_ids)

counts <- read.delim(gzfile(counts_path), check.names = FALSE, stringsAsFactors = FALSE)
if (!identical(names(counts)[[1L]], "GeneID")) stop("Expected `GeneID` as the first count column.")
if (nrow(counts) != 39376L || ncol(counts) != 13L) stop("Expected 39,376 genes across 12 samples.")
if (anyDuplicated(counts$GeneID)) stop("NCBI count matrix contains duplicated GeneIDs.")
sample_columns <- counts[-1L]
valid_counts <- vapply(sample_columns, function(x) all(is.finite(x) & x >= 0 & x == floor(x)), logical(1))
if (!all(valid_counts)) stop("NCBI count matrix contains non-integer, negative, or missing values.")
if (!all(selected_ids %in% names(counts))) stop("One or more declared samples are absent from the count matrix.")

soft <- readLines(gzfile(soft_path), warn = FALSE)
if (!any(soft == "^SERIES = GSE231397")) stop("SOFT file is not GSE231397.")

sample_field <- function(sample_id, prefix) {
  start <- match(paste0("^SAMPLE = ", sample_id), soft)
  if (is.na(start)) stop("Sample absent from SOFT: ", sample_id)
  later_samples <- which(seq_along(soft) > start & grepl("^\\^SAMPLE = ", soft))
  end <- if (length(later_samples)) later_samples[[1L]] - 1L else length(soft)
  block <- soft[start:end]
  values <- sub(prefix, "", block[grepl(paste0("^", prefix), block)])
  if (!length(values)) stop("Missing SOFT field for ", sample_id, ": ", prefix)
  values
}

expected_titles <- c(
  GSM7270183 = "MCF7-Veh biol rep 1",
  GSM7270188 = "MCF7-Veh biol rep 2",
  GSM7270193 = "MCF7-Veh biol rep 3",
  GSM7270186 = "MCF7-E2 biol rep 1",
  GSM7270190 = "MCF7-E2 biol rep 2",
  GSM7270195 = "MCF7-E2 biol rep 3"
)
expected_treatments <- c(
  GSM7270183 = "DMSO", GSM7270188 = "DMSO", GSM7270193 = "DMSO",
  GSM7270186 = "1 nM E2", GSM7270190 = "1 nM E2", GSM7270195 = "1 nM E2"
)
for (sample_id in selected_ids) {
  title <- sample_field(sample_id, "!Sample_title = ")[[1L]]
  treatment <- sample_field(sample_id, "!Sample_characteristics_ch1 = treatment: ")[[1L]]
  if (!identical(title, expected_titles[[sample_id]])) stop("Unexpected SOFT title for ", sample_id)
  if (!identical(treatment, expected_treatments[[sample_id]])) stop("Unexpected SOFT treatment for ", sample_id)
}

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
required_metadata <- c("sample_id", "condition", "biological_replicate", "treatment")
if (!identical(names(metadata), required_metadata)) stop("Unexpected metadata columns.")
if (!identical(metadata$sample_id, selected_ids)) stop("Metadata sample order or membership changed.")
if (!identical(metadata$condition, c(rep("vehicle", 3L), rep("estradiol", 3L)))) stop("Metadata conditions changed.")
if (!identical(metadata$biological_replicate, rep(1:3, 2L))) stop("Metadata replicate labels changed.")

prepared <- counts[, c("GeneID", selected_ids), drop = FALSE]
names(prepared)[[1L]] <- "gene_id"
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(prepared, output_path, row.names = FALSE, quote = FALSE)
cat("Verified 39,376 unique GeneIDs across 12 NCBI count columns; selected six independently verified samples.\n")
