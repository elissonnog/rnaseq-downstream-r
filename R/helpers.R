required_core_packages <- function() {
  c("rmarkdown", "knitr", "DESeq2", "SummarizedExperiment", "ggplot2")
}

required_enrichment_packages <- function() {
  c("AnnotationDbi", "clusterProfiler", "ReactomePA", "org.Hs.eg.db")
}

is_pandoc_available <- function() {
  nzchar(Sys.which("pandoc")) ||
    (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available())
}

check_dependencies <- function(include_enrichment = FALSE) {
  packages <- required_core_packages()
  required_for <- rep("core analysis", length(packages))
  if (isTRUE(include_enrichment)) {
    optional <- setdiff(required_enrichment_packages(), packages)
    packages <- c(packages, optional)
    required_for <- c(required_for, rep("optional enrichment", length(optional)))
  }
  available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  data.frame(
    dependency = c(packages, "Pandoc"),
    kind = c(rep("R package", length(packages)), "renderer"),
    required_for = c(required_for, "core report"),
    available = c(unname(available), is_pandoc_available()),
    row.names = NULL
  )
}

load_analysis_config <- function(path) {
  if (!file.exists(path)) stop("Configuration file not found: ", path)
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  if (!exists("analysis_config", envir = env, inherits = FALSE)) {
    stop("Configuration must define an `analysis_config` list.")
  }
  config <- get("analysis_config", envir = env, inherits = FALSE)
  if (!is.list(config)) stop("`analysis_config` must be a list.")
  required <- c(
    "project_title", "data_label", "counts_path", "metadata_path", "contrasts_path",
    "output_dir", "design", "min_total_count", "fdr",
    "pca_top_genes", "heatmap_top_genes", "enrichment"
  )
  missing <- setdiff(required, names(config))
  if (length(missing)) stop("Missing configuration fields: ", paste(missing, collapse = ", "))
  if (!is.numeric(config$min_total_count) || length(config$min_total_count) != 1L || config$min_total_count < 0) {
    stop("`min_total_count` must be one non-negative number.")
  }
  if (!is.numeric(config$fdr) || length(config$fdr) != 1L || config$fdr <= 0 || config$fdr >= 1) {
    stop("`fdr` must be one number strictly between 0 and 1.")
  }
  for (field in c("pca_top_genes", "heatmap_top_genes")) {
    value <- config[[field]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < 2 || value != as.integer(value)) {
      stop("`", field, "` must be one integer greater than or equal to 2.")
    }
  }
  enrichment_required <- c("enabled", "species", "gene_id_type", "organism_db")
  enrichment_missing <- setdiff(enrichment_required, names(config$enrichment))
  if (length(enrichment_missing)) {
    stop("Missing enrichment configuration fields: ", paste(enrichment_missing, collapse = ", "))
  }
  if (isTRUE(config$enrichment$enabled) &&
      (!identical(config$enrichment$species, "Homo sapiens") ||
       !identical(config$enrichment$organism_db, "org.Hs.eg.db"))) {
    stop("Optional enrichment currently supports only Homo sapiens with `org.Hs.eg.db`.")
  }
  config
}

read_and_validate_csv_header <- function(path, input_name) {
  header <- read.csv(
    path, header = FALSE, nrows = 1L, check.names = FALSE,
    stringsAsFactors = FALSE, colClasses = "character"
  )
  header <- unname(unlist(header[1L, ], use.names = FALSE))
  if (anyNA(header) || any(trimws(header) == "")) {
    stop(input_name, " column headers must be non-blank.")
  }
  if (anyDuplicated(header)) {
    duplicates <- unique(header[duplicated(header)])
    stop(input_name, " column headers must be unique; duplicated: ", paste(duplicates, collapse = ", "))
  }
  invisible(header)
}

read_count_matrix <- function(path) {
  if (!file.exists(path)) stop("Count matrix not found: ", path)
  read_and_validate_csv_header(path, "Count matrix")
  tab <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"gene_id" %in% names(tab)) stop("Count matrix must contain a `gene_id` column.")
  if (ncol(tab) < 3L) stop("Count matrix must contain at least two sample columns.")
  if (anyNA(tab$gene_id) || any(tab$gene_id == "") || anyDuplicated(tab$gene_id)) {
    stop("`gene_id` values must be non-empty and unique.")
  }
  sample_columns <- setdiff(names(tab), "gene_id")
  non_numeric <- sample_columns[!vapply(tab[sample_columns], is.numeric, logical(1))]
  if (length(non_numeric)) stop("Non-numeric count columns: ", paste(non_numeric, collapse = ", "))
  counts <- as.matrix(tab[sample_columns])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- tab$gene_id
  if (anyNA(counts) || any(!is.finite(counts))) stop("Counts must be finite and non-missing.")
  if (any(counts < 0) || any(abs(counts - round(counts)) > 1e-8)) {
    stop("Counts must be non-negative integers.")
  }
  round(counts)
}

read_sample_metadata <- function(path) {
  if (!file.exists(path)) stop("Metadata not found: ", path)
  read_and_validate_csv_header(path, "Metadata")
  metadata <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"sample_id" %in% names(metadata)) stop("Metadata must contain a `sample_id` column.")
  if (anyNA(metadata$sample_id) || any(metadata$sample_id == "") || anyDuplicated(metadata$sample_id)) {
    stop("`sample_id` values must be non-empty and unique.")
  }
  rownames(metadata) <- metadata$sample_id
  metadata$sample_id <- NULL
  metadata
}

read_contrasts <- function(path) {
  if (!file.exists(path)) stop("Contrast table not found: ", path)
  read_and_validate_csv_header(path, "Contrast table")
  contrasts <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("label", "factor", "numerator", "denominator")
  missing <- setdiff(required, names(contrasts))
  if (length(missing)) stop("Contrast table missing columns: ", paste(missing, collapse = ", "))
  contrasts <- contrasts[required]
  if (!nrow(contrasts)) stop("Contrast table must contain at least one comparison.")
  if (anyNA(contrasts) || any(contrasts == "") || anyDuplicated(contrasts$label)) {
    stop("Contrast values must be complete and labels must be unique.")
  }
  if (any(contrasts$numerator == contrasts$denominator)) {
    labels <- contrasts$label[contrasts$numerator == contrasts$denominator]
    stop("Contrast numerator and denominator must differ; invalid labels: ", paste(labels, collapse = ", "))
  }
  stems <- safe_file_stem(contrasts$label)
  normalized_stems <- tolower(stems)
  if (anyDuplicated(normalized_stems)) {
    collisions <- unique(stems[duplicated(normalized_stems) | duplicated(normalized_stems, fromLast = TRUE)])
    stop("Contrast labels produce colliding output filenames after sanitization: ", paste(collisions, collapse = ", "))
  }
  contrasts
}

validate_analysis_inputs <- function(counts, metadata, contrasts, design) {
  if (!identical(colnames(counts), rownames(metadata))) {
    missing_metadata <- setdiff(colnames(counts), rownames(metadata))
    extra_metadata <- setdiff(rownames(metadata), colnames(counts))
    if (length(missing_metadata) || length(extra_metadata)) {
      stop(
        "Count and metadata sample IDs differ. Missing metadata: ",
        paste(missing_metadata, collapse = ", "), "; extra metadata: ",
        paste(extra_metadata, collapse = ", ")
      )
    }
    metadata <- metadata[colnames(counts), , drop = FALSE]
  }
  design_formula <- tryCatch(as.formula(design), error = function(e) NULL)
  if (is.null(design_formula)) stop("Invalid design formula: ", design)
  design_terms <- stats::terms(design_formula)
  if (attr(design_terms, "response") != 0L) stop("The design must not include a response variable.")
  term_labels <- attr(design_terms, "term.labels")
  if (!length(term_labels) || any(!grepl("^[A-Za-z.][A-Za-z0-9._]*$", term_labels))) {
    stop("Only additive categorical design terms are supported; interactions, transformations, and continuous covariates are not supported.")
  }
  design_variables <- all.vars(design_formula)
  absent <- setdiff(design_variables, names(metadata))
  if (length(absent)) stop("Design variables absent from metadata: ", paste(absent, collapse = ", "))
  if (anyNA(metadata[design_variables])) stop("Design variables must not contain missing values.")
  non_categorical <- design_variables[!vapply(
    metadata[design_variables],
    function(x) is.character(x) || is.factor(x) || is.logical(x),
    logical(1)
  )]
  if (length(non_categorical)) {
    stop("Design variables must be categorical for this workflow: ", paste(non_categorical, collapse = ", "))
  }
  for (column in design_variables) {
    metadata[[column]] <- factor(metadata[[column]])
    if (nlevels(metadata[[column]]) < 2L) stop("Design variable must contain at least two levels: ", column)
  }
  for (column in unique(contrasts$factor)) {
    if (!column %in% names(metadata)) stop("Contrast factor absent from metadata: ", column)
    if (!column %in% design_variables) stop("Contrast factor absent from design formula: ", column)
    levels_present <- unique(as.character(metadata[[column]]))
    rows <- contrasts$factor == column
    requested <- unique(c(contrasts$numerator[rows], contrasts$denominator[rows]))
    absent_levels <- setdiff(requested, levels_present)
    if (length(absent_levels)) {
      stop("Contrast levels absent from `", column, "`: ", paste(absent_levels, collapse = ", "))
    }
  }
  metadata
}

safe_file_stem <- function(x) {
  stem <- gsub("[^A-Za-z0-9._-]+", "_", x)
  stem <- gsub("^_+|_+$", "", stem)
  ifelse(nchar(stem), stem, "contrast")
}

prepare_output_directories <- function(output_dir) {
  paths <- c(
    root = output_dir,
    tables = file.path(output_dir, "tables"),
    figures = file.path(output_dir, "figures")
  )
  for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  paths
}

prepare_deseq_dataset <- function(counts, metadata, design, min_total_count) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Package `DESeq2` is required.")
  keep <- rowSums(counts) >= as.numeric(min_total_count)
  if (sum(keep) < 2L) stop("Fewer than two genes remain after low-count filtering.")
  DESeq2::DESeqDataSetFromMatrix(
    countData = counts[keep, , drop = FALSE],
    colData = metadata,
    design = as.formula(design)
  )
}

fit_deseq_model <- function(dds) {
  DESeq2::DESeq(dds, quiet = TRUE)
}

transform_for_visualization <- function(dds) {
  DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
}

run_deseq_contrast <- function(dds, contrast_row, fdr) {
  result <- DESeq2::results(
    dds,
    contrast = c(contrast_row$factor, contrast_row$numerator, contrast_row$denominator),
    alpha = fdr
  )
  tab <- as.data.frame(result)
  tab$gene_id <- rownames(tab)
  tab$significant <- !is.na(tab$padj) & tab$padj < fdr
  tab$direction <- ifelse(
    tab$significant & tab$log2FoldChange > 0, "up",
    ifelse(tab$significant & tab$log2FoldChange < 0, "down", "not_significant")
  )
  tab[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significant", "direction")]
}

make_pca_plot <- function(transformed, metadata, top_genes = 500L, plot_title = "PCA") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.")
  matrix <- SummarizedExperiment::assay(transformed)
  gene_variance <- apply(matrix, 1L, stats::var)
  selected <- names(sort(gene_variance, decreasing = TRUE))[seq_len(min(length(gene_variance), top_genes))]
  pca <- stats::prcomp(t(matrix[selected, , drop = FALSE]), center = TRUE, scale. = FALSE)
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  data <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L],
    metadata[rownames(pca$x), , drop = FALSE],
    check.names = FALSE
  )
  color_variable <- if ("condition" %in% names(data)) "condition" else names(metadata)[1L]
  ggplot2::ggplot(data, ggplot2::aes(x = PC1, y = PC2, color = .data[[color_variable]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = sample_id), vjust = -0.8, size = 3.5, check_overlap = TRUE, show.legend = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.18)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.16)) +
    ggplot2::labs(
      title = plot_title,
      x = sprintf("PC1 (%.1f%%)", variance[1L]),
      y = sprintf("PC2 (%.1f%%)", variance[2L]),
      color = color_variable
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.margin = ggplot2::margin(12, 16, 12, 16))
}

make_volcano_plot <- function(result_table, label, fdr, data_label = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.")
  plot_data <- result_table
  plot_data$minus_log10_padj <- -log10(pmax(plot_data$padj, .Machine$double.xmin))
  ggplot2::ggplot(plot_data, ggplot2::aes(x = log2FoldChange, y = minus_log10_padj, color = direction)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(up = "#B2182B", down = "#2166AC", not_significant = "#BDBDBD")) +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = "dashed", linewidth = 0.4) +
    ggplot2::labs(
      title = if (is.null(data_label)) label else paste(data_label, label, sep = ": "),
      x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

save_ggplot <- function(plot, path, width = 7, height = 5) {
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = 150)
  invisible(path)
}

save_correlation_heatmap <- function(transformed, path, plot_title = "Sample correlation") {
  matrix <- SummarizedExperiment::assay(transformed)
  correlation <- stats::cor(matrix, method = "pearson")
  grDevices::png(path, width = 1400, height = 1200, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  stats::heatmap(correlation, symm = TRUE, margins = c(9, 9), main = plot_title)
  invisible(path)
}

save_top_gene_heatmap <- function(transformed, metadata, path, top_genes = 30L, result_table = NULL, plot_title = "Top-ranked genes") {
  matrix <- SummarizedExperiment::assay(transformed)
  if (!is.null(result_table)) {
    ranking <- result_table$padj
    ranking[is.na(ranking)] <- Inf
    ranked_genes <- result_table$gene_id[order(ranking, -abs(result_table$log2FoldChange), na.last = TRUE)]
    selected <- intersect(ranked_genes, rownames(matrix))
    selected <- selected[seq_len(min(length(selected), top_genes))]
  } else {
    gene_variance <- apply(matrix, 1L, stats::var)
    selected <- names(sort(gene_variance, decreasing = TRUE))[seq_len(min(length(gene_variance), top_genes))]
  }
  display <- matrix[selected, , drop = FALSE]
  display <- t(scale(t(display)))
  display[!is.finite(display)] <- 0
  grDevices::png(path, width = 1500, height = 1600, res = 170)
  on.exit(grDevices::dev.off(), add = TRUE)
  stats::heatmap(display, scale = "none", margins = c(9, 7), main = plot_title)
  invisible(path)
}

run_optional_enrichment <- function(gene_ids, universe_ids, direction, contrast_label, config) {
  if (!identical(config$enrichment$species, "Homo sapiens") ||
      !identical(config$enrichment$organism_db, "org.Hs.eg.db")) {
    stop("Optional enrichment currently supports only Homo sapiens with `org.Hs.eg.db`.")
  }
  packages <- required_enrichment_packages()
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Enrichment packages missing: ", paste(missing, collapse = ", "))
  orgdb_name <- config$enrichment$organism_db
  if (!requireNamespace(orgdb_name, quietly = TRUE)) stop("Organism database missing: ", orgdb_name)
  orgdb <- getExportedValue(orgdb_name, orgdb_name)
  keytype <- config$enrichment$gene_id_type
  map_entrez <- function(ids) {
    AnnotationDbi::mapIds(
      orgdb, keys = unique(ids), column = "ENTREZID",
      keytype = keytype, multiVals = "first"
    )
  }
  selected_map <- map_entrez(gene_ids)
  universe_map <- map_entrez(universe_ids)
  genes_entrez <- unique(stats::na.omit(unname(selected_map)))
  universe_entrez <- unique(stats::na.omit(unname(universe_map)))
  coverage <- data.frame(
    contrast = contrast_label,
    direction = direction,
    selected_unique_source_ids = length(selected_map),
    selected_mapped_source_ids = sum(!is.na(selected_map)),
    selected_unique_entrez_ids = length(genes_entrez),
    universe_unique_source_ids = length(universe_map),
    universe_mapped_source_ids = sum(!is.na(universe_map)),
    universe_unique_entrez_ids = length(universe_entrez),
    stringsAsFactors = FALSE
  )
  if (!length(genes_entrez)) {
    return(list(go_bp = data.frame(), reactome = data.frame(), mapping_coverage = coverage))
  }
  go <- clusterProfiler::enrichGO(
    gene = genes_entrez, universe = universe_entrez, OrgDb = orgdb,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE
  )
  reactome <- ReactomePA::enrichPathway(
    gene = genes_entrez, universe = universe_entrez, organism = "human",
    pAdjustMethod = "BH", readable = TRUE
  )
  list(
    go_bp = as.data.frame(go),
    reactome = as.data.frame(reactome),
    mapping_coverage = coverage
  )
}
