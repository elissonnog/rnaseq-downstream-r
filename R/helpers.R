required_core_packages <- function() {
  c("rmarkdown", "knitr", "DESeq2", "SummarizedExperiment", "ggplot2", "openssl")
}

required_enrichment_packages <- function() {
  c("AnnotationDbi", "clusterProfiler", "ReactomePA", "org.Hs.eg.db")
}

is_pandoc_available <- function() {
  nzchar(Sys.which("pandoc")) ||
    (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available())
}

check_dependencies <- function(config = NULL, include_enrichment = FALSE) {
  packages <- required_core_packages()
  required_for <- rep("core analysis", length(packages))
  if (isTRUE(include_enrichment)) {
    optional <- setdiff(required_enrichment_packages(), packages)
    packages <- c(packages, optional)
    required_for <- c(required_for, rep("optional enrichment", length(optional)))
  }
  if (!is.null(config) && isTRUE(config$gene_labels$enabled)) {
    label_packages <- setdiff(c("AnnotationDbi", config$gene_labels$organism_db), packages)
    packages <- c(packages, label_packages)
    required_for <- c(required_for, rep("configured gene labels", length(label_packages)))
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
    "pca_top_genes", "heatmap_top_genes", "gene_labels", "enrichment"
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
  label_required <- c("enabled", "gene_id_type", "organism_db", "column")
  label_missing <- setdiff(label_required, names(config$gene_labels))
  if (length(label_missing)) {
    stop("Missing gene-label configuration fields: ", paste(label_missing, collapse = ", "))
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
    values <- unique(as.character(metadata[[column]]))
    if (length(values) && all(grepl("^[+-]?[0-9.]+$", values))) {
      stop("Design variables encoded only as numeric-looking values are rejected; use explicit categorical labels: ", column)
    }
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

classify_result_status <- function(input_row_sums, retained, pvalue, padj, significant) {
  status <- rep("unavailable", length(input_row_sums))
  status[input_row_sums == 0] <- "all_zero"
  status[input_row_sums > 0 & !retained] <- "low_count_filtered"
  status[retained & is.finite(pvalue) & !is.finite(padj)] <- "independent_filtered"
  status[retained & is.finite(padj) & !significant] <- "tested_nonsignificant"
  status[retained & is.finite(padj) & significant] <- "significant"
  status
}

complete_result_ledger <- function(result_table, input_counts, id_column, min_total_count,
                                   increase_label, decrease_label) {
  input_ids <- rownames(input_counts)
  ledger <- data.frame(
    id = input_ids,
    input_total_count = rowSums(input_counts),
    stringsAsFactors = FALSE
  )
  names(ledger)[1L] <- id_column
  result_table$retained_for_model <- TRUE
  ledger <- merge(ledger, result_table, by = id_column, all.x = TRUE, sort = FALSE)
  ledger <- ledger[match(input_ids, ledger[[id_column]]), , drop = FALSE]
  retained <- !is.na(ledger$retained_for_model)
  significant <- retained & is.finite(ledger$padj) & ledger$padj < attr(result_table, "fdr")
  ledger$test_status <- classify_result_status(
    ledger$input_total_count, retained, ledger$pvalue, ledger$padj, significant
  )
  ledger$significant <- ifelse(
    ledger$test_status %in% c("significant", "tested_nonsignificant"),
    ledger$test_status == "significant", NA
  )
  ledger$direction <- NA_character_
  ledger$direction[ledger$test_status == "significant" & ledger$log2FoldChange > 0] <- increase_label
  ledger$direction[ledger$test_status == "significant" & ledger$log2FoldChange < 0] <- decrease_label
  ledger$retained_for_model <- retained
  ledger
}

fit_deseq_model <- function(dds) {
  DESeq2::DESeq(dds, quiet = TRUE)
}

transform_for_visualization <- function(dds) {
  DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
}

run_deseq_contrast <- function(dds, contrast_row, fdr, input_counts = NULL,
                               min_total_count = 0) {
  result <- DESeq2::results(
    dds,
    contrast = c(contrast_row$factor, contrast_row$numerator, contrast_row$denominator),
    alpha = fdr,
    independentFiltering = TRUE
  )
  tab <- as.data.frame(result)
  tab$gene_id <- rownames(tab)
  attr(tab, "fdr") <- fdr
  settings <- list(
    alpha = fdr,
    p_adjust_method = "BH",
    independent_filtering = TRUE,
    filter_threshold = unname(S4Vectors::metadata(result)$filterThreshold %||% NA_real_),
    bh_scope = "one DESeq2 contrast among retained rows passing the selected independent-filter threshold with finite raw p-values"
  )
  if (is.null(input_counts)) {
    tab$input_total_count <- NA_real_
    tab$retained_for_model <- TRUE
    tab$test_status <- classify_result_status(
      rep(1, nrow(tab)), rep(TRUE, nrow(tab)), tab$pvalue, tab$padj,
      is.finite(tab$padj) & tab$padj < fdr
    )
    tab$significant <- ifelse(
      tab$test_status %in% c("significant", "tested_nonsignificant"),
      tab$test_status == "significant", NA
    )
    tab$direction <- NA_character_
    tab$direction[tab$test_status == "significant" & tab$log2FoldChange > 0] <- "up"
    tab$direction[tab$test_status == "significant" & tab$log2FoldChange < 0] <- "down"
  } else {
    tab <- complete_result_ledger(tab, input_counts, "gene_id", min_total_count, "up", "down")
  }
  tab <- tab[, c("gene_id", "input_total_count", "retained_for_model", "baseMean", "log2FoldChange",
                 "lfcSE", "stat", "pvalue", "padj", "test_status", "significant", "direction")]
  attr(tab, "results_settings") <- settings
  tab
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

summarize_contrast_result <- function(result_table, contrast_label) {
  settings <- attr(result_table, "results_settings")
  data.frame(
    contrast = contrast_label,
    input_rows = nrow(result_table),
    retained_rows = sum(result_table$retained_for_model),
    finite_pvalue = sum(is.finite(result_table$pvalue)),
    finite_padj = sum(is.finite(result_table$padj)),
    significant = sum(result_table$test_status == "significant"),
    increased = sum(result_table$direction == "up", na.rm = TRUE),
    decreased = sum(result_table$direction == "down", na.rm = TRUE),
    all_zero = sum(result_table$test_status == "all_zero"),
    low_count_filtered = sum(result_table$test_status == "low_count_filtered"),
    unavailable = sum(result_table$test_status == "unavailable"),
    independent_filtered = sum(result_table$test_status == "independent_filtered"),
    tested_nonsignificant = sum(result_table$test_status == "tested_nonsignificant"),
    alpha = settings$alpha,
    p_adjust_method = settings$p_adjust_method,
    independent_filtering = settings$independent_filtering,
    filter_threshold = settings$filter_threshold,
    bh_scope = settings$bh_scope,
    stringsAsFactors = FALSE
  )
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
  plot_data$result_class <- ifelse(plot_data$test_status == "significant", plot_data$direction,
                                   plot_data$test_status)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = log2FoldChange, y = minus_log10_padj, color = result_class)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(
      up = "#B2182B", down = "#2166AC", tested_nonsignificant = "#BDBDBD",
      independent_filtered = "#7F7F7F", unavailable = "#4D4D4D"
    ), na.value = "#4D4D4D") +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = "dashed", linewidth = 0.4) +
    ggplot2::labs(
      title = if (is.null(data_label)) label else paste(data_label, label, sep = ": "),
      x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

make_ma_plot <- function(result_table, label, fdr, data_label = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.")
  plot_data <- result_table[
    is.finite(result_table$baseMean) & result_table$baseMean > 0 &
      is.finite(result_table$log2FoldChange),
    , drop = FALSE
  ]
  plot_data$result_class <- ifelse(plot_data$test_status == "significant", plot_data$direction,
                                   plot_data$test_status)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = baseMean, y = log2FoldChange, color = result_class)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.35, na.rm = TRUE) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c(
      up = "#B2182B", down = "#2166AC", tested_nonsignificant = "#BDBDBD",
      independent_filtered = "#7F7F7F", unavailable = "#4D4D4D"
    ), na.value = "#4D4D4D") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::labs(
      title = if (is.null(data_label)) label else paste(data_label, label, sep = ": "),
      subtitle = sprintf("Significance: adjusted p-value < %.3g; fold changes are unshrunk", fdr),
      x = "Mean normalized count", y = "log2 fold change", color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

select_significant_genes <- function(result_table, fdr, max_genes) {
  qualifying <- !is.na(result_table$padj) & result_table$padj < fdr &
    is.finite(result_table$log2FoldChange)
  candidates <- result_table[qualifying, , drop = FALSE]
  if (!nrow(candidates)) return(character())
  ordering <- order(
    candidates$padj,
    -abs(candidates$log2FoldChange),
    as.character(candidates$gene_id),
    na.last = TRUE
  )
  candidates$gene_id[ordering][seq_len(min(nrow(candidates), max_genes))]
}

map_gene_display_labels <- function(gene_ids, config) {
  labels <- as.character(gene_ids)
  mapping <- config$gene_labels
  if (!isTRUE(mapping$enabled) || identical(mapping$gene_id_type, mapping$column)) return(labels)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace(mapping$organism_db, quietly = TRUE)) {
    stop("Configured gene-label mapping requires AnnotationDbi and ", mapping$organism_db, ".")
  }
  orgdb <- getExportedValue(mapping$organism_db, mapping$organism_db)
  if (!mapping$gene_id_type %in% AnnotationDbi::keytypes(orgdb)) {
    stop("Configured gene-label key type is unavailable: ", mapping$gene_id_type)
  }
  symbols <- suppressMessages(AnnotationDbi::mapIds(
    orgdb,
    keys = labels,
    column = mapping$column,
    keytype = mapping$gene_id_type,
    multiVals = "first"
  ))
  symbols <- unname(symbols[labels])
  has_symbol <- !is.na(symbols) & nzchar(symbols)
  labels[has_symbol] <- paste0(symbols[has_symbol], " (", labels[has_symbol], ")")
  labels
}

save_significant_gene_heatmap <- function(
    transformed, metadata, path, top_genes, result_table, fdr, config,
    plot_title = "Significant genes across replicates") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.")
  matrix <- SummarizedExperiment::assay(transformed)
  selected <- select_significant_genes(result_table, fdr, top_genes)
  selected <- intersect(selected, rownames(matrix))
  if (!length(selected)) {
    empty_plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = sprintf("No genes passed adjusted p-value < %.3g", fdr), size = 5) +
      ggplot2::labs(title = plot_title) +
      ggplot2::theme_void(base_size = 11)
    save_ggplot(empty_plot, path, width = 8, height = 4)
    return(invisible(list(selected_genes = character(), available_significant = 0L)))
  }

  significant_available <- sum(!is.na(result_table$padj) & result_table$padj < fdr)
  display <- matrix[selected, , drop = FALSE]
  display <- sweep(display, 1L, rowMeans(display), FUN = "-")
  labels <- map_gene_display_labels(selected, config)
  rownames(display) <- labels
  sample_order <- rownames(metadata)
  display <- display[, sample_order, drop = FALSE]

  heatmap_data <- as.data.frame(as.table(display), stringsAsFactors = FALSE)
  names(heatmap_data) <- c("gene", "sample_id", "row_centered_vst")
  heatmap_data$sample_id <- factor(heatmap_data$sample_id, levels = sample_order)
  heatmap_levels <- c(rev(labels), "Condition")
  heatmap_data$gene <- factor(heatmap_data$gene, levels = heatmap_levels)
  condition <- if ("condition" %in% names(metadata)) as.character(metadata[sample_order, "condition"]) else rep("sample", length(sample_order))
  annotation <- data.frame(
    sample_id = factor(sample_order, levels = sample_order),
    gene = factor(rep("Condition", length(sample_order)), levels = heatmap_levels),
    condition = condition,
    stringsAsFactors = FALSE
  )

  plot <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = sample_id, y = gene, fill = row_centered_vst)) +
    ggplot2::geom_tile() +
    ggplot2::geom_point(
      data = annotation,
      ggplot2::aes(x = sample_id, y = gene, color = condition),
      inherit.aes = FALSE,
      shape = 15, size = 5
    ) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
    ggplot2::labs(
      title = plot_title,
      subtitle = sprintf("%d of %d genes passing adjusted p-value < %.3g; deterministic ranking by adjusted p-value and effect size", length(selected), significant_available, fdr),
      x = NULL, y = NULL, fill = "Row-centered\nVST", color = "Condition"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.margin = ggplot2::margin(18, 12, 12, 12)
    )
  save_ggplot(plot, path, width = 9, height = max(6, 0.28 * length(selected) + 2.5))
  invisible(list(selected_genes = selected, available_significant = significant_available))
}

save_ggplot <- function(plot, path, width = 7, height = 5) {
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = 150)
  invisible(path)
}

save_dispersion_plot <- function(dds, path, plot_title = "DESeq2 dispersion estimates") {
  grDevices::png(path, width = 1200, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  DESeq2::plotDispEsts(dds, main = plot_title)
  invisible(path)
}

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  as.character(openssl::sha256(file(path)))
}

input_manifest <- function(paths) {
  data.frame(
    input = names(paths),
    path = normalizePath(unname(paths), mustWork = TRUE),
    sha256 = vapply(unname(paths), sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
}

software_versions <- function(packages) {
  data.frame(
    software = c("R", packages, "Pandoc"),
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      vapply(packages, function(package) {
        if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else NA_character_
      }, character(1)),
      if (is_pandoc_available()) as.character(rmarkdown::pandoc_version()) else NA_character_
    ),
    stringsAsFactors = FALSE
  )
}

code_manifest <- function(paths) {
  data.frame(
    file = unname(paths), sha256 = vapply(paths, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
}

save_correlation_heatmap <- function(transformed, path, plot_title = "Sample correlation") {
  matrix <- SummarizedExperiment::assay(transformed)
  correlation <- stats::cor(matrix, method = "pearson")
  plot_data <- as.data.frame(as.table(correlation), stringsAsFactors = FALSE)
  names(plot_data) <- c("sample_x", "sample_y", "pearson_correlation")
  plot_data$sample_x <- factor(plot_data$sample_x, levels = colnames(correlation))
  plot_data$sample_y <- factor(plot_data$sample_y, levels = rev(rownames(correlation)))
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(sample_x, sample_y, fill = pearson_correlation)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", pearson_correlation)), size = 3) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(-1, 1)) +
    ggplot2::labs(title = plot_title, x = NULL, y = NULL, fill = "Pearson r") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_ggplot(plot, path, width = 7, height = 6)
  invisible(correlation)
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

empty_enrichment_table <- function() {
  data.frame(
    ID = character(), Description = character(), GeneRatio = character(), BgRatio = character(),
    pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(), geneID = character(),
    Count = integer(), stringsAsFactors = FALSE
  )
}

map_enrichment_ids <- function(gene_ids, universe_ids, map_function) {
  selected_source <- unique(as.character(gene_ids[!is.na(gene_ids) & nzchar(as.character(gene_ids))]))
  universe_source <- unique(as.character(universe_ids[!is.na(universe_ids) & nzchar(as.character(universe_ids))]))
  selected_map <- if (length(selected_source)) map_function(selected_source) else setNames(character(), character())
  universe_map <- if (length(universe_source)) map_function(universe_source) else setNames(character(), character())
  selected_entrez <- unique(as.character(stats::na.omit(unname(selected_map))))
  universe_entrez <- unique(as.character(stats::na.omit(unname(universe_map))))
  status <- if (!length(universe_source)) {
    "empty_universe"
  } else if (!length(universe_entrez)) {
    "universe_unmapped"
  } else if (!length(selected_source)) {
    "no_selected_genes"
  } else if (!length(selected_entrez)) {
    "selected_unmapped"
  } else if (length(selected_entrez) < length(selected_source)) {
    "partial_selected_mapping"
  } else {
    "ready"
  }
  list(
    selected_source = selected_source, universe_source = universe_source,
    selected_entrez = selected_entrez, universe_entrez = universe_entrez, status = status,
    selected_mapped = sum(!is.na(selected_map) & nzchar(as.character(selected_map))),
    universe_mapped = sum(!is.na(universe_map) & nzchar(as.character(universe_map)))
  )
}

run_optional_enrichment <- function(gene_ids, universe_ids, direction, contrast_label, config) {
  if (!identical(config$enrichment$species, "Homo sapiens") ||
      !identical(config$enrichment$organism_db, "org.Hs.eg.db")) {
    stop("Optional enrichment currently supports only Homo sapiens with `org.Hs.eg.db`.")
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Enrichment requires AnnotationDbi.")
  orgdb_name <- config$enrichment$organism_db
  if (!requireNamespace(orgdb_name, quietly = TRUE)) stop("Organism database missing: ", orgdb_name)
  orgdb <- getExportedValue(orgdb_name, orgdb_name)
  keytype <- config$enrichment$gene_id_type
  if (!keytype %in% AnnotationDbi::keytypes(orgdb)) {
    stop("Configured enrichment key type is unavailable in ", orgdb_name, ": ", keytype)
  }
  map_entrez <- function(ids) {
    ids <- unique(as.character(ids))
    tryCatch(
      suppressMessages(AnnotationDbi::mapIds(
        orgdb, keys = ids, column = "ENTREZID",
        keytype = keytype, multiVals = "first"
      )),
      error = function(error) setNames(rep(NA_character_, length(ids)), ids)
    )
  }
  mapped <- map_enrichment_ids(gene_ids, universe_ids, map_entrez)
  coverage <- data.frame(
    contrast = contrast_label,
    direction = direction,
    status = mapped$status,
    source_keytype = keytype,
    selected_unique_source_ids = length(mapped$selected_source),
    selected_mapped_source_ids = mapped$selected_mapped,
    selected_unique_entrez_ids = length(mapped$selected_entrez),
    universe_unique_source_ids = length(mapped$universe_source),
    universe_mapped_source_ids = mapped$universe_mapped,
    universe_unique_entrez_ids = length(mapped$universe_entrez),
    stringsAsFactors = FALSE
  )
  if (!mapped$status %in% c("ready", "partial_selected_mapping")) {
    return(list(go_bp = empty_enrichment_table(), reactome = empty_enrichment_table(), mapping_coverage = coverage))
  }
  analysis_packages <- c("clusterProfiler", "ReactomePA")
  missing <- analysis_packages[!vapply(analysis_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Enrichment analysis packages missing: ", paste(missing, collapse = ", "))
  go <- clusterProfiler::enrichGO(
    gene = mapped$selected_entrez, universe = mapped$universe_entrez, OrgDb = orgdb,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE
  )
  reactome <- ReactomePA::enrichPathway(
    gene = mapped$selected_entrez, universe = mapped$universe_entrez, organism = "human",
    pAdjustMethod = "BH", readable = TRUE
  )
  list(
    go_bp = as.data.frame(go),
    reactome = as.data.frame(reactome),
    mapping_coverage = coverage
  )
}
