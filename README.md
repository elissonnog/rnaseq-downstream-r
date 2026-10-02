# Reproducible bulk RNA-seq downstream analysis in R

A configuration-driven DESeq2 workflow for turning a gene-by-sample count matrix into an auditable downstream analysis. It connects sample-level quality assessment, normalization, declared differential-expression contrasts, publication-ready figures, and optional human functional enrichment in one reproducible HTML report.

This workflow model was developed during my postdoctoral research at Van Andel Institute. The reusable version in this repository is adapted from my original analysis code and refactored into a portable report, validation helpers, and a synthetic test fixture. The demonstration data and figures below are entirely synthetic and are not Van Andel Institute experimental results.

## Questions this workflow answers

- Do samples separate by the biological groups declared in the metadata, and are any samples discordant in PCA or sample correlation?
- Which genes change for each explicitly declared numerator-versus-denominator comparison?
- How large and statistically supported are those changes after Benjamini-Hochberg correction?
- Which top-ranked genes drive the expression pattern across samples?
- For human data, which GO Biological Process and Reactome terms are over-represented among significant up- and down-regulated genes when enrichment is enabled?

The workflow does not perform read alignment, transcript quantification, or raw-read QC. Those are upstream steps; this repository begins with integer gene counts.

## Synthetic demonstration

![PCA of the synthetic demonstration data](figures/synthetic-pca.png)

*Synthetic demonstration only. PCA of the 30 most variable variance-stabilized genes shows the planted condition-level separation used to exercise sample-level QC. It is not a biological result.*

![Volcano plot of the synthetic demonstration contrast](figures/synthetic-volcano.png)

*Synthetic demonstration only. The volcano plot summarizes the planted `treatment_vs_control` contrast by log2 fold change and adjusted p-value. It verifies model, classification, and plotting behavior and has no biological interpretation.*

## Implemented analysis

| Stage | What the workflow does | Main deliverable |
| --- | --- | --- |
| Input validation | Checks unique headers, integer counts, sample matching, categorical additive designs, declared contrast levels, and safe output names | Early, readable failures before model fitting |
| Filtering and normalization | Removes genes below the configured total-count threshold and estimates DESeq2 size factors and dispersions | Normalized count matrix and run metadata |
| Sample-level QC | Applies a variance-stabilizing transformation and evaluates PCA, library size, detected genes, and sample correlation | PCA and correlation heatmap |
| Differential expression | Fits the configured design and evaluates every row of the contrast table | Complete DESeq2 table per comparison |
| Result visualization | Classifies genes at the configured FDR and plots effect size, significance, and top-ranked expression patterns | Volcano plot and heatmap per comparison |
| Optional enrichment | Separately tests significant up- and down-regulated human genes against the tested-gene universe | GO-BP, Reactome, and ID-mapping tables |

## Requirements

- R
- Pandoc
- CRAN: `rmarkdown`, `knitr`, `ggplot2`
- Bioconductor: `DESeq2` and its dependency `SummarizedExperiment`

Optional human functional enrichment additionally uses `AnnotationDbi`, `clusterProfiler`, `ReactomePA`, and `org.Hs.eg.db`.

## Inputs and configuration

- `counts.csv`: `gene_id` followed by non-negative integer sample columns.
- `metadata.csv`: one row per sample, with `sample_id` and categorical design variables.
- `contrasts.csv`: `label`, `factor`, `numerator`, and `denominator`.
- An R configuration file modeled on `config/example_config.R`, defining paths, design, filtering threshold, FDR, plot sizes, and optional enrichment settings.

For a study with condition and batch effects, the central configuration would look like:

```r
analysis_config <- list(
  project_title = "Study name",
  data_label = "Internal study data",
  counts_path = "data/counts.csv",
  metadata_path = "data/metadata.csv",
  contrasts_path = "data/contrasts.csv",
  output_dir = "output/study_name",
  design = "~ batch + condition",
  min_total_count = 10L,
  fdr = 0.05,
  pca_top_genes = 500L,
  heatmap_top_genes = 30L,
  enrichment = list(
    enabled = FALSE,
    species = "Homo sapiens",
    gene_id_type = "SYMBOL",
    organism_db = "org.Hs.eg.db"
  )
)
```

Comparisons are never inferred from column order. Each requested comparison is declared in `contrasts.csv`, for example:

```csv
label,factor,numerator,denominator
treatment_vs_control,condition,treatment,control
```

## Run the synthetic example

```sh
git clone https://github.com/elissonnog/rnaseq-downstream-r.git
cd rnaseq-downstream-r
Rscript tests/test_helpers.R
Rscript run_analysis.R config/example_config.R
```

The synthetic fixture can be regenerated deterministically with:

```sh
Rscript example/synthetic/generate_fixture.R
```

To inspect dependency availability without running the analysis:

```sh
Rscript run_analysis.R --check-dependencies
```

## Outputs

The configured output directory receives:

- `rnaseq_downstream.html`
- `tables/run_metadata.csv`
- `tables/normalized_counts.csv`
- `tables/<contrast>_deseq2.csv`
- `figures/pca.png`
- `figures/sample_correlation.png`
- `figures/<contrast>_volcano.png`
- `figures/<contrast>_top_genes_heatmap.png`
- enrichment tables and identifier-mapping summaries when enrichment is enabled

Each report records the configured design, thresholds, declared comparisons, data-status label, and R session information. A real worked example should be added only when its data provenance and release permissions are confirmed.

No license has been selected.
