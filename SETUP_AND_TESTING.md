# Setup, custom data, and testing

Run every command from the repository root. Configuration paths may be absolute; relative paths are resolved from the repository root.

## Install requirements

Install a current R release and external Pandoc. Confirm both executables are on `PATH`:

```sh
R --version
pandoc --version
```

Pandoc installation instructions and binaries are available from the [official Pandoc site](https://pandoc.org/installing.html). Common package-manager commands are `brew install pandoc` on macOS and `sudo apt install pandoc` on Debian/Ubuntu.

Install the core R packages:

```r
install.packages(c("BiocManager", "rmarkdown", "knitr", "ggplot2", "openssl"))
BiocManager::install(c("DESeq2", "SummarizedExperiment"))
```

Gene-label mapping in the public configuration additionally requires:

```r
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))
```

Optional human enrichment requires `AnnotationDbi`, `org.Hs.eg.db`, `clusterProfiler`, and `ReactomePA`:

```r
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "clusterProfiler", "ReactomePA"))
```

No workflow command installs packages automatically. The last verified environment used R 4.6.1, Pandoc 3.12, rmarkdown 2.32, knitr 1.52, ggplot2 4.0.3, openssl 2.4.2, DESeq2 1.52.0, SummarizedExperiment 1.42.0, AnnotationDbi 1.74.0, and org.Hs.eg.db 3.23.1. This is a tested-version record, not a lockfile.

## Custom-data contract

Use three CSV files:

- Counts: the first column must be `gene_id`; every remaining column is a sample. Gene IDs must be unique and non-empty, and headers must be unique and nonblank. Counts must be finite, nonmissing, nonnegative integers, with at least two sample columns. Do not include annotation columns because they would be interpreted as samples.
- Metadata: must contain a unique, non-empty `sample_id` column. Its values must exactly match the count-column names; order may differ and is realigned safely. Every design variable must be complete and categorical, contain at least two levels, and use descriptive labels rather than numeric or numeric-looking codes.
- Contrasts: must contain `label,factor,numerator,denominator`. Labels must be unique and must not collide after filename sanitization. `factor` must be a variable in the design, numerator and denominator must be different, and both requested levels must exist in metadata.

Only additive categorical designs such as `~ condition` or `~ batch + condition` are supported. Interactions, transformations, continuous covariates, and response variables are rejected.

The synthetic files are small schema examples:

- [counts.csv](example/synthetic/counts.csv)
- [metadata.csv](example/synthetic/metadata.csv)
- [contrasts.csv](example/synthetic/contrasts.csv)

## Configure and run custom data

1. Copy the complete example configuration; do not build a partial list.

```sh
cp config/example_config.R config/my_analysis.R
```

2. Edit `config/my_analysis.R`:

   - Set `project_title`, `data_label`, the three input paths, and a new `output_dir`.
   - Set the additive `design`, `min_total_count`, `fdr`, and plotting limits.
   - For gene labels, either keep `gene_labels$enabled = FALSE`, or declare the input `gene_id_type`, installed `organism_db`, and a valid output `column`.
   - Keep enrichment disabled unless the data are human, identifiers match `enrichment$gene_id_type`, and all optional packages above are installed. The implemented enrichment path currently supports `species = "Homo sapiens"` with `org.Hs.eg.db`.

3. Check the exact configuration. A missing required dependency produces a nonzero exit and nothing is installed:

```sh
Rscript --vanilla run_analysis.R --check-dependencies config/my_analysis.R
```

4. Render the analysis:

```sh
Rscript --vanilla run_analysis.R config/my_analysis.R
```

The configured output directory receives the HTML report, complete result ledgers, normalized counts, diagnostics, input/code hashes, resolved configuration, software versions, warnings, and session information.

## Tests and public-example regeneration

```sh
Rscript --vanilla tests/test_helpers.R
Rscript --vanilla tests/test_core_analysis.R
Rscript --vanilla tests/test_preflight.R
Rscript --vanilla tests/test_render.R
```

To regenerate the public report and compact portfolio files:

```sh
Rscript --vanilla run_analysis.R config/gse231397_e2_vs_vehicle.R
Rscript --vanilla scripts/export_portfolio_results.R
Rscript --vanilla scripts/render_portfolio_snapshot.R
```

Empty and wholly unmapped enrichment inputs are integration-tested without `clusterProfiler` or `ReactomePA`; partial mapping is unit-tested. A live enrichment run still requires all optional packages.
