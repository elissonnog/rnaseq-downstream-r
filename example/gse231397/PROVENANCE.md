# GSE231397 public-study example provenance

This example reanalyzes six public MCF-7 RNA-seq samples from GEO accession [GSE231397](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231397): three DMSO vehicle biological replicates and three biological replicates treated with 1 nM 17beta-estradiol (E2) for 24 hours. The comparison is `E2_vs_vehicle` with design `~ condition`. Replicate numbers are descriptive labels only; the analysis does not model or infer pairing.

## Official source files

- NCBI-generated raw gene counts: `GSE231397_raw_counts_GRCh38.p13_NCBI.tsv.gz`
  - URL: <https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE231397&format=file&file=GSE231397_raw_counts_GRCh38.p13_NCBI.tsv.gz>
  - SHA-256: `877dee87d397d0abed6325252596c06ec64ac0fba2f00c0da58210f5659e3d85`
  - Size: 576,699 bytes
  - Verified content: 39,376 unique GeneIDs by 12 sample columns; all counts are non-negative integers.
- GEO family SOFT metadata: `GSE231397_family.soft.gz`
  - URL: <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231397/soft/GSE231397_family.soft.gz>
  - SHA-256: `e298fea9b1431b8acda9a1a105c17d17a8a5d61104a72b4d0cebfcc3a4295875`

The downloader verifies these checksums and the expected count-archive size before preparing the six-sample matrix. The count archive and prepared matrix are intentionally ignored by Git and are not redistributed from this repository.

## Sample verification

The SOFT titles and treatment fields identify vehicle samples `GSM7270183`, `GSM7270188`, and `GSM7270193`, and estradiol samples `GSM7270186`, `GSM7270190`, and `GSM7270195`. Each group contains three biological replicates.

NCBI's standardized RNA-seq count resource reports HISAT2 alignment and featureCounts quantification against GRCh38.p13, NCBI annotation release 109.20190905. See [NCBI RNA-seq counts](https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html). This repository starts from those NCBI-generated counts; it does not reproduce raw-read processing or the exact analysis pipeline used by the study authors.

Study citation: Min et al., *Proceedings of the National Academy of Sciences* (2024), [doi:10.1073/pnas.2321344121](https://doi.org/10.1073/pnas.2321344121).

The upstream GEO data are not presented here as CC0. Their terms and attribution remain separate from this repository, and no repository license grant applies to the downloaded data.
