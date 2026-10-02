#!/usr/bin/env sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
data_dir="$script_dir/data"
counts_file="$data_dir/GSE231397_raw_counts_GRCh38.p13_NCBI.tsv.gz"
soft_file="$data_dir/GSE231397_family.soft.gz"
prepared_counts="$data_dir/counts.csv"

counts_url='https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE231397&format=file&file=GSE231397_raw_counts_GRCh38.p13_NCBI.tsv.gz'
soft_url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231397/soft/GSE231397_family.soft.gz'
counts_sha256='877dee87d397d0abed6325252596c06ec64ac0fba2f00c0da58210f5659e3d85'
soft_sha256='e298fea9b1431b8acda9a1a105c17d17a8a5d61104a72b4d0cebfcc3a4295875'

mkdir -p "$data_dir"

download_if_missing() {
  destination=$1
  url=$2
  if [ ! -f "$destination" ]; then
    curl -L --fail --silent --show-error --output "$destination" "$url"
  fi
}

sha256_file() {
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  elif command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    echo "A SHA-256 utility (shasum or sha256sum) is required." >&2
    exit 1
  fi
}

verify_file() {
  path=$1
  expected=$2
  observed=$(sha256_file "$path")
  if [ "$observed" != "$expected" ]; then
    echo "SHA-256 mismatch for $path" >&2
    echo "expected: $expected" >&2
    echo "observed: $observed" >&2
    exit 1
  fi
}

download_if_missing "$counts_file" "$counts_url"
download_if_missing "$soft_file" "$soft_url"
verify_file "$counts_file" "$counts_sha256"
verify_file "$soft_file" "$soft_sha256"

counts_bytes=$(wc -c < "$counts_file" | tr -d ' ')
if [ "$counts_bytes" != "576699" ]; then
  echo "Unexpected count archive size: $counts_bytes bytes" >&2
  exit 1
fi

Rscript --vanilla "$script_dir/prepare_inputs.R" \
  "$counts_file" \
  "$soft_file" \
  "$script_dir/metadata.csv" \
  "$prepared_counts"

echo "Prepared six-sample count matrix: $prepared_counts"
