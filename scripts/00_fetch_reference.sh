#!/usr/bin/env bash
set -euo pipefail

# Fetch original reference genome (GCF_000755505.1) from NCBI
# Outputs:
#   ref/lferri.fa
#   ref/lferri.gtf

source "$(dirname "$0")/_common.sh"

ACC="${ACC:-GCF_000755505.1}"
OUT_REF_DIR="$ROOT_DIR/ref"
mkdir -p "$OUT_REF_DIR"

require datasets
require unzip
require gffread

log "Downloading reference $ACC via NCBI datasets..."
work="$ROOT_DIR/ref/ncbi_dataset_tmp"
rm -rf "$work"; mkdir -p "$work"

(
  cd "$work"
  datasets download genome accession "$ACC" --include genome,gff3,gtf
  unzip -o ncbi_dataset.zip
)

# Copy genome FASTA
fa_src=$(ls "$work"/ncbi_dataset/data/$ACC/*_genomic.fna | head -n 1 || true)
if [[ -z "${fa_src}" || ! -s "${fa_src}" ]]; then
  echo "ERROR: could not find *_genomic.fna in datasets download" >&2
  exit 2
fi
cp -f "$fa_src" "$OUT_REF_DIR/lferri.fa"

# Copy/convert GTF
if [[ -f "$work"/ncbi_dataset/data/$ACC/genomic.gtf ]]; then
  cp -f "$work"/ncbi_dataset/data/$ACC/genomic.gtf "$OUT_REF_DIR/lferri.gtf"
else
  gff_src=$(ls "$work"/ncbi_dataset/data/$ACC/*_genomic.gff | head -n 1 || true)
  if [[ -z "${gff_src}" || ! -s "${gff_src}" ]]; then
    echo "ERROR: could not find *_genomic.gff for conversion to GTF" >&2
    exit 3
  fi
  gffread "$gff_src" -T -o "$OUT_REF_DIR/lferri.gtf"
fi

log "Wrote: $OUT_REF_DIR/lferri.fa"
log "Wrote: $OUT_REF_DIR/lferri.gtf"
