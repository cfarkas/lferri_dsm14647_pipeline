#!/usr/bin/env bash
set -euo pipefail

# Download PacBio native subreads (bax.h5) from ENA using the run accessions.
#
# Default runs: ERR2028495–ERR2028504 (DSM14647 PacBio).
#
# Outputs:
#   DSM14647_PacBio_native/native_urls.txt
#   DSM14647_PacBio_native/*.bax.h5 ...

source "$(dirname "$0")/_common.sh"

OUTDIR="${OUTDIR:-$ROOT_DIR/DSM14647_PacBio_native}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

require curl
require aria2c
require awk

# You can override by exporting RUNS="ERR... ERR..."
RUNS=(${RUNS:-ERR2028504 ERR2028503 ERR2028502 ERR2028501 ERR2028500 ERR2028499 ERR2028498 ERR2028497 ERR2028496 ERR2028495})

log "Building ENA native URL list -> native_urls.txt"
: > native_urls.txt
for r in "${RUNS[@]}"; do
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${r}&result=read_run&fields=submitted_ftp&format=tsv" \
  | tail -n +2 \
  | cut -f2 \
  | tr ';' '\n' \
  | sed 's#^#ftp://#' >> native_urls.txt

  log "Added URLs for ${r}"
done

log "Total URLs: $(wc -l < native_urls.txt)"
log "Downloading with aria2c (tune -j/-s/-x if needed)"
aria2c -i native_urls.txt -c -s16 -x16 -j5 -d .

log "DONE: PacBio native files in $OUTDIR"
