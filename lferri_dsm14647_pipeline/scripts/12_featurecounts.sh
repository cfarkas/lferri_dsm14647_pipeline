#!/usr/bin/env bash
set -euo pipefail

# Generate gene-level counts using featureCounts.
#
# Inputs:
#   L_ferri_bakta/L_ferriphilum.fixed.gtf
#   align/*.bam
# Outputs:
#   counts/gene_counts.txt
#   counts/gene_counts.matrix.tsv

source "$(dirname "$0")/_common.sh"

THREADS="${THREADS:-32}"
GTF="${GTF:-$ROOT_DIR/L_ferri_bakta/L_ferriphilum.fixed.gtf}"
OUTDIR="${OUTDIR:-$ROOT_DIR/counts}"
BAM_GLOB="${BAM_GLOB:-$ROOT_DIR/align/*.bam}"

require featureCounts

mkdir -p "$OUTDIR"

if [[ ! -s "$GTF" ]]; then
  log "Missing GTF: $GTF"
  log "Run: scripts/09_format_gtf.sh"
  exit 1
fi

# Check BAMs
shopt -s nullglob
bams=( $BAM_GLOB )
shopt -u nullglob

if [[ ${#bams[@]} -eq 0 ]]; then
  log "No BAMs found matching: $BAM_GLOB"
  log "Run: scripts/11_align_rnaseq.sh"
  exit 2
fi

log "Running featureCounts on ${#bams[@]} BAMs"
featureCounts \
  -T "$THREADS" \
  -p -B -C \
  -s 0 \
  -a "$GTF" \
  -t CDS \
  -g gene_id \
  -o "$OUTDIR/gene_counts.txt" \
  "${bams[@]}"

# Clean the counts matrix to a simple TSV
awk 'BEGIN{OFS="\t"}
     NR==2{
        h="GeneID";
        for(i=7;i<=NF;i++){
          gsub(/^.*\//,"",$i);
          gsub(/\.bam$/,"",$i);
          h=h"\t"$i;
        }
        print h;
        next
     }
     NR>2{
        l=$1;
        for(i=7;i<=NF;i++){
          l=l"\t"$i;
        }
        print l
     }' \
  "$OUTDIR/gene_counts.txt" > "$OUTDIR/gene_counts.matrix.tsv"

log "Wrote: $OUTDIR/gene_counts.matrix.tsv"
