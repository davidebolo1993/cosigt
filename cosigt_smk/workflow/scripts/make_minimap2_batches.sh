#!/bin/bash
# Split an assembly .fai into one contig-name list per PanSN sample, which is
# the unit minimap2 alignment is then parallelised over.
#
#   make_minimap2_batches.sh <assembly.fai> <output_dir>
#
# A single awk pass writing straight to the per-sample files, rather than a grep
# over the whole index plus an append per contig.
set -euo pipefail

input_file_index="$1"
output_dir="$2"

mkdir -p "$output_dir"

# PanSN names are sample#haplotype#contig, so the sample is the first field.
awk -F '\t' -v out="$output_dir" '
{
  split($1, parts, "#")
  print $1 > (out "/" parts[1] ".txt")
}
' "$input_file_index"
