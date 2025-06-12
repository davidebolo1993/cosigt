#!/bin/bash

# Usage: check_flanks_and_uniqueness.sh <bedpe_file> <region_string> <flank_size>

bedpe_file="$1"
region_bed="$2"
flank_size="$3"

# Extract start and end from the region string (safe even with underscores in chr)
target_start=$(cat $region_bed| cut -f 2)
target_end=$(cat $region_bed| cut -f 3)

awk -v ts="$target_start" -v te="$target_end" -v size="$flank_size" '
BEGIN {
  l_start = ts
  l_end   = ts + size
  r_start = te - size
  r_end   = te
}
{
  if ($5 <= l_start && $6 >= l_end &&
      $5 <= r_start && $6 >= r_end) {
    print
  }
}
' "$bedpe_file"

awk '{ print $1 }' "$bedpe_file" | sort -u | wc -l | awk '{ exit ($1 == 1 ? 0 : 1) }'