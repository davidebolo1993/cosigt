#!/bin/bash

# Usage: check_flanks.sh <bedpe_file> <region_string> <flank_pct>

bedpe_file="$1"
region="$2"
flank_pct="$3"

# Extract start and end from the region string (safe even with underscores in chr)
rev_coords=$(echo "$region" | rev | cut -d"_" -f1,2 | rev)
target_start=$(echo "$rev_coords" | cut -d"_" -f1)
target_end=$(echo "$rev_coords" | cut -d"_" -f2)

awk -v ts="$target_start" -v te="$target_end" -v pct="$flank_pct" '
BEGIN {
  len = te - ts
  flank = len * pct
  l_start = ts
  l_end   = ts + flank
  r_start = te - flank
  r_end   = te
}
{
  if ($5 <= l_start && $6 >= l_end &&
      $5 <= r_start && $6 >= r_end) {
    print
  }
}
' "$bedpe_file"
