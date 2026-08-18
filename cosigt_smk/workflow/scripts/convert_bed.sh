#!/bin/bash
# Reshape `pangene --bed` output into the gggenes input plotgggenes.r expects.
#
#   convert_bed.sh <bed|->
#
# Column 4 of the pangene BED is gene:transcript, of which only the gene part is
# kept; strand becomes 1/0. A single awk pass, header included, rather than
# cat into awk into awk.
set -euo pipefail

bed_file="$1"

awk -F '\t' '
BEGIN { OFS = "\t"; print "molecule", "gene", "start", "end", "strand" }
{
  split($4, g, ":")
  print $1, g[1], $2, $3, ($6 == "+" ? 1 : 0)
}
' "$bed_file"
