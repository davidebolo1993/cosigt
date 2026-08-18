#!/bin/bash
# Merge projected bedpe blocks and keep only the contigs that span the region.
#
#   bedpe_merge_filter.sh <bedpe|-> <merge_distance> <region> <flank_size>
#
# Replaces the former bedpe_merge.sh + bedpe_filter.sh pair, which ran as two
# passes over two gzip round-trips. Behaviour is unchanged:
#
#   merge   consecutive rows sharing query and target names are combined when
#           either gap is within <merge_distance>, taking the min start and max
#           end on both sides. Input must already be sorted, as bedtools sort
#           leaves it.
#   filter  a merged block is kept only if its target interval covers both
#           flanks of the region, i.e. spans start..start+flank and
#           end-flank..end.
#
# <region> is chrom_start_end; the coordinates are the last two underscore
# separated fields, so chromosome names may themselves contain underscores.
set -euo pipefail

bedpe_file="$1"
distance="$2"
region="$3"
flank_size="$4"

awk -v d="$distance" -v region="$region" -v size="$flank_size" '
function abs(x) { return (x < 0) ? -x : x }

# Emit the block held in the accumulator, if it spans both flanks.
function flush() {
  if (qname != "" &&
      tstart <= l_start && tend >= l_end &&
      tstart <= r_start && tend >= r_end) {
    print qname, qstart, qend, tname, tstart, tend, ".", "0", ".", "."
  }
}

BEGIN {
  OFS = "\t"
  n = split(region, parts, "_")
  ts = parts[n - 1] + 0
  te = parts[n] + 0
  l_start = ts
  l_end   = ts + size
  r_start = te - size
  r_end   = te
}

{
  if ($1 == qname && $4 == tname &&
      (abs($2 - qend) <= d || abs(qstart - $3) <= d)) {

    # Extend query block
    qstart = (qstart < $2) ? qstart : $2
    qend   = (qend   > $3) ? qend   : $3

    # Extend target block
    tstart = (tstart < $5) ? tstart : $5
    tend   = (tend   > $6) ? tend   : $6

  } else {
    flush()
    # New block
    qname = $1; qstart = $2; qend = $3
    tname = $4; tstart = $5; tend = $6
  }
}

END { flush() }
' "$bedpe_file"
