#!/bin/bash

bedpe_file="$1"
distance="$2"

awk -v d="$distance" '
function abs(x) { return (x < 0) ? -x : x }

function flush() {
  if (qname != "") {
    print qname, qstart, qend, tname, tstart, tend, ".", "0", ".", "."
  }
}

BEGIN { OFS = "\t" }

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
