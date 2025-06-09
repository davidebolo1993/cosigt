#!/bin/bash

bedpe_file="$1"
distance="$2"

awk -v d="$distance" '
function flush() {
  if (qname != "") {
    print qname, qstart, qend, tname, tstart, tend, ".", "0", "+", "+" 
  }
}
BEGIN { OFS = "\t" }
{
  if ($1 == qname && $4 == tname &&
      $2 - qend <= d && $5 - tend <= d) {
    # Extend current block
    qend = ($3 > qend ? $3 : qend)
    tend = ($6 > tend ? $6 : tend)
  } else {
    # Flush previous block
    flush()
    # Start new block
    qname = $1
    qstart = $2
    qend = $3
    tname = $4
    tstart = $5
    tend = $6
  }
}
END { flush() }
' "$bedpe_file"