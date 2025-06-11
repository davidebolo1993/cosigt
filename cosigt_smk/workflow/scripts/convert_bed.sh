#!/bin/bash

bed_file="$1"

awk -F'\t' '{
  split($4, g, ":");
  strand = ($6 == "+" ? 1 : 0);
  print $1 "\t" g[1] "\t" $2 "\t" $3 "\t" strand
}' "$bed_file" | 
awk 'BEGIN {OFS="\t"; print "molecule","gene","start","end","strand"} 1'