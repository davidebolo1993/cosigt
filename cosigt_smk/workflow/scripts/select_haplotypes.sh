#!/bin/bash
# Read the haplotype names cosigt predicted out of a genotype table and emit
# them as a samtools faidx region list, checking each exists in the FASTA index.
#
#   select_haplotypes.sh <fasta.fai> <genotype.tsv> <columns> [pansn_prefix]
#
#   <columns>       "all" for every haplotype.N column in header order, or a
#                   single column name such as haplotype.1
#   [pansn_prefix]  when given, every FASTA entry whose name starts with it is
#                   emitted first, ahead of the haplotypes
#
# Replaces three near-identical awk programs that were embedded in cosigt.smk,
# where every brace had to be doubled for Snakemake's formatter.
#
# Exit codes: 2 no reference path matched the prefix, 3 no matching haplotype
# column in the table, 4 a predicted haplotype is absent from the FASTA index.
set -euo pipefail

fai="$1"
genotype="$2"
columns="$3"
pansn="${4:-}"

awk -F '\t' -v want="$columns" -v pansn="$pansn" '
FNR == NR {
  len[$1] = $2
  if (pansn != "" && index($1, pansn) == 1) refs[++n_refs] = $1
  next
}
FNR == 1 {
  for (i = 1; i <= NF; i++) {
    if (want == "all") {
      if ($i ~ /^haplotype\.[0-9]+$/) cols[++n_cols] = i
    } else if ($i == want) {
      cols[++n_cols] = i
    }
  }
  next
}
NF > 0 {
  # cosigt writes a single data row; keep the last one seen regardless
  for (i = 1; i <= NF; i++) last[i] = $i
}
END {
  if (pansn != "" && n_refs < 1) {
    print "No reference path starts with " pansn > "/dev/stderr"
    exit 2
  }
  if (n_cols < 1) {
    if (want == "all") {
      print "No haplotype.N columns found in " FILENAME > "/dev/stderr"
    } else {
      print "Missing " want " in genotype table " FILENAME > "/dev/stderr"
    }
    exit 3
  }
  for (i = 1; i <= n_refs; i++) print refs[i]
  for (i = 1; i <= n_cols; i++) {
    hap = last[cols[i]]
    if (!(hap in len)) {
      print "Haplotype " hap " not found in FASTA index" > "/dev/stderr"
      exit 4
    }
    print hap
  }
}
' "$fai" "$genotype"
