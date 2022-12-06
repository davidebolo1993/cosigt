#!/bin/bash

g=$1
base=$2
t=$3
scratch=$4

#chop
echo "odgi chop"
odgi chop -i $g -c 32 -o - -t $t | odgi view -i - -g > $base.gfa

#haplotype binary matrix
echo "odgi build"
odgi build -g $base.gfa -t $t -o - | odgi paths -i - -H | cut -f 1,4- | pigz > $base.paths.tsv.gz

#index for giraffe
echo "vg autoindex"
vg autoindex -w giraffe -g $base.gfa -t $t -p $base -T $scratch

