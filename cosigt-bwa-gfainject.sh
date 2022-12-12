#!/bin/bash

base=$1
in=$2
out=$3
t=$4

mkdir -p $(dirname $out)

#run bwa
echo "bwa mem"
bwa mem -t $t $base.fa $in >$out.sam
samtools view -b $out.sam >$out.bam

#gfainject
echo "gfainject"
gfainject --gfa $base.gfa --bam $out.bam >$out.gaf

#get node coverage
echo "gafpack"
gafpack -g $base.gfa -a $out.gaf | pigz > $out.gafpack.gz

#calculate genotype
echo "cosigt.py"
cosigt.py $base.paths.tsv.gz $out.gafpack.gz $out
