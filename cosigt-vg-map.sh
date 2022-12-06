#!/bin/bash

base=$1
in=$2
out=$3
t=$4

mkdir -p $(dirname $out)

#run giraffe
echo "vg map"
vg map -x $base.xg -g $base.gcsa -f $in -% -t $t > $out.gaf

#get node coverage
echo "gafpack"
gafpack -g $base.gfa -a $out.gaf | pigz > $out.gafpack.gz

#calculate genotype
echo "cosigt.py"
cosigt.py $base.paths.tsv.gz $out.gafpack.gz $out
