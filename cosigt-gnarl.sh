#!/bin/bash

base=$1
in=$2
out=$3

mkdir -p $(dirname $out)

#run giraffe
echo "vg giraffe"
vg giraffe -Z $base.giraffe.gbz -m  $base.min -d $base.dist -f $in -o gaf > $out.gaf

#get node coverage
echo "gafpack"
gafpack -g $base.gfa -a $out.gaf | pigz > $out.gafpack.gz

#calculate genotype
echo "cosigt.py"
./cosigt-gnarl.py $base.paths.tsv.gz $out.gafpack.gz $out
