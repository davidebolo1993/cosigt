#!/bin/bash

#initialize config
echo -e "samples: config/samples.tsv"> config/config.yaml

#link cram files and populate samples.tsv
cram_base=$(readlink -f $1)

if [ $# -eq 4 ]; then

	cram=$(ls $cram_base/*.*am)

else

	cram=$(ls $cram_base/*.*am | grep -v -f $5)

fi

echo -e "sample_id\tcram" > config/samples.tsv
mkdir -p resources/cram

for c in $cram; do 

	filename=$(basename -- "$c")
	ln -sf $c resources/cram/$filename
	idx=$(ls $c"."*i)
	ext="${idx##*.}"
	ln -sf $c"."$ext resources/cram/$filename"."$ext
	filename=$(echo $filename | cut -d "." -f 1)
	echo -e "$filename\t$c" >> config/samples.tsv

done

#add reference
mkdir -p resources/ref
ref=$(readlink -f $2)
ref_dir=$(dirname $ref)
filename=$(basename -- "$ref")
ln -sf $ref resources/ref/$filename
echo -e "reference: resources/ref/$filename" >> config/config.yaml

#add graph
mkdir -p resources/graph
graph=$(readlink -f $3)
graph_dir=$(dirname $graph)
filename=$(basename -- "$graph")
ln -sf $graph resources/graph/$filename 
echo -e "graph: resources/graph/$filename" >> config/config.yaml

#add region
echo -e "region: \"$4\"" >> config/config.yaml

#export singularity bind paths
echo $cram_base","$ref_dir","$graph_dir > singularity_bind_paths.csv

