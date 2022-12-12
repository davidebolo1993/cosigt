#!/bin/bash

usage() { echo "Usage: $0 [-i <make index>] [-g <graph.gfa>] [-b <alignment.bam/.cram>] [-c <coordinates.string>] [-t <threads.int>] [-l <label.string>] [-r <reference.fasta>]" 1>&2; exit 1; }

while getopts ":g:b:c:t:l:r:i" opt; do
    case "${opt}" in
        g)
            g=${OPTARG}
            ;;
        b)
            b=${OPTARG}
            ;;
	    c)
            c=${OPTARG}
            ;;
	    t)
	        t=${OPTARG}
	        ;;
	    l)
	        l=${OPTARG}
	        ;;
	    r)
	        r=${OPTARG}
	        ;;
	    i)
	        make_index=true
	        ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${g}" ]; then
	echo "missing gfa file" && usage
elif [ -z "${b}" ]; then
	echo "missing bam file" && usage
#elif [ -z "${c}" ]; then
#	echo "missing coordinates" && usage
elif [ -z "${t}" ]; then
	t="1"
elif [ -z "${l}" ]; then
	l="output"
fi

#check if this is .cram or .bam
filename=$(basename -- "$b")
extension="${filename##*.}"

echo extension is $extension

echo "gfa file: " $g
echo "bam file: " $b
echo "coordinates - samtools format: " $c
echo "computing threads: " $t
echo "label: " $l
echo "reference: " $r
mkdir -p $l

if [ $extension == "cram" ]; then

	if [ -z "${r}" ]; then

        	echo "reference file is mandatory if providing .cram" && usage
	fi

	samtools view -O bam -o $l/region.bam -T $r -@ $t $b $c
	samtools index -@ $t $l/region.bam
    echo "samtools fastq"
    samtools fastq -@ $t $l/region.bam | pigz > $l/region.fastq.gz
elif [ $extension == "gz" ]; then
    # nothing to do
    cp $b $l/region.fa.gz
else
	samtools view -O bam -o $l/region.bam -@ $t $b $c
	samtools index -@ $t $l/region.bam
    echo "samtools fastq"
    samtools fastq -@ $t $l/region.bam | pigz > $l/region.fastq.gz
fi

if [ $make_index ]; then

#chop
echo "odgi chop"
odgi chop -i $g -c 32 -o - | odgi view -i - -g > $l/z.gfa

#haplotype binary matrix
echo "odgi build"
odgi build -g $l/z.gfa -o - | odgi paths -i - -H | cut -f 1,4- | pigz > $l/z.paths.tsv.gz

#index for giraffe
echo "vg autoindex"
vg autoindex -w giraffe -g $l/z.gfa -t $t -p $l/index -T $l

#run giraffe
echo "vg giraffe"
vg giraffe -Z $l/index.giraffe.gbz -m  $l/index.min -d $l/index.dist -f $l/region.fa.gz -o gaf > $l/x.gaf

#get node coverage
echo "gafpack"
gafpack -g $l/z.gfa -a $l/x.gaf | pigz > $l/x.gafpack.gz

#calculate genotype
echo "genotype.py"
cosigt $l/z.paths.tsv.gz $l/x.gafpack.gz $l
