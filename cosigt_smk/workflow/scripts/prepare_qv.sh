#!/bin/bash

input_table=$1
input_fasta=$2
outdir=$3
mkdir -p $outdir
#index
samtools faidx $input_fasta
#read table line-by-line
#excluding header
tail -n +2 $input_table | while read line; do
    sample=$(echo "$line" | cut -d$'\t' -f 1)
    h1p=$(echo "$line" | cut -d$'\t' -f 2)
    h2p=$(echo "$line" | cut -d$'\t' -f 3)
    h1t=$(echo "$line" | cut -d$'\t' -f 4)
    h2t=$(echo "$line" | cut -d$'\t' -f 5)
    mkdir -p $outdir/$sample
    if [ $h1t = "missing" ] || [ $h2t = "missing" ]; then
        echo -e "missing.fasta\tmissing.fasta" >> "$outdir/$sample/ids.tsv"
    else
        fs=$(echo "$h1p" "$h2p" "$h1t" "$h2t" | sort | uniq)
        for f in $fs; do
            #sanitize
            sf=$(echo $f | sed 's/[^a-zA-Z0-9._-]/_/g')
            samtools faidx "$input_fasta" "$f" > "$outdir/$sample/$sf.fasta"
        done
        #sanitize
        sf1p=$(echo "$h1p" | sed 's/[^a-zA-Z0-9._-]/_/g')
        sf2p=$(echo "$h2p" | sed 's/[^a-zA-Z0-9._-]/_/g')
        sf1t=$(echo "$h1t" | sed 's/[^a-zA-Z0-9._-]/_/g')
        sf2t=$(echo "$h2t" | sed 's/[^a-zA-Z0-9._-]/_/g')

	    if [ -f "$outdir/$sample/ids.tsv" ]; then
		    rm "$outdir/$sample/ids.tsv"
	    fi

        echo -e "${sf1p}.fasta\t${sf1t}.fasta" >> "$outdir/$sample/ids.tsv"
        echo -e "${sf2p}.fasta\t${sf2t}.fasta" >> "$outdir/$sample/ids.tsv"
        echo -e "${sf1p}.fasta\t${sf2t}.fasta" >> "$outdir/$sample/ids.tsv"
        echo -e "${sf2p}.fasta\t${sf1t}.fasta" >> "$outdir/$sample/ids.tsv"
    fi
done
