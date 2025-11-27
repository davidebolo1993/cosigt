#!/bin/bash

indir=$1
samples=$(ls -d $indir/* | rev | cut -d "/" -f 1 | rev)

for sample in $samples; do

    outdir="$indir/$sample"
    outfile="$indir/$sample/qv.tsv"
    mkdir -p $outdir

    if [ -f "$outdir/qv.tmp.tsv" ]; then
	    rm "$outdir/qv.tmp.tsv"
    fi

    cat "$outdir/ids.tsv" | while read line; do
        fa1=$(echo "$line" | cut -d$'\t' -f 1)
        fa2=$(echo "$line" | cut -d$'\t' -f 2)
        if [ $fa1 = "missing.fasta" ]; then
            echo -e "missing\tmissing\t0" >> "$outdir/qv.tmp.tsv"
        else
            QV=$(compute_qv "$outdir/$fa1" "$outdir/$fa2")
            echo -e "$fa1\t$fa2\t$QV" >> "$outdir/qv.tmp.tsv"
        fi
    done

    #first combination
    sum1=$(head -n 2 "$outdir/qv.tmp.tsv" | awk '{sum += $3} END {printf "%.2f", sum}')
    #second combination
    sum2=$(tail -n 2 "$outdir/qv.tmp.tsv" | awk '{sum += $3} END {printf "%.2f", sum}')

    #which is best?
    if awk -v s1="$sum1" -v s2="$sum2" 'BEGIN {exit !(s1 > s2)}'; then
        echo -e "$sample\t$(head -n 1 $outdir/qv.tmp.tsv | cut -f 3)\t$(head -n 2 $outdir/qv.tmp.tsv | tail -1 | cut -f 3)" > $outfile
    else
        echo -e "$sample\t$(head -n 3 $outdir/qv.tmp.tsv| tail -1 | cut -f 3)\t$(tail -1 $outdir/qv.tmp.tsv | cut -f 3)" > $outfile
    fi

done