#!/bin/bash

indir=$1
sample=$(basename $indir)
outfile=$2

if [ -f "$indir/qv.tmp.tsv" ]; then
	rm "$indir/qv.tmp.tsv"
fi

cat "$indir/ids.tsv" | while read line; do
    fa1=$(echo "$line" | cut -d$'\t' -f 1)
    fa2=$(echo "$line" | cut -d$'\t' -f 2)
    if [ $fa1 = "missing.fasta" ]; then
        echo -e "missing\tmissing\t0" >> "$indir/qv.tmp.tsv"
    else
        QV=$(compute_qv "$indir/$fa1" "$indir/$fa2")
        echo -e "$fa1\t$fa2\t$QV" >> "$indir/qv.tmp.tsv"
    fi
done

#first combination
sum1=$(head -n 2 "$indir/qv.tmp.tsv" | awk '{sum += $3} END {printf "%.2f", sum}')
#second combination
sum2=$(tail -n 2 "$indir/qv.tmp.tsv" | awk '{sum += $3} END {printf "%.2f", sum}')

#which is best?
if awk -v s1="$sum1" -v s2="$sum2" 'BEGIN {exit !(s1 > s2)}'; then
    echo -e "$sample\t$(head -n 1 $indir/qv.tmp.tsv | cut -f 3)\t$(head -n 2 $indir/qv.tmp.tsv | tail -1 | cut -f 3)" > $outfile
else
    echo -e "$sample\t$(head -n 3 $indir/qv.tmp.tsv| tail -1 | cut -f 3)\t$(tail -1 $indir/qv.tmp.tsv | cut -f 3)" > $outfile
fi
