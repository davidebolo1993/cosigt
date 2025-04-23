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
    stretcher -asequence "$indir/$fa1" -bsequence "$indir/$fa2" -outfile "$indir/out.txt"
    length=$(grep "^# Length:" "$indir/out.txt" | awk '{print $3}')
    identity=$(grep "^# Identity:" "$indir/out.txt" | sed -E 's/.*Identity:[[:space:]]*([0-9]+)\/[0-9]+.*/\1/')
    delta=$((length - identity))
    delta_max=$(echo - | awk -v d="$delta" 'BEGIN { printf "%.4f", (d < 0.5) ? 0.5 : d }')
    QV=$(awk -v delta="$delta_max" -v len="$length" 'BEGIN {
        ratio = delta / len;
        qv = -10 * log(ratio) / log(10);
        printf "%.2f", qv
    }')
    echo -e "$fa1\t$fa2\t$QV" >> "$indir/qv.tmp.tsv"
done
rm "$indir/out.txt"

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
