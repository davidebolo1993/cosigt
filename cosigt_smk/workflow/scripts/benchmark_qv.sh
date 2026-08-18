#!/bin/bash
# Score cosigt's predicted haplotypes against each sample's true haplotypes with
# edlib, and emit one row per sample for this region.
#
#   benchmark_qv.sh <prepared_dir> <region> <gene_name> <mode> <threads> <output.tsv>
#
# QV_sum_pred is the achieved score: the two predicted haplotypes against the two
# true ones, taking whichever of the two assignments scores higher, so haplotype
# order does not matter.
#
# In leave_all_out the true haplotypes are absent from the genotyping graph, so a
# perfect score is impossible by construction and the achieved value alone cannot
# separate a bad prediction from a panel that simply held nothing close. Each true
# haplotype is therefore also compared against every haplotype the genotyping
# graph did offer, and the best match for each is taken:
#
#   QV_sum_best  = best(true_1 vs panel) + best(true_2 vs panel)
#   QVfrac      = QV_sum_pred / QV_sum_best
#
# QVfrac near 1 means cosigt picked about the best pair available to it, so a low
# absolute QV is the panel's limitation rather than the genotyper's.
#
# compute_qv reports QV = -10*log10(edit_distance/alignment_length) with the edit
# distance floored at 0.5, so an exact match yields a length-dependent ceiling
# rather than infinity. error_rate is recovered as 10^(-QV/10), inverting that
# definition exactly.
set -euo pipefail

indir="$1"
region="$2"
gene="$3"
mode="$4"
threads="$5"
output="$6"

seqdir="$indir/sequences"
paneldir="$seqdir/panel"

printf 'sample\tregion\tgene_name\thap_1_pred\thap_2_pred\tcluster_1_pred\tcluster_2_pred\thap_1_true\thap_2_true\tQV_1_pred\tQV_2_pred\tQV_sum_pred\terror_rate_1_pred\terror_rate_2_pred\thap_1_best\thap_2_best\tQV_1_best\tQV_2_best\tQV_sum_best\tQVfrac\n' > "$output"

# Best-scoring panel haplotype for one true haplotype, as "<qv>\t<name>".
# Every candidate is scored; the panel index maps file number to PanSN name.
best_against() {
    local truth="$1"
    cut -f1 "$paneldir/index.tsv" \
        | xargs -P "$threads" -I{} sh -c \
            'printf "%s\t%s\n" "$(compute_qv "'"$paneldir"'/{}.fasta" "'"$truth"'")" "{}"' \
        | sort -g -r -k1,1 \
        | head -n 1 \
        | awk -F'\t' -v idx="$paneldir/index.tsv" '
            { qv = $1; n = $2 }
            END {
                while ((getline line < idx) > 0) {
                    split(line, f, "\t")
                    if (f[1] == n) { name = f[2]; break }
                }
                printf "%s\t%s\n", qv, name
            }'
}

tail -n +2 "$indir/manifest.tsv" | while IFS=$'\t' read -r sample status pred1 pred2 clust1 clust2 true1 true2; do
    if [ "$status" != "ok" ]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' \
            "$sample" "$region" "$gene" "$pred1" "$pred2" "$clust1" "$clust2" >> "$output"
        continue
    fi

    sdir="$seqdir/$sample"
    # achieved: both assignments of the two predictions to the two truths
    a1=$(compute_qv "$sdir/pred1.fasta" "$sdir/true1.fasta")
    a2=$(compute_qv "$sdir/pred2.fasta" "$sdir/true2.fasta")
    b1=$(compute_qv "$sdir/pred1.fasta" "$sdir/true2.fasta")
    b2=$(compute_qv "$sdir/pred2.fasta" "$sdir/true1.fasta")

    if [ "$mode" = "leave_all_out" ]; then
        IFS=$'\t' read -r m1 name1 < <(best_against "$sdir/true1.fasta")
        IFS=$'\t' read -r m2 name2 < <(best_against "$sdir/true2.fasta")
    else
        m1="NA"; m2="NA"; name1="NA"; name2="NA"
    fi

    awk -v s="$sample" -v r="$region" -v g="$gene" \
        -v p1="$pred1" -v p2="$pred2" -v c1="$clust1" -v c2="$clust2" \
        -v t1="$true1" -v t2="$true2" \
        -v a1="$a1" -v a2="$a2" -v b1="$b1" -v b2="$b2" \
        -v m1="$m1" -v m2="$m2" -v n1="$name1" -v n2="$name2" '
        BEGIN {
            OFS = "\t"
            if (a1 + a2 >= b1 + b2) { q1 = a1; q2 = a2; u1 = t1; u2 = t2 }
            else                    { q1 = b1; q2 = b2; u1 = t2; u2 = t1 }
            e1 = 10 ^ (-q1 / 10)
            e2 = 10 ^ (-q2 / 10)
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.8f\t%.8f",
                   s, r, g, p1, p2, c1, c2, u1, u2, q1, q2, q1 + q2, e1, e2
            if (m1 == "NA") {
                printf "\tNA\tNA\tNA\tNA\tNA\tNA\n"
            } else {
                mx = m1 + m2
                printf "\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n",
                       n1, n2, m1, m2, mx, (mx > 0 ? (q1 + q2) / mx : 0)
            }
        }' >> "$output"
done
