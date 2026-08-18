#!/bin/bash
# Extract the sequences the benchmark needs from the orientation-normalised
# region FASTA: each sample's predicted and true haplotypes, and - for
# leave_all_out - every haplotype the genotyping graph offered, which is the
# candidate set the oracle picks from.
#
#   benchmark_prepare.sh <flipped.fasta> <panel.fai> <mode> <outdir> <genotype.tsv>...
#
# <mode> is leave_zero_out or leave_all_out.
#
# Writes:
#   <outdir>/manifest.tsv                     one row per sample
#   <outdir>/sequences/<sample>/{pred,true}{1,2}.fasta
#   <outdir>/sequences/panel/<n>.fasta        leave_all_out only
#   <outdir>/sequences/panel/index.tsv        <n> -> haplotype name
#
# The comparison is defined for two predicted against two true haplotypes, so
# regions whose ploidy is not 2, and samples whose own haplotypes are not in the
# truth graph exactly twice, are recorded with a status and skipped.
set -euo pipefail

fasta="$1"
panel_fai="$2"
mode="$3"
outdir="$4"
shift 4

seqdir="$outdir/sequences"
mkdir -p "$seqdir"
manifest="$outdir/manifest.tsv"
printf 'sample\tstatus\thap_1_pred\thap_2_pred\tcluster_1_pred\tcluster_2_pred\thap_1_true\thap_2_true\n' > "$manifest"

[ -f "$fasta.fai" ] || samtools faidx "$fasta"

# The oracle may only choose haplotypes that existed in the genotyping graph.
# Numbered files avoid having to sanitise PanSN names into filenames.
if [ "$mode" = "leave_all_out" ]; then
    paneldir="$seqdir/panel"
    mkdir -p "$paneldir"
    : > "$paneldir/index.tsv"
    n=0
    while IFS= read -r name; do
        n=$((n + 1))
        samtools faidx "$fasta" "$name" > "$paneldir/$n.fasta"
        printf '%s\t%s\n' "$n" "$name" >> "$paneldir/index.tsv"
    done < <(cut -f1 "$panel_fai")
fi

for tsv in "$@"; do
    # cosigt writes exactly one data row: sample id, haplotype.N, optional cluster.N
    read -r sample pred1 pred2 clust1 clust2 < <(
        awk -F '\t' '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == "haplotype.1") h1 = i
                    if ($i == "haplotype.2") h2 = i
                    if ($i == "cluster.1")   c1 = i
                    if ($i == "cluster.2")   c2 = i
                }
                next
            }
            NF > 0 {
                s = $1; p1 = (h1 ? $h1 : "NA"); p2 = (h2 ? $h2 : "NA")
                k1 = (c1 ? $c1 : "NA"); k2 = (c2 ? $c2 : "NA")
            }
            END { print s, p1, p2, k1, k2 }
        ' OFS='\t' "$tsv"
    )

    if [ "$pred1" = "NA" ] || [ "$pred2" = "NA" ]; then
        printf '%s\tnot_diploid\t%s\t%s\t%s\t%s\tNA\tNA\n' \
            "$sample" "$pred1" "$pred2" "$clust1" "$clust2" >> "$manifest"
        continue
    fi

    # True haplotypes: the sample's own paths. PanSN names are
    # sample#haplotype#contig, so anchor on the first field to avoid matching a
    # sample whose name is a prefix of another. Collected with a plain loop
    # rather than mapfile, which needs bash 4.
    true1=""; true2=""; n_true=0
    while IFS= read -r name; do
        n_true=$((n_true + 1))
        case "$n_true" in
            1) true1="$name" ;;
            2) true2="$name" ;;
        esac
    done < <(cut -f1 "$fasta.fai" | awk -F'#' -v s="$sample" '$1 == s' | sort)

    if [ "$n_true" -ne 2 ]; then
        printf '%s\tmissing\t%s\t%s\t%s\t%s\tNA\tNA\n' \
            "$sample" "$pred1" "$pred2" "$clust1" "$clust2" >> "$manifest"
        continue
    fi

    sdir="$seqdir/$sample"
    mkdir -p "$sdir"
    samtools faidx "$fasta" "$pred1" > "$sdir/pred1.fasta"
    samtools faidx "$fasta" "$pred2" > "$sdir/pred2.fasta"
    samtools faidx "$fasta" "$true1" > "$sdir/true1.fasta"
    samtools faidx "$fasta" "$true2" > "$sdir/true2.fasta"

    printf '%s\tok\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$pred1" "$pred2" "$clust1" "$clust2" "$true1" "$true2" >> "$manifest"
done
