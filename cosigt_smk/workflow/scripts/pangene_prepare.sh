#!/bin/bash
# A contribution by chiara.paleni@fht.org

#   pangene_prepare.sh <alleles.fasta.gz> <proteins.faa.gz> <pansn_prefix> <outdir> [threads]
set -euo pipefail

INPUT_ASM="$1"
INPUT_PROTEINS="$2"
REF_PATH="$3"
OUTPUT_DIR="$4"
THREADS="${5:-1}"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# split gzipped multi-fasta into single bgzip-compressed fasta files
zcat "$INPUT_ASM" | awk -v out=$OUTPUT_DIR '/^>/ {
    if (NR > 1) close(cmd)
    filename = out"/"substr($1, 2) ".fasta"
    cmd = "bgzip > " filename ".gz"
    print | cmd
    next
}
{
    print | cmd
}'

# run miniprot for each single assembly, in parallel: the invocations are
# independent and each writes its own PAF
ls ./*.fasta.gz \
  | sed 's|^\./||; s|\.fasta\.gz$||' \
  | xargs -P "$THREADS" -I{} sh -c \
      'miniprot -u "{}.fasta.gz" "'"$INPUT_PROTEINS"'" | bgzip > "{}.paf.gz"'

# function to flip PAF (works on uncompressed input, so uncompress first)
flip_paf() {
    local input_file="$1"
    zcat "$input_file" | awk 'BEGIN{FS=OFS="\t"} {
        temp_start = $8
        temp_end = $9
        temp_strand = $5

        $8 = $7 - temp_end
        $9 = $7 - temp_start

        if (temp_strand == "+") $5 = "-"
        else if (temp_strand == "-") $5 = "+"

        print
    }'
}

# get orientation profile
get_orientation_profile() {
    local file="$1"
    zcat "$file" | awk '$5 != "*" {print $1, $5, $3}' | sort | uniq -c | awk '{print $2, $3, $1}'
}

# Match the PanSN prefix at the start of the contig name. `ls | grep` matched
# anywhere in the path and returned every hit, so more than one match produced a
# multi-line REF_PAF that broke the profile call below.
REF_PAF=""
for candidate in *.paf.gz; do
    case "$(basename "$candidate" .paf.gz)" in
        "$REF_PATH"*) REF_PAF="$candidate"; break ;;
    esac
done
if [ -z "$REF_PAF" ]; then
    echo "No contig starting with the reference prefix '$REF_PATH' in $INPUT_ASM" >&2
    exit 1
fi

# reference profile
get_orientation_profile "$REF_PAF" > reference.profile

# process other PAFs 
for paf_file in *.paf.gz; do
    base_name="${paf_file%.paf.gz}"
    get_orientation_profile "$paf_file" > current.profile

    comparison_result=$(awk '
        NR==FNR {
            ref_counts[$1"_"$2] = $3;
            all_gene_ids[$1] = 1;
            next;
        }
        {
            current_counts[$1"_"$2] = $3;
            all_gene_ids[$1] = 1;
        }
        END {
            matches = 0; mismatches = 0;
            for (gene_id in all_gene_ids) {
                ref_plus = ref_counts[gene_id"_+"] + 0;
                cur_plus = current_counts[gene_id"_+"] + 0;
                ref_minus = ref_counts[gene_id"_-"] + 0;
                cur_minus = current_counts[gene_id"_-"] + 0;
                matches += (ref_plus == cur_plus) + (ref_minus == cur_minus);
                mismatches += (ref_plus != cur_plus) + (ref_minus != cur_minus);
            }
            print (mismatches > matches ? "flip" : "keep");
        }
    ' reference.profile current.profile)

    # Both branches emit bgzip. Previously flip wrote plain gzip while keep
    # cp'd the bgzip input, so one rule produced two container formats
    # depending on the orientation call.
    if [[ "$comparison_result" == "flip" ]]; then
        flip_paf "$paf_file" | bgzip > "${base_name}_oriented.paf.gz"
    else
        zcat "$paf_file" | bgzip > "${base_name}_oriented.paf.gz"
    fi

    rm "$paf_file" current.profile
done

rm reference.profile