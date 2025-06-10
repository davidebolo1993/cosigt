#!/bin/bash
# this is mainly a contribution by chiara.paleni@fht.org

INPUT_ASM="$1"
INPUT_PROTEINS="$2"
REF_PATH="$3"
OUTPUT_DIR="$4"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# split multi-fasta into single-fasta files
awk -v out=$OUTPUT_DIR '/^>/ {
    if (NR > 1) close(filename)
    filename = substr($1, 2) ".fasta"
    print > filename
    next
}
{
    print > filename
}' "$INPUT_ASM"

# run miniprot for each single assembly
for fasta_file in *.fasta; do
    base_name="${fasta_file%.fasta}"
    miniprot -u "$fasta_file" "$INPUT_PROTEINS" > "${base_name}.paf"
done

# re-orient PAF files
# actual function for PAF-flipping
flip_paf() {
    local input_file="$1"
    awk 'BEGIN{FS=OFS="\t"} {
        # Store original values
        temp_start = $8
        temp_end = $9
        temp_strand = $5
        
        # Flip coordinates: new_start = t.len - old_end, new_end = t.len - old_start
        $8 = $7 - temp_end
        $9 = $7 - temp_start
        
        # Flip strand
        if (temp_strand == "+") $5 = "-"
        else if (temp_strand == "-") $5 = "+"
        
        print
    }' "$input_file"
}

# get orientation profile - unchanged, as it correctly extracts gene, strand, count
get_orientation_profile() {
    local file="$1"
    awk '$5 != "*" {print $1, $5, $3}' "$file" | sort | uniq -c | awk '{print $2, $3, $1}'
}

REF_PAF=$(ls *paf | grep "$REF_PATH")

# get reference orientation profile
get_orientation_profile "$REF_PAF" > reference.profile

# process all the other paf files
for paf_file in *.paf; do
    base_name="${paf_file%.paf}"
    get_orientation_profile "$paf_file" > current.profile

    comparison_result=$(awk '
        # Process reference.profile first
        NR==FNR {
            # Store count for each gene_id_strand combination
            ref_counts[$1"_"$2] = $3;
            # Keep track of all unique gene_ids encountered
            all_gene_ids[$1] = 1;
            next;
        }
        # Process current.profile
        {
            # Store count for each gene_id_strand combination
            current_counts[$1"_"$2] = $3;
            # Keep track of all unique gene_ids encountered
            all_gene_ids[$1] = 1;
        }
        END {
            matches = 0;
            mismatches = 0;

            # Iterate through all unique gene IDs found in either profile
            for (gene_id in all_gene_ids) {
                # Check for '+' strand
                ref_plus_key = gene_id"_+"
                current_plus_key = gene_id"_+"
                
                ref_plus_count = (ref_plus_key in ref_counts) ? ref_counts[ref_plus_key] : 0;
                current_plus_count = (current_plus_key in current_counts) ? current_counts[current_plus_key] : 0;
                
                if (ref_plus_count == current_plus_count) {
                    matches++;
                } else {
                    mismatches++;
                }

                # Check for '-' strand
                ref_minus_key = gene_id"_-";
                current_minus_key = gene_id"_-";
                
                ref_minus_count = (ref_minus_key in ref_counts) ? ref_counts[ref_minus_key] : 0;
                current_minus_count = (current_minus_key in current_counts) ? current_counts[current_minus_key] : 0;

                if (ref_minus_count == current_minus_count) {
                    matches++;
                } else {
                    mismatches++;
                }
            }

            # Decision logic based on R script: if mismatches are more than matches, flip
            if (mismatches > matches) {
                print "flip";
            } else {
                print "keep";
            }
        }
    ' reference.profile current.profile)
    
    # Create oriented PAF file
    if [[ "$comparison_result" == "flip" ]]; then
        flip_paf "$paf_file" > "${base_name}_oriented.paf"
    else
        cp "$paf_file" "${base_name}_oriented.paf"
    fi

    rm "$paf_file"
    rm current.profile
done

# clean
rm reference.profile