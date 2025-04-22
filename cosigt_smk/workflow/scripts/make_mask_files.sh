#!/bin/bash

input_file=$1 
output_dir=$2 
target_length=$3

mkdir -p "$output_dir"

declare -a node_ids
declare -A node_lengths

while read -r node length; do
    node_ids+=("$node")
    node_lengths["$node"]=$length
done < "$input_file"

total_nodes=${#node_ids[@]}
index=0
file_count=0

while (( index < total_nodes )); do
    sum=0
    mask=()
    
    #until target length
    while (( index < total_nodes )); do
        node="${node_ids[index]}"
        length=${node_lengths[$node]}
        
        if (( sum + length > target_length )); then
            break
        fi
        
        sum=$((sum + length))
        mask+=("$node")
        ((index++))
    done

    output_file="$output_dir/mask_$(printf "%04d" $file_count).tsv"
    awk -v nodes="${mask[*]}" 'BEGIN {
        split(nodes, arr, " ");
    }
    {
        for (i in arr) {
            if ($1 == arr[i]) {
                print 1;
                next;
            }
        }
        print 0;
    }' "$input_file" > "$output_file"

    ((file_count++))
done