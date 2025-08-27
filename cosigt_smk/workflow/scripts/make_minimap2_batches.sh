#!/bin/bash

input_file_index=$1
output_dir=$2 
mkdir -p "$output_dir"
cut -f 1 $input_file_index | cut -d "#" -f 1 | sort | uniq | while read f; do
    input_file=$(echo "${input_file_index%.*}")
    grep -w $f $input_file_index | cut -f 1 | while read g; do
        echo $g >> $output_dir/$f.txt
    done
done