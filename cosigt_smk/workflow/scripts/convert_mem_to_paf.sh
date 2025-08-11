#!/bin/bash

target_info_file="$1"

awk '
# This block runs first, processing the target_info_file
FNR==NR {
    target_info[$1] = $2;
    next;
}

# This block runs for every line from stdin
{
    query_name = $1;
    query_start = $2;
    query_end = $3;
    
    mem_length = query_end - query_start;
    
    for (i = 6; i <= NF; i++) {
        split($i, parts, ":");
        
        num_parts = length(parts);
        strand = parts[num_parts - 1];
        target_start = parts[num_parts];
        
        target_name = parts[1];
        for (j = 2; j < num_parts - 1; j++) {
            target_name = target_name ":" parts[j];
        }

        target_end = target_start + mem_length;
        query_length = 150;
        target_length = target_info[target_name];
        if (target_length == "") {
            target_length = 0;
        }
        
        print query_name,
              query_length,
              query_start,
              query_end,
              strand,
              target_name,
              target_length,
              target_start,
              target_end,
              mem_length,
              mem_length,
              255,
              "cg:Z:" mem_length "=";
    }
}
' OFS='\t' "$target_info_file" -
