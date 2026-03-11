#!/bin/bash

# Input arguments with flags
while getopts "f:e:" opt; do
    case $opt in
        f) folder_name="$OPTARG" ;;
        e) entity_name="$OPTARG" ;;
        *) 
            echo "Usage: $0 -f <folder_name> -e <entity_name>"
            exit 1
            ;;
    esac
done

# Check if both arguments are provided
if [ -z "$folder_name" ] || [ -z "$entity_name" ]; then
    echo "Usage: $0 -f <folder_name> -e <entity_name>"
    exit 1
fi

# File paths
file_paths=(../../data/MGnify/"${folder_name}"/mgy_"${folder_name}"_*.tsv.gz)
output_file="../../data/MGnify/${folder_name}/delimiters.csv"

# Write header to output file
echo "filename,block_number,first_${entity_name}_id,last_${entity_name}_id" > "$output_file"

# Process each file
for file_path in "${file_paths[@]}"; do
    file_name=$(basename "$file_path")
    echo "Processing $file_name ..."
    
    # Use zcat and awk to process the file in 0.5M (500,000) line blocks.
    zcat "$file_path" | awk -v block_size=500000 -v fname="$file_name" '{
        # If starting a new block, store the first key
        if ( (NR - 1) % block_size == 0 ) {
            first = $1;
            block = int((NR-1) / block_size) + 1;
        }
        last = $1;  # Always update the last key of the current block
        # When we hit the end of a block, print out the block info
        if ( NR % block_size == 0 ) {
            print fname "," block "," first "," last;
        }
    }
    END {
        # If the file ended in the middle of a block, print that block info as well.
        if (NR % block_size != 0) {
            block = int(NR / block_size) + 1;
            print fname "," block "," first "," last;
        }
    }' >> "$output_file"
done

echo "Output saved to $output_file"