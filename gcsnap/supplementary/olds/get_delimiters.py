import os
import argparse
import gzip
import csv
import glob

# python3 get_delimiters.py --input-dir ../../../data/MGnify_2023_02 --folder-name contig_map --entity-name contig
# python3 get_delimiters.py --input-dir ../../../data/MGnify_2023_02 --folder-name seq_metadata --entity-name sequence

def parse_args():
    """
    Parses command-line arguments for processing MGnify files.
    """
    parser = argparse.ArgumentParser(
        description="Process gzipped TSV files to create a 'delimiters.csv' index."
    )
    
    # We use new flag names to be more "Pythonic" like the first script
    parser.add_argument(
        "-d", "--input-dir", 
        type=str, 
        required=True,
        help="Base MGnify data directory (e.g., ../../data/MGnify)"
    )
    parser.add_argument(
        "-f", "--folder-name", 
        type=str, 
        required=True,
        help="Subfolder name to process (e.g., contig_map)"
    )
    parser.add_argument(
        "-e", "--entity-name", 
        type=str, 
        required=True,
        help="Name of the ID entity (e.g., protein, contig)"
    )
    parser.add_argument(
        "-s", "--block-size", 
        type=int, 
        default=500_000,
        help="Number of lines per block (default: 500,000)"
    )

    return parser.parse_args()

def process_file_blocks(file_path, file_name, block_size):
    """
    Reads a single gzipped file and yields its block information.
    
    This function is a generator, making it memory-efficient.
    It reads the file line-by-line and yields a summary
    for each block.
    """
    first_id = None
    last_id = None
    block_number = 1
    
    # 'rt' = Read as Text (handles decompression automatically)
    with gzip.open(file_path, 'rt', encoding='utf-8') as f:
        for i, line in enumerate(f):
            # 1. Check if this is the start of a new block
            if i % block_size == 0:
                # If this isn't the very first line, we just finished a block.
                # Yield the info for the block that just ended.
                if i > 0:
                    yield (file_name, block_number, first_id, last_id)
                    block_number += 1
                
                # Start the new block
                first_id = line.split('\t')[0].strip()
            
            # 2. Always update the last_id for the current block
            last_id = line.split('\t')[0].strip()

    # 3. END block: After the loop, yield the final (potentially partial) block
    if first_id is not None:
        yield (file_name, block_number, first_id, last_id)

def main():
    
    args = parse_args()
    
    # 1. Define file paths based on arguments
    #    This is the Python equivalent of:
    #    ../../data/MGnify/"${folder_name}"/
    process_dir = os.path.join(args.input_dir, args.folder_name)
    
    #    This is the Python equivalent of:
    #    mgy_"${folder_name}"_*.tsv.gz
    search_pattern = os.path.join(process_dir, f"mgy_{args.folder_name}_*.tsv.gz")
    
    #    This is the Python equivalent of:
    #    ../../data/MGnify/${folder_name}/delimiters.csv
    output_file = os.path.join(process_dir, "delimiters.csv")
    
    # Use glob to find all files matching the pattern (like the * in Bash)
    file_paths = glob.glob(search_pattern)

    if not file_paths:
        print(f"No files found matching pattern: {search_pattern}")
        return

    print(f"Found {len(file_paths)} files to process.")

    # 2. Write header and process files
    #    'w' = write (overwrites existing file), newline='' is standard for csv
    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as f_out:
            # Use the csv module for robust CSV writing
            writer = csv.writer(f_out)
            
            # Write the header row (like the 'echo "..." > $output_file')
            header = ["filename", "block_number", f"first_{args.entity_name}_id", f"last_{args.entity_name}_id"]
            writer.writerow(header)
            
            # 3. Process each file (like the 'for file_path in ...' loop)
            for file_path in file_paths:
                file_name = os.path.basename(file_path)
                print(f"Processing {file_name} ...")
                
                # Get all blocks from our generator and write them to the file
                # (This replaces the zcat | awk >> $output_file)
                for block_data in process_file_blocks(file_path, file_name, args.block_size):
                    writer.writerow(block_data)

        print(f"\nOutput saved to {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()