#!/usr/bin/env python3
import os
import argparse
import glob
import pandas as pd
import gzip
import csv

MAKESNAP_DESCRIPTION = "Get the metadta files delimiters and convert them to Parquet format."
DELIMITER_SIZE = 500_000

def convert_tsv_gz_to_parquet(input_file, output_file, entity):
    """
    Reads a gzipped TSV file and writes it as a Parquet file.
    Assumes the TSV file has no header and two columns:
    - entity_id (first column)
    - data (second column)
    """

    if os.path.exists(output_file):
        print(f"Parquet file already exists for {input_file}. Skipping...")
        return

    try:
        print(f"Processing {input_file}...")
        # Read TSV file with gzip compression. Assuming no header.
        df = pd.read_csv(input_file, sep='\t', header=None, compression='gzip', 
                         names=[f'{entity}_id', 'data'], low_memory=False)
        # Write DataFrame to Parquet.
        df.to_parquet(output_file, index=False)
        print(f"Saved Parquet to {output_file}")
    except Exception as e:
        print(f"Failed to process {input_file}: {e}")

def convert_folder(folder, entity):
    """
    Converts all .tsv.gz files in the given folder to Parquet.
    Output Parquet files will be stored in the same folder.
    """
    # Get list of .tsv.gz files in the folder.
    tsv_files = glob.glob(os.path.join(folder, "*.tsv.gz"))
    if not tsv_files:
        print("No .tsv.gz files found in the specified folder.")
        return

    for tsv_file in tsv_files:
        base = os.path.splitext(os.path.basename(tsv_file))[0]
        # Remove the extra extension from '.tsv.gz' to just get the base name.
        if base.endswith('.tsv'):
            base = base[:-4]
        output_file = os.path.join(folder, f"{base}.parquet")

        if os.path.exists(output_file):
            print(f"Parquet file already exists for {tsv_file}. Skipping...")
            continue
        else:
            convert_tsv_gz_to_parquet(tsv_file, output_file, entity)

def process_file_blocks(file_path, file_name, block_size=DELIMITER_SIZE):
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

def get_delimiters(database, folder, entity):
    """
    Get the delimiters from the database.
    """

    target_dir = os.path.join(database, folder)
    output_file = os.path.join(target_dir, "delimiters.csv")

    if os.path.exists(output_file):
        print(f"Delimiters already exist for {folder}. Skipping...")
        return

    search_pattern = os.path.join(target_dir, f"mgy_{folder}_*.tsv.gz")

    file_paths = glob.glob(search_pattern)

    if not file_paths:
        print(f"No files found matching pattern: {search_pattern}")
        return

    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as f_out:
            # Use the csv module for robust CSV writing
            writer = csv.writer(f_out)
            
            # Write the header row (like the 'echo "..." > $output_file')
            header = ["filename", "block_number", f"first_{entity}_id", f"last_{entity}_id"]
            writer.writerow(header)
            
            # 3. Process each file (like the 'for file_path in ...' loop)
            for file_path in file_paths:
                file_name = os.path.basename(file_path)
                print(f"Processing {file_name} ...")
                
                # Get all blocks from our generator and write them to the file
                # (This replaces the zcat | awk >> $output_file)
                for block_data in process_file_blocks(file_path, file_name, DELIMITER_SIZE):
                    writer.writerow(block_data)

        print(f"\nOutput saved to {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")

def main():

    parser = argparse.ArgumentParser( description=MAKESNAP_DESCRIPTION )
    parser.add_argument("--database", default='../../../data/MGnify_2023_02',help="Folder containing .tsv.gz files to convert.")
    args = parser.parse_args()

    if not os.path.isdir(args.database):
        print(f"The folder {args.database} does not exist or is not a directory.")
        return

    ###################################################################################
    ### PARQUET CONVERSION
    ######################

    targets =[  { "folder": "contig_map", "entity": "contig" },
                { "folder": "seq_metadata", "entity": "sequence" } ]

    for target in targets:

        get_delimiters(args.database, target['folder'], target['entity'])
        convert_folder(os.path.join(args.database, target['folder']), target['entity'])


if __name__ == "__main__":
    main()
