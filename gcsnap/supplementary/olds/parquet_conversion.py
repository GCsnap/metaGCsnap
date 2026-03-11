#!/usr/bin/env python3
import os
import argparse
import glob
import pandas as pd

def convert_tsv_gz_to_parquet(input_file, output_file, entity):
    """
    Reads a gzipped TSV file and writes it as a Parquet file.
    Assumes the TSV file has no header and two columns:
    - entity_id (first column)
    - data (second column)
    """
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

def convert_folder(args):
    """
    Converts all .tsv.gz files in the given folder to Parquet.
    Output Parquet files will be stored in the same folder.
    """
    # Get list of .tsv.gz files in the folder.
    tsv_files = glob.glob(os.path.join(args.folder, "*.tsv.gz"))
    if not tsv_files:
        print("No .tsv.gz files found in the specified folder.")
        return

    for tsv_file in tsv_files:
        base = os.path.splitext(os.path.basename(tsv_file))[0]
        # Remove the extra extension from '.tsv.gz' to just get the base name.
        if base.endswith('.tsv'):
            base = base[:-4]
        output_file = os.path.join(args.folder, f"{base}.parquet")
        convert_tsv_gz_to_parquet(tsv_file, output_file, args.entity)

def main():

    # python3 parquet_conversion.py --folder ../../../data/MGnify_2023_02/contig_map --entity contig
    # python3 parquet_conversion.py --folder ../../../data/MGnify_2023_02/seq_metadata --entity sequence

    parser = argparse.ArgumentParser(
        description="Convert all .tsv.gz files in a folder to Parquet format."
    )
    parser.add_argument("--folder", help="Folder containing .tsv.gz files to convert.")
    parser.add_argument("--entity", help="Folder containing .tsv.gz files to convert.")
    args = parser.parse_args()

    if not os.path.isdir(args.folder):
        print(f"The folder {args.folder} does not exist or is not a directory.")
        return

    convert_folder(args)

if __name__ == "__main__":
    main()
