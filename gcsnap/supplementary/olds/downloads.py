import os
import requests
from urllib.parse import urljoin
import argparse

def download_file(file_url, file_path):
    with requests.get(file_url, stream=True) as r:
        if r.status_code == 200:
            total_size = int(r.headers.get('content-length', 0))
            block_size = 1024
            downloaded = 0
            with open(file_path, 'wb') as file:
                for data in r.iter_content(block_size):
                    file.write(data)
                    downloaded += len(data)
                    if total_size > 0:
                        percent = 100 * downloaded / total_size
                        print(f"\rDownloading {file_path}: {percent:.2f}% complete", end="")
            print(f"\nDownloaded {file_path}")
        else:
            print(f"Failed to download {file_path}")

def download_files(base_url, download_dir, file_count, file_pattern):
    os.makedirs(download_dir, exist_ok=True)
    for i in range(1, file_count + 1):
        file_name = file_pattern.format(i)
        file_path = os.path.join(download_dir, file_name)
        if os.path.exists(file_path):
            print(f"File {file_name} already exists. Skipping download.")
            continue
        file_url = urljoin(base_url, file_name)
        download_file(file_url, file_path)
    print("All files processed.")

def parse_args():
    """
    Parses command-line arguments for MGnify version and output directory.

    Returns:
        argparse.Namespace: Parsed arguments namespace.
    """

    parser = argparse.ArgumentParser(description="Download MGnify cluster files for a specific version.")
    parser.add_argument("MGnify-version", type=str, help="MGnify peptide database version (e.g., 2024_04)")
    parser.add_argument("output-dir", type=str, help="Output directory to save downloaded files")

    return parser.parse_args()

def main():

    args = parse_args()

    base_url = f"https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/{args.MGnify_version}/"
    
    # contig maps
    download_dir = os.path.join(args.output_dir, "contig_map/")
    file_count = 16  # Manually set the range for the file count
    file_pattern = "mgy_contig_map_{}.tsv.gz"  # Pattern for the file names
    download_files(base_url, download_dir, file_count, file_pattern)

    # seq metadata
    download_dir = os.path.join(args.output_dir, "seq_metadata/")
    file_count = 25  # Manually set the range for the file count
    file_pattern = "mgy_seq_metadata_{}.tsv.gz"  # Pattern for the file names
    download_files(base_url, download_dir, file_count, file_pattern)

    # Download mgy_clusters.fa.gz
    download_dir = args.output_dir
    file_pattern = "mgy_clusters.fa.gz"
    download_file(urljoin(base_url, file_pattern), os.path.join(download_dir, file_pattern))

if __name__ == "__main__":
    main()
