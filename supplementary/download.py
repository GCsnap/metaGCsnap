import os
import re
import requests
import argparse
from urllib.parse import urljoin

def list_files_at_url(url, pattern):
    """
    Fetches a URL and scrapes it for href links matching a regex pattern.
    
    Args:
        url (str): The URL of the directory page to scrape.
        pattern (str): A regex pattern to match filenames.
    
    Returns:
        list: A list of matching filenames, or None if the request fails.
    """
    print(f"Checking for files at {url}...")
    try:
        with requests.get(url) as r:
            if r.status_code == 200:
                # Find all links that match the file pattern
                # e.g., href="mgy_contig_map_1.tsv.gz"
                filenames = re.findall(pattern, r.text)
                if not filenames:
                    print("Warning: No files found matching the pattern.")
                    return []
                # Return a list of unique filenames
                return sorted(list(set(filenames)))
            else:
                print(f"Failed to list files. Status code: {r.status_code}")
                return None
    except requests.RequestException as e:
        print(f"Error connecting to {url}: {e}")
        return None

def download_file(file_url, file_path):
    """
    Downloads a single file with a progress bar.
    """
    print(f"Requesting {file_url}...")
    try:
        with requests.get(file_url, stream=True, timeout=30) as r:
            r.raise_for_status() # Will raise an HTTPError for bad responses
            total_size = int(r.headers.get('content-length', 0))
            block_size = 1024
            downloaded = 0
            
            with open(file_path, 'wb') as file:
                for data in r.iter_content(block_size):
                    file.write(data)
                    downloaded += len(data)
                    if total_size > 0:
                        percent = 100 * downloaded / total_size
                        # Use a carriage return to keep the progress bar on one line
                        print(f"\rDownloading {os.path.basename(file_path)}: {percent:.2f}% complete", end="")
            
            # Print a newline after download is complete
            print(f"\nDownloaded {os.path.basename(file_path)}")
        
    except requests.HTTPError as e:
        print(f"\nFailed to download {os.path.basename(file_path)}. HTTP Error: {e}")
    except requests.RequestException as e:
        print(f"\nFailed to download {os.path.basename(file_path)}. Error: {e}")

def download_taxonomy_db(download_job):
    """
    Downloads a taxonomy database file.
    """
    print(f"Downloading taxonomy database from {download_job['url']}...")

    kraken_gtdb_url = 'https://genome-idx.s3.amazonaws.com/kraken/k2_gtdb_genome_reps_20250609.tar.gz'

    os.makedirs(download_job['dir'], exist_ok=True)

    file_path = os.path.join(download_job['dir'], download_job['prefix'])

    download_file(kraken_gtdb_url, file_path)

def argument_parsing():

    parser = argparse.ArgumentParser(
        description="Download MGnify dataset files from the EBI FTP server."
    )

    parser.add_argument(
        "--MGnify-url",
        type=str,
        default="https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/",
        help="The URL of the MGnify dataset. Default: https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/"
    )

    parser.add_argument(
        "--out-dir",
        type=str,
        required=True,
        help="The base output directory to store the 'MGnify' folder. (Required)"
    )

    parser.add_argument(
        "--MGnify-version",
        type=str,
        required=True,
        help="The MGnify dataset version to download (e.g., '2024_04'). Default: 2024_04"
    )

    parser.add_argument(
        "--protein-db",
        action="store_true",
        default=False,
        help="Flag to use the local protein database directory (default: False). If present, set to True."
    )

    parser.add_argument(
        "--taxonomy-db",
        action="store_true",
        default=False,
        help="Flag to use the local protein database directory (default: False). If present, set to True."
    )

    parser.add_argument(
        "--taxonomy-db-dir",
        default=None,
        help="The local taxonomy database directory. Default: None"
    )
    
    return parser.parse_args()

def get_jobs(args):

    # Define the download jobs. We will filter the master list for each job.
    download_jobs = dict()
    download_jobs['MGnify'] = [
        {
            "dir": os.path.join(args.out_dir, f"MGnify_{args.MGnify_version}", "contig_map"),
            "prefix": "mgy_contig_map_"
        },
        {
            "dir": os.path.join(args.out_dir, f"MGnify_{args.MGnify_version}", "seq_metadata"),
            "prefix": "mgy_seq_metadata_"
        },
    ]

    if args.protein_db:
        print('Including protein files')
        download_jobs['MGnify'].append( {
            "dir": os.path.join(args.out_dir, f"MGnify_{args.MGnify_version}"), # Will be saved in the root MGnify dir
            "prefix": "mgy_clusters.fa.gz"
        } )

    if args.taxonomy_db:
        
        download_jobs['GTDB'] = { "dir": args.taxonomy_db_dir, # Will be saved in the taxonomy db dir
                                  "prefix": "k2_gtdb_genome_reps_20250609.tar.gz" }
    
    return download_jobs

def main():

    args = argument_parsing()

    os.makedirs(os.path.join(args.out_dir, f"MGnify_{args.MGnify_version}"), exist_ok=True)

    base_url = os.path.join(args.MGnify_url, args.MGnify_version , '/')
    
    # This single regex will find all files we are interested in.
    # It looks for links (href) to files starting with "mgy_" and ending in ".gz"
    file_regex_pattern = r'href="(mgy_[^"]+\.gz)"'

    # Get the complete list of relevant files from the server ONCE
    all_remote_files = list_files_at_url(base_url, file_regex_pattern)
    print(base_url)

    if all_remote_files is None:
        print(f"Could not fetch file list from server: {base_url}. Check version and connection. Exiting.")
        return

    print(f"\nFound {len(all_remote_files)} total matching files on server.")

    # Define the download jobs. We will filter the master list for each job.
    download_jobs = get_jobs(args)

    for job in download_jobs['MGnify']:
        
        download_dir = job["dir"]
        file_prefix = job["prefix"]
        
        # Ensure the target directory exists
        os.makedirs(download_dir, exist_ok=True)
        
        # Filter the master file list based on the job's prefix
        files_for_this_job = [f for f in all_remote_files if f.startswith(file_prefix)]
        
        if not files_for_this_job:
            print(f"\nNo remote files found for prefix '{file_prefix}'.")
            continue

        print(f"\n--- Starting Job: '{file_prefix}' ---")
        print(f"Found {len(files_for_this_job)} files to download to {download_dir}")

        for file_name in files_for_this_job:
            file_path = os.path.join(download_dir, file_name)
            
            # Check if file already exists
            if os.path.exists(file_path):
                print(f"File {file_name} already exists. Skipping download.")
                continue
            
            # Construct the full URL and download
            file_url = urljoin(base_url+'/', file_name)
            download_file(file_url, file_path)

    print(base_url)
    if args.taxonomy_db:

        download_taxonomy_db(download_jobs['GTDB'])

    print("\nAll download jobs processed.")

if __name__ == "__main__":
    main()