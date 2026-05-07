# metaGCsnap

metaGCsnap is an extended version of [GCsnap 2.0](https://github.com/GCsnap/gcsnap2desktop) that supports multiple genomic-context data providers in a single run: **NCBI**, **UniProt**, **MGnify**, and a **local database**. Targets from different providers are merged into a unified genomic-context view, enabling direct comparison of genomic contexts across resources. An overview of metaGCsnap workflow is shown below.


![](./metaGCsnap.png)

---

## Table of Contents

1. [Installation](#installation)
2. [Providers setup](#setting-up-the-mgnify-database)
4. [Usage](#usage)
5. [Configuration](#configuration)
6. [Credits](#credits)

---

## 1] Installation

**conda or mamba must be installed before running the install script.**

### Step 1 — clone and create the base environment

```bash
git clone https://github.com/GCsnap/metaGCsnap.git
cd metaGCsnap
bash install_providers.sh --base
conda activate gcsnap
```

This creates the `gcsnap` conda environment with all core dependencies (Python 3.11, MMseqs2, Biopython, Bokeh, Matplotlib, Networkx, PaCMAP, Scikit-learn, Pandas, Rich, Requests, and more) and registers the `GCsnap` CLI. No data provider is available yet.

### Step 2 — install provider dependencies

Run the script again with the flag for the provider(s) you need:

```bash
bash install_providers.sh --ncbi      # NCBI / UniProt provider
bash install_providers.sh --mgnify    # MGnify metagenomics provider
bash install_providers.sh --local     # local database provider
bash install_providers.sh --complete  # all three providers
```

Providers can be added at any time after the base install. The base environment is preserved and only the provider-specific packages are added.

### Step 3 — configure provider configuration fields

Open `config.yaml` and fill in the paths required by the providers you installed. For a detailed description of the fields, see the config section below:

| Provider | Keys to set |
|----------|-------------|
| NCBI     | `ncbi-user-email`, `ncbi-api-key` |
| MGnify   | `MGnify-path`, `MGnify-proxies` (optional) |
| local    | `gff-path`, `db-path` |

For instructions on how to get the ncbi-api, see https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys.
### Windows note

MMseqs2 has no Windows conda package. Install without `--local` and download the static binary from https://mmseqs.com/latest/. Pass the path via `--mmseqs-executable-path` when running GCsnap.

---

## 2] Providers setup

To access the data of each repositroy, you are required to run some scripts.  

## 2.1] Setting Up the MGnify Database

Before using the MGnify provider, you need to download the relevant MGnify database files.

### Choose a version

- `2024_04` — most recent version; best for exhaustive local searches. Matches the [MGnify online phmmer search](https://www.ebi.ac.uk/metagenomics/sequence-search/search/phmmer). Requires considerable storage (~103 GB after conversion).
- `2023_02` — corresponds to the [ESM Atlas](https://esmatlas.com/about). Use this to investigate have a correspondence to protein structures.

### Download the raw files

As a concrete example, we will consider the most recent MGnify version, which is 2024_04

```bash
python3 gcsnap/supplementary/MGnify/download.py \
    --MGnify-version 2024_04 \
    --out-dir your/out/dir
```

Both `--MGnify-version` and `--out-dir` are required. Given the limitations of MGnify API we will need to host locally a set of files to match different IDs. This creates `your/out/dir/MGnify_2024_04/` with the relevant data. Separating this step from the conversion step is useful when working on HPC systems with distinct download and compute nodes. This step is sufficient to run metaGCsnap, once the target MGnify IDs are provided. Such IDs can be obtained from the putput of the online phmmer search. For a concrete example, you can see the tutorial.

### Convert to Parquet

To speed up ID searches, convert the `.tsv.gz` files to `.parquet`:

```bash
python3 gcsnap/supplementary/MGnify/make_MGnify.py \
    --database your/out/dir/MGnify_2024_04
```

Once the conversion is complete, you can delete the original `.tsv.gz` files to save space. (we will take care of ths aspect later)

### Expected database layout

After setup, `your/out/dir/MGnify_2024_04/` should look like:

```
MGnify_2024_04/ (~103 GB)
├── contig_map/           (~32 GB)
│   ├── delimiters.csv
│   └── mgy_contig_map_*.parquet
└── seq_metadata/         (~71 GB)
    ├── delimiters.csv
    └── mgy_seq_metadata_*.parquet
```

The folder structure is the same regardless of the MGnify version used.

---

## 2.2] Setting Up a Local Database

The local provider works with a pre-built SQLite database of assemblies and sequences. Setup scripts are in `gcsnap/supplementary/local/`.

```bash
python3 gcsnap/supplementary/local/db_create_assemblies.py --config config.yaml
```

Set `gff-path` and `db-path` in `config.yaml` to point to your GFF annotation folder and the resulting database, respectively.

---

## Provider-specific config edits

Depending on the installed providers, you will have to specify 

| Provider | Mandatory field | Description |
|----------|----------------|-------------|
| MGnify | `MGnify-path` | Path to your local `MGnify_<version>/` folder |
| MGnify | `MGnify-proxies` | http/s proxies for API calls |
| NCBI | `ncbi-user-email` | Your email address, required by NCBI Entrez |
| NCBI | `ncbi-api-key` | Your NCBI API key — raises rate limit from 3 to 10 req/s |
| local | `gff-path` | Path to your folder of `.gff.gz` annotation files |
| local | `db-path` | Path to your GCsnap SQLite database |

## Usage
GCsnap requires at least one provider target file. Provider flags can be combined freely to run a mixed NCBI + MGnify job.

```bash
GCsnap --ncbi-targets    path/to/ncbi_ids.txt
GCsnap --mgnify-targets  path/to/mgyp_ids.txt
GCsnap --local-targets   path/to/local_ids.txt

# Combined run
GCsnap --ncbi-targets path/to/ncbi_ids.txt --mgnify-targets path/to/mgyp_ids.txt
```

Each target file is a plain-text file with one identifier per line (`#` lines are ignored).

All optional arguments can be set in `config.yaml` or passed directly on the CLI (e.g. `--n-cpu 8`). CLI values take precedence when `overwrite-config: true`.

```bash
GCsnap --help   # full list of arguments and current defaults
```

### Resume from a previous run

If a run is interrupted, re-running the same command will automatically skip providers and steps whose output files are already present on disk.

---


## Advanced topics

### Taxonomic assignment

By default, metaGCsnap clusters contigs using SourMash (ANI estimation). To enable actual taxonomic assignment, download the GTDB Kraken2 index (requires ~0.5 TB):

```bash
python3 gcsnap/supplementary/MGnify/download.py \
    --MGnify-version 2024_04 \
    --out-dir your/out/dir \
    --taxonomy-db \
    --taxonomy-db-dir your/taxdb/dir
```

To the include taxonomic profiling in metaGCsnap workflow, edit the following configuration flags:

| Provider | Mandatory field | Description |
|----------|----------------|-------------|
| MGnify | `genome-classification` | taxonomy
| MGnify | `kraken-path` | Path to your GTDB Kraken2 index |

The taxonomy and protein search directories can be stored separately. Set `kraken-path` in `config.yaml` to `your/taxdb/dir` to activate Kraken2-based taxonomy.

### Local sequence search against MGnify proteins
By default, metaGCsnap identifies MGnify targets by matching MGYP identifiers against the locally hosted Parquet index. For large-scale or custom queries that require searching by sequence rather than by ID, you can build a local MMseqs2 index against the full MGnify protein catalogue.

First, download the protein FASTA (approx. 74 GB for 2024_04) by re-running the download script with the --local-proteins flag:

```bash
python3 gcsnap/supplementary/MGnify/download.py \
    --MGnify-version 2024_04 \
    --out-dir your/out/dir \
    --local-proteins
```

This command will skip the download of already present files. Then submit the provided SLURM script to build the MMseqs2 index (edit resource settings as needed for your cluster):

sbatch gcsnap/supplementary/MGnify/make_mmseqs.sh your/out/dir
This generates an index at your/out/dir/mmseqsDBs/mgyc. For most users the 90% identity cluster representative set is recommended. Note that local sequence search can require up to ~350 GB of RAM.


## Configuration

Edit `config.yaml` to set default values. Key options:

| Option | Description | Default |
|--------|-------------|---------|
| `out-label` | Output directory name (defaults to input filename) | `default` |
| `n-cpu` | Number of CPU cores | `4` |
| `n-flanking5` | Flanking genes on the 5′ end | `4` |
| `n-flanking3` | Flanking genes on the 3′ end | `4` |
| `collect-only` | Collect genomic contexts only, skip comparisons | `false` |
| `exclude-partial` | Exclude partial genomic context blocks | `true` |
| `max-evalue` | E-value cutoff for homology (protein families) | `0.001` |
| `min-coverage` | Minimum alignment coverage for homology | `0.7` |
| `genome-classification` | Contig classification method: `binning` or `taxonomy` | `binning` |
| `MGnify-path` | Path to the MGnify database folder | — |
| `kraken-path` | Path to the Kraken2 GTDB index (for taxonomy) | — |
| `gff-path` | Path to GFF annotation folder (local provider) | — |
| `db-path` | Path to local GCsnap SQLite database | — |
| `tmp-folder` | Temporary folder for MMseqs2 files | `./tmp` |
| `operon-cluster-advanced` | Enable advanced PaCMAP operon clustering | `false` |
| `max-family-freq` | Max family frequency for advanced clustering | `20` |
| `min-family-freq` | Min family frequency for advanced clustering | `2` |
| `min-family-freq-accross-contexts` | Min family frequency within a context type | `30` |
| `n-max-operons` | Max number of top-populated operon types shown | `30` |
| `interactive` | Generate interactive HTML output | `true` |
| `genomic-context-cmap` | Matplotlib colormap for syntenic blocks | `Spectral` |
| `gc-legend-mode` | Legend mode: `species` or `ncbi_code` | `species` |
| `out-format` | Figure format: `png`, `svg`, or `pdf` | `png` |
| `min-coocc` | Minimum co-occurrence to connect two genes in graphs | `0.3` |
| `sort-mode` | Context sort mode: `taxonomy`, `as_input`, `tree`, `operon`, `operon cluster` | `taxonomy` |
| `overwrite-config` | Allow CLI values to overwrite `config.yaml` | `false` |
| `mmseqs-executable-path` | Path to MMseqs2 binary if not in conda env | — |
| `foldseek-executable-path` | Path to Foldseek binary if not in conda env | — |
| `sourmash-executable-path` | Path to Sourmash binary if not in conda env | — |

For arguments ending in `-path` to executables, run `which mmseqs` (or similar) and paste the result into `config.yaml`.

---

## Credits

metaGCsnap is built on top of [GCsnap and GCsnap 2.0](https://www.sciencedirect.com/science/article/pii/S0022283621001443). The original desktop and cluster versions can be found at [gcsnap2desktop](https://github.com/GCsnap/gcsnap2desktop) and [gcsnap2cluster](https://github.com/GCsnap/gcsnap2cluster).

If you use metaGCsnap, please cite the original GCsnap paper:

> J. Pereira, GCsnap: interactive snapshots for the comparison of protein-coding genomic contexts, *J. Mol. Biol.* (2021) 166943. https://doi.org/10.1016/j.jmb.2021.166943

metaGCsnap is being developed at the Biozentrum of the University of Basel (Schwede group).
