#!/bin/bash
#SBATCH --job-name=mmseqs_createdb
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --qos=1week
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G

set -euo pipefail

########################
# 1. Parse arguments
########################

OUTROOT="$1"
INPUT_FASTA="$1/mgy_clusters.fa.gz"
DB_NAME="DB"

########################
# 2. Layout of directories
########################

DB_DIR="$OUTROOT/mmseqsDBs/mgyc"
TMP_DIR="$OUTROOT/tmp"

mkdir -p "$DB_DIR" "$TMP_DIR"

DB_PREFIX="$DB_DIR/$DB_NAME"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

echo "===== JOB INFO ====="
echo "Input FASTA    : $INPUT_FASTA"
echo "Output root    : $OUTROOT"
echo "DB prefix      : $DB_PREFIX"
echo "Tmp directory  : $TMP_DIR"
echo "Threads        : $THREADS"
echo "===================="
echo

########################
# 3. Activate mamba env
########################

if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi

mamba activate metaGCsnap
echo "Activated mamba env: metaGCsnap"
echo

########################
# 4. Build DB
########################

echo "Running: mmseqs createdb ..."
mmseqs createdb \
    "$INPUT_FASTA" \
    "$DB_PREFIX" \
    "$TMP_DIR" \
    --threads "$THREADS" \
    --createdb-mode 1

echo "createdb done."
echo

########################
# 5. Create index (for faster search)
########################

echo "Running: mmseqs createindex ..."
mmseqs createindex \
    "$DB_PREFIX" \
    "$TMP_DIR" \
    --threads "$THREADS"

STATUS=$?
echo "createindex exit status: $STATUS"

if [ "$STATUS" -ne 0 ]; then
    echo "ERROR: mmseqs createindex failed. Tmp dir kept at: $TMP_DIR"
    exit "$STATUS"
fi

# Optionally clean tmp after everything:
# rm -rf "$TMP_DIR"

echo "All done."
echo "DB + index: $DB_DIR"
echo "Tmp files : $TMP_DIR"
echo "SLURM logs: ./logs"
