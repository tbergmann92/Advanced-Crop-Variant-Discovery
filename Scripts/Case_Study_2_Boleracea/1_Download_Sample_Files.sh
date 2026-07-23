#!/bin/bash -l

# Script to download SRA files from NCBI
# Written by Thomas Bergmann

set -euo pipefail

# Set variables
TEMP_DIR="/PATH/TMP/DIR" # a temp dir for collecting data
SRA_DIR="/PATH/SRA/DIR" # path to your sra dir
FASTQ_DIR="/PATH/DATA/DIR" # path to where to store fastq data
ACCESSIONS="/PATH/ACCESSIONS_TXT" # list of sra accessions to download

# Read sample into an array
mapfile -t SAMPLES < "$ACCESSIONS"

# Get sample ID for this task if you use slurm - parallise the download
SAMPLE_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# Prepare directories
mkdir -p $FASTQ_DIR
mkdir -p $SRA_DIR

echo ""
echo "$(date) Downloading --> $SAMPLE_ID"
echo ""

prefetch "$SAMPLE_ID" -p --max-size u

echo ""
echo "$(date) Finished!"
echo ""
echo "Validating file ..."
echo ""

vdb-validate "$SRA_DIR/${SAMPLE_ID}.sra"

echo ""
echo "$(date) - Done!"
echo ""
echo "Extracting FASTQ from SRA - $SAMPLE_ID"
echo ""

fasterq-dump -e 32 -p -t "$TEMP_DIR" --outdir "$FASTQ_DIR" "${SRA_DIR}/${SAMPLE_ID}.sra"

echo "$(date) - Done!"
echo ""
echo "Compressing the FASTQ"
echo ""

pigz -p 64 "$FASTQ_DIR/${SAMPLE_ID}_1.fastq"
pigz -p 64 "$FASTQ_DIR/${SAMPLE_ID}_2.fastq"

echo ""
echo "$(date) - Done!"
echo ""
