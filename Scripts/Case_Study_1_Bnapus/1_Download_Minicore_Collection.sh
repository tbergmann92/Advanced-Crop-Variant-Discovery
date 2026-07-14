#!/bin/bash -l

# Fail fast on errors and unset variables.
set -euo pipefail

# Paths
FASTQ_DIR="/PATH/TO/DATA/DIR"
ACCESSIONS="/PATH/TO/DATA/LIST/Download_List_Minicore.txt"
LOG_DIR="/PATH/TO/DATA/DIR/"

mkdir -p "$LOG_DIR"
mkdir -p "$FASTQ_DIR"

# Get the correct line from the file - array job for SLURM
TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set}
LINE=$(awk -v id="$TASK_ID" 'NR==id' "$ACCESSIONS")

# Parse columns
IFS=$'\t' read -r ACCESSION FILE1 MD5_1 URL1 FILE2 MD5_2 URL2 <<< "$LINE"

# Clean filenames
FILE1=$(echo "$FILE1" | tr -d '\r' | cut -d' ' -f1)
FILE2=$(echo "$FILE2" | tr -d '\r' | cut -d' ' -f1)
URL1=$(echo "$URL1" | tr -d '\r')
URL2=$(echo "$URL2" | tr -d '\r')

echo "Processing $ACCESSION"

# Download function

download_and_check () {
	local url=$1
	local outfile=$2
	local expected_md5=$3
	local calculated_md5
	
	if [[ -f "$outfile" ]]; then
		echo "File exists, checking MD5: $outfile"
		local calculated_md5
		read -r calculated_md5 _ < <(md5sum "$outfile")

		if [[ "$calculated_md5" == "$expected_md5" ]]; then
			echo "MD5 OK, skipping download: $outfile"
			return
		else
			echo "MD5 mismatch, re-downloading: $outfile"
			rm -f "$outfile"
		fi
	fi
	
	echo "$(date) Downloading $outfile"
	wget -c --tries=3 --waitretry=5 -O "$outfile" "$url"

	echo "Checking MD5 for $outfile"
	read -r calculated_md5 _ < <(md5sum "$outfile")

	if [[ "$calculated_md5" != "$expected_md5" ]]; then
		echo "ERROR: MD5 mismatch for $outfile"
		echo "Expected: $expected_md5"
		echo "Got: $calculated_md5"
		exit 1
	else
		echo "MD5 OK for $outfile"
	fi
}

# Run downloads
download_and_check "$URL1" "$FASTQ_DIR/$FILE1" "$MD5_1" \
	> "$LOG_DIR/${ACCESSION}_f1.log" 2>&1 &
download_and_check "$URL2" "$FASTQ_DIR/$FILE2" "$MD5_2" \
	> "$LOG_DIR/${ACCESSION}_r2.log" 2>&1 &
wait

echo "Finished $ACCESSION"





 
