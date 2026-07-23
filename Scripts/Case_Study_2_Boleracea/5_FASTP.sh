#!/bin/bash -l

# run FASTP on raw sequencing data
# Thomas Bergmann

# read the file containing the FASTQ pairs
FASTQ_LIST="PATH/FASTQ_PAIRS_LIST" # pointing to a txt file that contains fastq pairs (tab separated)
QC="PATH/CLEAN/DATA" # where to put the clean data

# get the corresponding line for the current job array index (if run on an HPC with slurm)
READ1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FASTQ_LIST | awk '{print $1}')
READ2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FASTQ_LIST | awk '{print $2}')

# define output file names
OUT1=$(basename "$READ1" .fastq.gz).clean.fastq.gz
OUT2=$(basename "$READ2" .fastq.gz).clean.fastq.gz

SAMPLE=${OUT1%%_*}
ACCESSION=$(echo "$OUT1")

# run fastp
fastp \
	-i "$READ1" \
	-I "$READ2" \
	-o "$QC/$OUT1" \
	-O "$QC/$OUT2" \
	--detect_adapter_for_pe \
	--correction \
	--thread 16

echo ""
echo "Processed $ACCESSION"
echo ""
echo "In: $READ1"
echo "In: $READ2"
echo ""
echo "Out1: $OUT1"
echo "Out2: $OUT2"
echo ""

