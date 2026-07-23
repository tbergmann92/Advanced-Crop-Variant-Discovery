#!/usr/bin/bash -l

# script to split a reference for GATK 
# Thomas Bergmann

set -euo pipefail

# Load the necessary modules
module load gatk4/4.2.5.0--hdfd78af_0

REF="/PATH/REFERENCE.fa" # your reference genome
OUTDIR="/PATH/REFERENCE/intervals" # where to put the intervals 

mkdir -p "$OUTDIR"

gatk SplitIntervals \
	-R "$REF" \
	-scatter 100 \
	--subdivision-mode INTERVAL_SUBDIVISION \
	-O "$OUTDIR"

# Unload the necessary modules
module unload gatk4/4.2.5.0--hdfd78af_0
