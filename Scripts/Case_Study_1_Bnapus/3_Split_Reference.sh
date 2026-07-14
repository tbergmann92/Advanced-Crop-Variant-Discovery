#!/usr/bin/bash -l

# Fail fast on errors and unset variables.
set -euo pipefail

# Load the necessary modules if available
module load gatk4/4.2.5.0--hdfd78af_0

REF="/PATH/TO/REFERENCE/DarmorV10_Chromosomes_Only.fa"
OUTDIR="/PATH/TO/REFERENCE/INTERVALS_DIR"

mkdir -p "$OUTDIR"

gatk SplitIntervals \
	-R "$REF" \
	-scatter 100 \
	--subdivision-mode INTERVAL_SUBDIVISION \
	-O "$OUTDIR"

# Unload the necessary modules
module unload gatk4/4.2.5.0--hdfd78af_0



