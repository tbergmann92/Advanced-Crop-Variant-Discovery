#!/bin/bash -l

# Script to index a reference genome compatible for GATK
# Thomas Bergmann

# load modules if you run on an HPC - adjust
module load gatk4/4.2.5.0--hdfd78af_0
module load bwa/0.7.17--h7132678_9
module load samtools/1.15--h3843a85_0

# set variables
OUT="/PATH/REFERENCE/DIR"
REF="/PATH/REFERENCE_FASTA"
REF_DIC="${OUT}/-REFERENCE_FASTA_NAME-.dict"  # replace with name of your genome fasta!

echo "Creating BWA index"

bwa index \
	-p $OUT/your_reference_genome.fa \  # replace with name of your genome fasta
	 $REF

echo "Done!"
echo ""
echo "Creating GATK sequence dictionary"
echo ""

gatk CreateSequenceDictionary \
	-R $REF \
	-O $REF_DIC

echo "Done!"
echo ""
echo "Indexing with SAMtools"
echo ""

samtools faidx $REF

echo "Done!"

# unload modules
module unload gatk4/4.2.5.0--hdfd78af_0
module unload bwa/0.7.17--h7132678_9
module unload samtools/1.15--h3843a85_0


