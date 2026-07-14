#!/bin/bash -l

# load modules if available
module load gatk4/4.2.5.0--hdfd78af_0
module load bwa/0.7.17--h7132678_9
module load samtools/1.15--h3843a85_0

# read the file containing the FASTQ pairs
OUT="/PATH/TO/REFERENCE/DIR"
REF="/PATH/TO/REFERENCE/FILE/DarmorV10_Chromosomes_Only.fa"
REF_DIC="${OUT}/DarmorV10_Chromosomes_Only.dict"

echo "Creating BWA index"

bwa index \
	-p $OUT/DarmorV10_Chromosomes_Only.fa \
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

module unload gatk4/4.2.5.0--hdfd78af_0
module unload bwa/0.7.17--h7132678_9
module unload samtools/1.15--h3843a85_0


