#!/bin/bash -l

# script to identify and transform Brassica SNP probes into a VCF
# Thomas Bergmann

# Load the necessary modules
module load bwa/0.7.17--h7132678_9
module load samtools/1.15--h3843a85_0
module load bcftools/1.15--haf5b3da_0
module load gatk4/4.2.5.0--hdfd78af_0

# Define variables and paths
BWA_INDEX="/PATH/BWA_INDEX" # path to your reference with BWA index
OUT_DIR="/PATH/OUTPUT" # where to put your results
ALLELES_TSV="$OUT_DIR/array_data/Brassica_60k_REF_ALT_Alleles.tsv" # available in 'Resources'on my github
PROBE_SEQS="$OUT_DIR/array_data/Brassica_60k_50bp_ProbeSeqs.fasta" # available in 'Resources'on my github

echo -e "\nStarting mapping...\n"

bwa aln \
	-n 0.01 \
	-o 0 \
	-l 32 \
	-k 2 \
	-t 16 \
	"$BWA_INDEX" \
	"$PROBE_SEQS" | \
bwa samse \
	"$BWA_INDEX" \
	- \
	"$PROBE_SEQS" | \
samtools sort \
	-@ 16 \
	-o "${OUT_DIR}/array_data/Brassica_60k_vs_W03.sorted.bam"

samtools index "${OUT_DIR}/array_data/Brassica_60k_vs_W03.sorted.bam"

# Generate mapping statistics
echo -e "\nGenerating mapping statistics...\n"

samtools flagstat "$OUT_DIR/array_data/Brassica_60k_vs_W03.sorted.bam" > "$OUT_DIR/array_data/Brassica_60k_vs_W03.flagstat.txt"
samtools stats "$OUT_DIR/array_data/Brassica_60k_vs_W03.sorted.bam" > "$OUT_DIR/array_data/Brassica_60k_vs_W03.stats.txt"
samtools idxstats "$OUT_DIR/array_data/Brassica_60k_vs_W03.sorted.bam" > "$OUT_DIR/array_data/Brassica_60k_vs_W03.idxstats.txt"

echo -e "\nMapping done!\n"
echo "Filtering BAM and generating VCF..."

# Write VCF header for B. oleracea 'W03'
# Header needs to be adjusted to the corresponding reference genome !!!
# Use samtools faidx or other tools to get these stats and adjust the header correspondingly
VCF="$OUT_DIR/array_data/known_sites_W03.vcf"
{
echo "##fileformat=VCFv4.2"
echo "##contig=<ID=C01,length=58570197>"
echo "##contig=<ID=C02,length=72559525>"
echo "##contig=<ID=C03,length=80735250>"
echo "##contig=<ID=C04,length=71180230>"
echo "##contig=<ID=C05,length=63147374>"
echo "##contig=<ID=C06,length=49313063>"
echo "##contig=<ID=C07,length=64564754>"
echo "##contig=<ID=C08,length=57383795>"
echo "##contig=<ID=C09,length=76892030>"
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
} > "$VCF"

# Extract, filter and write known sites into VCF
samtools view -h -F 0x904 "$OUT_DIR/array_data/Brassica_60k_vs_W03.sorted.bam" | \
awk -v ref_tab="$ALLELES_TSV" '
BEGIN {
  # Load REF and ALT from table
  while ((getline < ref_tab) > 0) {
    split($0, a, "\t")
    ref[a[1]] = a[2]
    alt[a[1]] = a[3]
  }
}
# Skip headers
/^@/ { next }

($0 ~ /X0:i:1/) && ($0 !~ /XA:Z:/) && $6 == "50M" {

  id = $1; flag = $2; chrom = $3; pos = $4; mapq = $5; seq = $10

  is_reverse = and(flag, 16)
  
  if (is_reverse && substr(seq, 1, 1) == "N") {
    snp_pos = pos
  } else if (!is_reverse && substr(seq, length(seq), 1) == "N") {
    snp_pos = pos + 49
  } else {
    next
  }

  if (id in ref) {
    print chrom, snp_pos, id, ref[id], alt[id], mapq, "PASS", "."
  }
}
' OFS='\t' >> "$VCF"

# Sort and index VCF
bcftools sort -Oz "$VCF" -o "${OUT_DIR}/known_sites_W03.sorted.vcf.gz"
bcftools index -t "${OUT_DIR}/known_sites_W03.sorted.vcf.gz"

echo -e "\nVCF written to $VCF"