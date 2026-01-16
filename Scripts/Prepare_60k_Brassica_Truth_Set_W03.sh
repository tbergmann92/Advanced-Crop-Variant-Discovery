#!/bin/bash -l

# Written by Thomas Bergmann

# Required inputs:
# $BAM > Sorted and indexed BAM of Brassica 60k ProbeSeqs mapped against W03
# $ALELLES_TSV > Brassica_60k_REF_ALT_Alleles.tsv

# Write VCF header
VCF="$OUT_DIR/known_sites_W03.vcf"
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
samtools view -h -q 40 -F 0x104 "$BAM" | \
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

# Process reads with exactly 50M CIGAR
$6 == "50M" {
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

echo -e "\nVCF written to $VCF"
