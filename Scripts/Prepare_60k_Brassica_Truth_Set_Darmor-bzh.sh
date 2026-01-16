#!/bin/bash -l

# Written by Thomas Bergmann

# Required inputs:
# $BAM > Sorted and indexed BAM of Brassica 60k ProbeSeqs mapped against Darmor-bzh
# $ALELLES_TSV > Brassica_60k_REF_ALT_Alleles.tsv

# Write VCF header
VCF="$OUT_DIR/known_sites_DarmorV10.vcf"
{
echo "##fileformat=VCFv4.2"
echo "##contig=<ID=A01,length=32958928>"
echo "##contig=<ID=A02,length=33432960>"
echo "##contig=<ID=A03,length=39685748>"
echo "##contig=<ID=A04,length=23101715>"
echo "##contig=<ID=A05,length=42112164>"
echo "##contig=<ID=A06,length=45146386>"
echo "##contig=<ID=A07,length=29390523>"
echo "##contig=<ID=A08,length=26309499>"
echo "##contig=<ID=A09,length=53549826>"
echo "##contig=<ID=A10,length=20778245>"
echo "##contig=<ID=C01,length=48239358>"
echo "##contig=<ID=C02,length=62297340>"
echo "##contig=<ID=C03,length=73669886>"
echo "##contig=<ID=C04,length=65837619>"
echo "##contig=<ID=C05,length=56382805>"
echo "##contig=<ID=C06,length=50218839>"
echo "##contig=<ID=C07,length=55656957>"
echo "##contig=<ID=C08,length=41681856>"
echo "##contig=<ID=C09,length=66465249>"
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
