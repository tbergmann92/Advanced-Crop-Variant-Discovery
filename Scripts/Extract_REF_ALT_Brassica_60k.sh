#!/bin/bash

# Script to extract REF and ALT alleles for each SNP from the Brassica 60k SNP array

awk '
/^>/ { 
    snp_id = substr($0, 2); next 
}
/\[[ACGTN]/ {
    match($0, /\[([ACGTN])\/([ACGTN])\]/, m);
    if (m[1] && m[2]) print snp_id "\t" m[1] "\t" m[2];
}
' Brassica_60k_Sequences_complete.fasta > Brassica_60k_REF_ALT_Alleles.tsv


