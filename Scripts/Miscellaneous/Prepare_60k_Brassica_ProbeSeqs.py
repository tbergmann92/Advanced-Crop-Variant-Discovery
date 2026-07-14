#!/usr/bin/python3

# Written by Thomas Bergmann
# This script converts the full-length sequences from the Brassica 60k SNP array into the 50 bp probe sequences

import sys

# Set global variables
seqDic = {}
count_back=49
counter_fw = 0
counter_rev = 0
directionDic = {}

# try open the multiFASTA
try:
    multiFASTA=open("Brassica_60k_Sequences_Complete.fasta", "r")
except IOError:
    print("Check the input file!")
    sys.exit(1)

# Loop through the multiFASTA and create dictionary and trimming
for line in multiFASTA:
    line = line.rstrip()
    # remove leading and trailing whitespace
    if line.startswith(">"):
        identifier = line[1:]
        seqDic[identifier] = ""
        directionDic[identifier] = ""
    else:
        # add the sequence itself
        start_bracket = line.find('[') # start position of bracket
        end_bracket = line.find(']') + 1 # end position of bracket + 1 to include the bracket
        # check whether you can count backwards
        if start_bracket >= count_back:
            counter_fw += 1
            start_seq = max(0, start_bracket - count_back) # count 49 backwards from start of bracket
            trimmed_seq = line[start_seq:end_bracket] # trim the sequence to 49 length + brackets
            seqDic[identifier] += trimmed_seq # add trimmed seq to dictionary
            directionDic[identifier] += "forward"
        # else count other direction
        else:
            counter_rev += 1
            end_seq = end_bracket + count_back 
            trimmed_seq = line[start_bracket:end_seq] # trim the sequence to 49 length + brackets
            seqDic[identifier] += trimmed_seq # add trimmed seq to dictionary
            directionDic[identifier] += "reverse"


# Close the file after reading
multiFASTA.close()

# IUPAC mapping for forward variants
IUPAC_fw = {
    "R": "[A/G]",
    "Y": "[C/T]",
    "S": "[G/C]",
    "W": "[A/T]",
    "K": "[G/T]",
    "M": "[A/C]"
}

# IUPAC mapping for reverse variants
IUPAC_rev = {
    "R": "[G/A]",
    "Y": "[T/C]",
    "S": "[C/G]",
    "W": "[T/A]",
    "K": "[T/G]",
    "M": "[C/A]"
}

complement = {"A" : "T",
              "T" : "A",
              "C" : "G",
              "G" : "C",
              "R" : "Y",
              "Y" : "R",
              "S" : "S",
              "W" : "W",
              "K" : "M",
              "M" : "K",
              "N" : "N"}


# Function to replace bracketed variants with IUPAC codes
def replace_with_iupac(seq, iupac_fw, iupac_rev):
    # Replace using forward variants
    for iupac_code, variant in iupac_fw.items():
        seq = seq.replace(variant, iupac_code)
    # Replace using reverse variants
    for iupac_code, variant in iupac_rev.items():
        seq = seq.replace(variant, iupac_code)
    
    return seq

# Apply IUPAC replacement for each sequence in seqDIc
for identifier, sequence in seqDic.items():
    #print(f"Original sequence for {identifier}: {sequence}") 
    # Replace the brackets with the iupac code for forward direction
    seqDic[identifier] = replace_with_iupac(sequence, IUPAC_fw, IUPAC_rev)
    #print(f"Converted sequence for {identifier}: {seqDic[identifier]}")
    # Replace the brackets with the iupac code for reverse direction
    if directionDic[identifier] == "reverse":
        #print(f"Original sequence for {identifier}: {sequence}") 
        seqDic[identifier] = seqDic[identifier].translate(str.maketrans(complement))[::-1]
        #print(f"Converted sequence for {identifier}: {seqDic[identifier]}")

# Write to FASTA
with open("Brassica_60K_50bp_Seqs.fasta", "w") as fasta_file:
    for identifier, sequence in seqDic.items():
        fasta_file.write(f">{identifier}\n{sequence}\n")

print("FASTA format output written to Brassica_60k_50bp_ProbeSeqs.fasta")
print("Number of forward sequences: {}".format(counter_fw))
print("Number of reverse sequences: {}".format(counter_rev))

