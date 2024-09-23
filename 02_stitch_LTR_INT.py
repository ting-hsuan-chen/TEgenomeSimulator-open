import os
import sys
from Bio import SeqIO

# Load files
original_te_fa = str(sys.argv[1])
ltr_retro_fa = str(sys.argv[2])
out_file_name = str(sys.argv[3])
ltr_te_seq = SeqIO.to_dict(SeqIO.parse(ltr_retro_fa ,"fasta"))
te_seq = SeqIO.to_dict(SeqIO.parse(original_te_fa ,"fasta"))

with open('header_LTR_trim.list', 'r') as file:
    ltr_list = [line.strip() for line in file]

with open('header_LTR.list', 'r') as file:
    ltr_full = [line.strip() for line in file]

with open('header_INT_trim.list', 'r') as file:
    int_list = [line.strip() for line in file]

with open('header_INT.list', 'r') as file:
    int_full = [line.strip() for line in file]

# Find matched TE family name
matched = {}
for i in range(len(ltr_list)):
    ltr_name = ltr_list[i]
    for j in range(len(int_list)):
        int_name = int_list[j]
        if ltr_name == int_name:
           matched[ltr_name] = (i, j)   
        
# Stitch LTR-INT-LTR and output fasta file
stitched_out = open("ltr_retro_stitched.fasta", "w")
stitched_list = open("ltr_retro_stitched.list", "w")
for key in matched:
    i = matched[key][0]
    j = matched[key][1]
    ltr_full_name = ltr_full[i]
    int_full_name = int_full[j]
    family_name = ltr_full_name.replace("_LTR#", "#").replace("-LTR#", "#").replace("LTR#", "#")
    ltr_seq = ltr_te_seq[ltr_full_name].seq
    int_seq = ltr_te_seq[int_full_name].seq
    stitched = ">" + str(family_name) + "\n" + str(ltr_seq) + str(int_seq) + str(ltr_seq) + "\n"
    stitched_out.write(stitched)
    stitched_list.write(str(str(ltr_full_name) + "\n"))
    stitched_list.write(str(str(int_full_name) + "\n"))
stitched_out.close()
stitched_list.close()

# Load ltr_retro_stitched fasta and list
stitched_seq = SeqIO.to_dict(SeqIO.parse("ltr_retro_stitched.fasta" ,"fasta"))

with open('ltr_retro_stitched.list', 'r') as file:
    stitched_list = [line.strip() for line in file]

# Remove old entries and sequences
for header in stitched_list:
    te_seq.pop(header)

# Merge the stiched LTR-INT_LTR with the te_seq
merged = stitched_seq | te_seq

# Output the merged fasta file
SeqIO.write(merged.values(), out_file_name, "fasta")
    