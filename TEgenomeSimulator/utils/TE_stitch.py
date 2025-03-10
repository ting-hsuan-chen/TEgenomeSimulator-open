from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os
from pathlib import Path
import subprocess

# This script check if a LTR sequence and an INT sequence of a 
# LTR retrotransposon are seperated in the TE library file. If 
# so, this script stitches these two sequences into LTR-INT-LTR 
# format and use the stitched secuqnce to replace the separted 
# LTR & INT seq. LTR or INT seq that don't have a paired mate 
# and other TE sequences are kept unchaged.

# Input and output file names
prefix = sys.argv[1]
input_fasta = sys.argv[2]
outdir = sys.argv[3]
name_fasta = os.path.basename(input_fasta)
final_out = str(outdir) + "/TEgenomeSimulator_" + str(prefix) + "_result"
output_fasta = f"{name_fasta}.stitched"

# Dictionary to store sequences by TE id
te_dict = {}

# Read the input FASTA file
for record in SeqIO.parse(input_fasta, "fasta"):
    header = record.id
    sequence = str(record.seq)
    
    # Extract TE id and type (INT or LTR)
    if "_INT#" in header:
        te_id, te_class = header.split("_INT#")
        te_type = "INT"
    elif "_LTR#" in header:
        te_id, te_class = header.split("_LTR#")
        te_type = "LTR"
    else:
        # Store sequences that don't need stitching directly
        te_dict[header] = record
        continue

    # Store sequences for stitching
    if te_id not in te_dict:
        te_dict[te_id] = {"INT": None, "LTR": None, "class": te_class}
    te_dict[te_id][te_type] = sequence

# List to store the modified sequences
new_records = []

# Process the stored sequences
for te_id, parts in te_dict.items():
    # If it's a dictionary, it means it's a potential pair
    if isinstance(parts, dict):
        # Both INT and LTR parts exist, stitch them
        if parts["INT"] and parts["LTR"]:
            stitched_seq = parts["LTR"] + parts["INT"] + parts["LTR"]
            new_header = f"{te_id}#{parts['class']}"
            
            # Create new SeqRecord
            new_record = SeqRecord(
                Seq(stitched_seq),
                id=new_header,
                description=""
            )
            new_records.append(new_record)
        else:
            # If not a complete pair, keep the original sequences
            if parts["INT"]:
                new_header = f"{te_id}_INT#{parts['class']}"
                new_record = SeqRecord(
                    Seq(parts["INT"]),
                    id=new_header,
                    description=""
                )
                new_records.append(new_record)
            if parts["LTR"]:
                new_header = f"{te_id}_LTR#{parts['class']}"
                new_record = SeqRecord(
                    Seq(parts["LTR"]),
                    id=new_header,
                    description=""
                )
                new_records.append(new_record)
    else:
        # Just a regular sequence, keep as is
        new_records.append(parts)

# Write the modified sequences to the output file
SeqIO.write(new_records, Path(final_out, output_fasta), "fasta")

print(f"Stitched FASTA saved as {output_fasta}")

# Index stitched TE library
input_file = Path(final_out, output_fasta)
log_file = f"{input_file}.faidx.log"

with open(log_file, "w") as log:
    result = subprocess.run(
        ["samtools", "faidx", input_file],
        stderr=log,
        stdout=subprocess.PIPE,
        text=True
        )
