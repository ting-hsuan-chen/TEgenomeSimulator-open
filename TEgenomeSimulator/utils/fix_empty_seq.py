import random
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os

# In the case of genome assembly not in the chromosome level, i.e. having multiple
# scaffolds, some short scaffolds may be found as entirely composed of TE sequences.
# These scaffolds are entirely masked as 'Xs' by RepeatMasker and subsequently 
# removed using `sed`, resulting in empty sequence under their headers in the fasta 
# file. This causes error when building index using samtools and the following 
# simulation. To solve this problem, this script checks the fasta file and adds two 
# nucleotides randomly to the empty sequences. 

def add_random_nucleotides(fasta_file, output_file):
    nucleotides = ['A', 'T', 'C', 'G']
    modified_sequences = []

    # Read the FASTA file
    records = SeqIO.parse(fasta_file, "fasta")

    for record in records:
        # Check if the sequence is empty
        if len(record.seq) == 0:
            # Randomly add two nucleotides
            random_seq = random.choice(nucleotides) + random.choice(nucleotides)
            record.seq = Seq(random_seq)  # Wrap the string in a Seq object
            print(f"Modified {record.id} with sequence: {random_seq}")

        modified_sequences.append(record)

    # Write the modified sequences to a new FASTA file
    SeqIO.write(modified_sequences, output_file, "fasta")
    print(f"Modified FASTA saved to {output_file}")

# Input and output file names
input_fasta = sys.argv[1]
name_fasta = os.path.basename(input_fasta)
output_fasta = f"{name_fasta}.emptfixed"

# Check the FASTA file
add_random_nucleotides(input_fasta, output_fasta)

