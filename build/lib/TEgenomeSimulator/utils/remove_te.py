def remove_te(input_fasta):
    """
    Remove all occurrences of 'X' from the input FASTA file and save the output.

    Args:
        input_fasta (str): Path to the input FASTA file.

    Returns:
        str: Path to the output file with 'X' removed.
    """
    print(f"\nStarting to remove TE nucleotides using remove_te() function.")

    # Generate the output file name
    output_file = f"../{os.path.basename(input_fasta)}.nonTE"

    # Read the input fasta and write the sequences without 'X'
    with open(output_file, "w") as outfile:
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Remove 'X' from the sequence
            record.seq = record.seq.replace('X', '')
            # Write the updated sequence to the output file
            SeqIO.write(record, outfile, "fasta")

    # Print the full path of the output file
    full_output_path = os.path.abspath(output_file)

    print(f"\nTE nucleotides were removed successfully. Output saved to {full_output_path}.")

# Remove TE nucleotides
masked_file = f"{out_dir}/{base_name}.masked.reformatted"
remove_te(masked_file)