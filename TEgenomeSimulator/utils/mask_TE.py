import subprocess
import argparse
import os
import sys
from Bio import SeqIO

def setup_repeatmasker_input(genome, repeat):
    """
    Set up the input files for RepeatMasker.

    Args:
        out (str): Working and output directory.
        genome (str): Genome name.
        ref (str): Path to the reference FASTA file.

    Returns:
        tuple: (namefa, namelib)
    """

    # Create symbolic links
    if os.path.exists(os.path.basename(genome)) or os.path.islink(os.path.basename(genome)):  
        os.remove(os.path.basename(genome))
    if os.path.exists(os.path.basename(repeat)) or os.path.islink(os.path.basename(repeat)):  
        os.remove(os.path.basename(repeat))
    os.symlink(genome, os.path.basename(genome))
    os.symlink(repeat, os.path.basename(repeat))

    # Extract file names
    namefa = os.path.basename(genome)
    namelib = os.path.basename(repeat)

    return namefa, namelib

def run_repeatmasker(threads, namelib, namefa, sif_path, outdir):
    """
    Run RepeatMasker with the given parameters.

    Args:
        threads (int): Number of threads to use.
        namelib (str): Path to the RepeatMasker library file.
        namefa (str): Path to the input FASTA file.
        sif_path (str): Path to the dfam-tetools.sif container.

    Returns:
        tuple: (stdout, stderr, returncode)
    """
    try: 
        # The command to run repeatmasker
        cmd = [
            "singularity", "exec", sif_path, "RepeatMasker",
            "-pa", str(threads),
            "-a", "-x", "-q", "-no_is", "-norna", "-nolow",
            "-div", "40",
            "-lib", namelib,
            "-cutoff", "225",
            namefa
            ]
    
            # Run the command and capture the output
        with open(f"{outdir}/TEgenomeSimulator.log", "a") as log_file:
            subprocess.run(cmd, check=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        print(f"\nRunning RepeatMasker successfully. Output logged to {outdir}/TEgenomeSimulator.log")
    
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred while running fix_empty_seq.py: {e}")
        sys.exit(1)



def convert_fasta_to_single_line(input_fasta):
    """
    Convert a multi-line FASTA file into a single-line FASTA file.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to save the formatted FASTA file.
    """
    base_name = os.path.basename(input_fasta)
    output_fasta = f"{base_name}.reformatted"

    with open(output_fasta, "w") as out_file:
        for record in SeqIO.parse(input_fasta, "fasta"):
            out_file.write(f">{record.id}\n{str(record.seq)}\n")

    print(f"\nRunning convert_fasta_to_single_line() successfully.")


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
    


def main():
    # parse argument
    parser = argparse.ArgumentParser()
    parser.add_argument("-t","--threads", type=str,
                    help="path to the repeat fasta file")
    parser.add_argument("-r","--repeat", type=str,
                    help="path to the repeat fasta file")
    parser.add_argument("-g", "--genome", type=str,
                    help="path to genome fasta file")
    parser.add_argument("-s", "--sif_path", type=str,
                    help="path to dfam-tetools.sif")
    parser.add_argument("-o", "--outdir", type=str,
                    help="output directory")
    
    # Parse arguments
    args = parser.parse_args()
    threads = args.threads
    te_lib = args.repeat
    genome_fa = args.genome
    sif_path = args.sif_path
    out = args.outdir
    
    # Change directory
    out_dir = os.path.join(out, "repeatmasker")
    os.makedirs(out_dir, exist_ok=True)  # Ensure directory exists
    os.chdir(out_dir)

    # Setup soft-link
    namefa, namelib = setup_repeatmasker_input(genome_fa, te_lib)

    # Run RepeatMasker
    run_repeatmasker(threads, namelib, namefa, sif_path, out)

    # Convert multi-line masked fasta file to single-line fasta
    base_name = os.path.basename(genome_fa)
    input_fasta = f"{base_name}.masked"
    convert_fasta_to_single_line(input_fasta)
    
    # Remove TE nucleotides
    masked_file = f"{out_dir}/{base_name}.masked.reformatted"
    remove_te(masked_file)

if __name__ == "__main__":
    main()
