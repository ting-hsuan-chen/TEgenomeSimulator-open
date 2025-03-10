import subprocess
import argparse
import os
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
    os.symlink(genome, os.path.basename(genome))
    os.symlink(repeat, os.path.basename(repeat))

    # Extract file names
    namefa = os.path.basename(genome)
    namelib = os.path.basename(repeat)

    return namefa, namelib

def run_repeatmasker(threads, namelib, namefa):
    """
    Run RepeatMasker with the given parameters.

    Args:
        threads (int): Number of threads to use.
        namelib (str): Path to the RepeatMasker library file.
        namefa (str): Path to the input FASTA file.

    Returns:
        tuple: (stdout, stderr, returncode)
    """

    cmd = [
        "RepeatMasker",
        "-pa", str(threads),
        "-a", "-x", "-q", "-no_is", "-norna", "-nolow",
        "-div", "40",
        "-lib", namelib,
        "-cutoff", "225",
        namefa
    ]

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Define file paths
    #base_name = os.path.basename(namefa)
    stdout_file = "repeatmasker.stdout.txt"
    stderr_file = "repeatmasker.stderr.txt"
    returncode_file = "repeatmasker.returncode.txt"

    # Save outputs to files
    with open(stdout_file, "w") as f:
        f.write(result.stdout)

    with open(stderr_file, "w") as f:
        f.write(result.stderr)

    with open(returncode_file, "w") as f:
        f.write(str(result.returncode))

    return result.stdout, result.stderr, result.returncod


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

def remove_te(masked_file):
    """
    Remove all occurrences of 'X' from the masked file and save the output.

    Args:
        masked_file (str): Path to the input masked file.

    Returns:
        str: Path to the output file with 'X' removed.
    """
    output_file = f"../{os.path.basename(masked_file)}.nonTE"

    with open(masked_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            outfile.write(line.replace("X", ""))


def main():
    # parse argument
    parser = argparse.ArgumentParser()
    parser.add_argument("-t","--threads", type=str,
                    help="path to the repeat fasta file")
    parser.add_argument("-r","--repeat", type=str,
                    help="path to the repeat fasta file")
    parser.add_argument("-g", "--genome", type=int,
                    help="path to genome fasta file")
    parser.add_argument("-o", "--outdir", type=str,
                    help="output directory")
    
    # Parse arguments
    args = parser.parse_args()
    threads = args.threads
    te_lib = args.repeat
    genome_fa = args.genome
    out = args.outdir
    
    # Change directory
    out_dir = os.path.join(out, "repeatmasker")
    os.makedirs(out_dir, exist_ok=True)  # Ensure directory exists
    os.chdir(out_dir)

    # Setup soft-link
    namefa, namelib = setup_repeatmasker_input(genome_fa, te_lib)

    # Run RepeatMasker
    run_repeatmasker(threads, namelib, namefa)

    # Convert multi-line masked fasta file to single-line fasta
    base_name = os.path.basename(genome_fa)
    input_fasta = f"{base_name}.masked"
    convert_fasta_to_single_line(input_fasta)
    
    # Remove TE nucleotides
    masked_file = f"{base_name}.masked.reformatted"
    remove_te(masked_file)

if __name__ == "__main__":
    main()
