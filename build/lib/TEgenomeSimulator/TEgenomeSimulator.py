import os
import sys
import argparse
import subprocess
from pathlib import Path
from Bio import SeqIO
from datetime import datetime
import subprocess

# Use Tee to duplicate stdout to log file
class TeeLogger:
    def __init__(self, log_file_path):
        self.terminal = sys.stdout
        self.log = open(log_file_path, "a")  # use "a" to append, assuming log already initialized

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

# Use Tee to duplicate stderr to log file
class TeeErrorLogger:
    def __init__(self, log_file_path):
        self.terminal_err = sys.__stderr__
        self.log = open(log_file_path, "a")

    def write(self, message):
        for line in message.rstrip().splitlines():
            timestamped_line = f"[{datetime.now()}] {line}\n"
            self.terminal_err.write(timestamped_line)
            self.log.write(timestamped_line)

    def flush(self):
        self.terminal_err.flush()
        self.log.flush()

# Print subprocess output to both console and log
def run_subprocess_and_tee_output(command, log_path):
    """Run a subprocess and tee its output to both console and log file."""
    with open(log_path, "a") as log_file:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        for line in process.stdout:
            print(line, end='')          # Goes to TeeLogger, so both log and console
            log_file.flush()             # Flush to write immediately

        process.wait()

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, command)

# Function to check if the mode is 1 (user-provided genome) or 0 (random genome)
def mode_check(value):
    if value not in ['0', '1', '2']:
        raise argparse.ArgumentTypeError("Mode must be either '0', '1' or '2'")
    return int(value)

# Function to call TE_stitch.py to combine LTR-INT-LTR for LTR retrotransposons
def run_TE_stitch(prefix, te_lib, final_out):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/TE_stitch.py')
    try:
        # The command to run the external script
        te_stitch_command = [
            'python3', script_path, 
            prefix,
            te_lib,
            final_out
        ]
        
        # Run the command and capture the output
        log_path = os.path.join(final_out, "TEgenomeSimulator.log")
        run_subprocess_and_tee_output(te_stitch_command, log_path)

        print(f"\nStitching TE successfully. Output logged to {final_out}/TEgenomeSimulator.log", flush=True)
    
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred while running TE_stitch.py: {e}", flush=True)
        sys.exit(1)


# Function to call mask_TE.py to mask and remove TE nucleodites
def run_mask_TE(threads, stitched_telib, unmasked_genome, final_out):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/mask_TE.py')
    try:
        # The command to run the external script
        mask_te_command = [
            'python3', script_path, 
            '-t', threads,
            '-r', stitched_telib,
            '-g', unmasked_genome,
            '-o', final_out
        ]
        
        # Run the command and capture the output
        log_path = os.path.join(final_out, "TEgenomeSimulator.log")
        run_subprocess_and_tee_output(mask_te_command, log_path)

        print(f"\nMasking and removing TE successfully. Output logged to {final_out}/TEgenomeSimulator.log", flush=True)
    
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred while running mask_TE.py: {e}", flush=True)
        sys.exit(1)

# Function to call fix_empty_seq.py to mask and remove TE nucleodites
def run_fix_empty_seq(non_te_genome, final_out):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/fix_empty_seq.py')
    try:
        # The command to run the external script
        fix_empty_seq_command = [
            'python3', script_path, 
            non_te_genome
        ]
        
        # Run the command and capture the output
        log_path = os.path.join(final_out, "TEgenomeSimulator.log")
        run_subprocess_and_tee_output(fix_empty_seq_command, log_path)

        print(f"\nChecking empty sequences successfully. Output logged to {final_out}/TEgenomeSimulator.log", flush=True)
    
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred while running fix_empty_seq.py: {e}", flush=True)
        sys.exit(1)


# Function to run the three steps required for "--to_mask" option: run_TE_stitch(), run_mask_TE(), and run_fix_empty_seq()
def to_mask(args, final_out):
    """
    Runs the TE masking pipeline, including TE stitching, masking the genome, 
    and fixing empty sequences.

    Args:
        args (Namespace): Parsed command-line arguments with required parameters.
        final_out (str): The output directory for final processed files.

    Returns:
        Path: Path to the fixed non-TE genome file.
    """
    # Specify the input and output of run_TE_stitch()
    stitched_te_fasta_name = f"{os.path.basename(args.repeat2)}.stitched"
    stitched_telib = Path(final_out, stitched_te_fasta_name)

    # Specify the input genome for run_mask_TE()
    unmasked_genome = str(args.genome)

    # Specify the input genome for run_fix_empty_seq()
    non_te_genome_name = f"{os.path.basename(unmasked_genome)}.masked.reformatted.nonTE"
    non_te_genome = Path(final_out, non_te_genome_name)

    # Specify the input genome for run_prep_config_custom()
    fixed_non_te_genome = Path(final_out, f"{os.path.basename(unmasked_genome)}.masked.reformatted.nonTE.emptfixed")

    # Run the required preparation steps
    run_TE_stitch(args.prefix, args.repeat2, final_out)
    run_mask_TE(str(args.threads), stitched_telib, unmasked_genome, final_out)
    run_fix_empty_seq(non_te_genome, final_out)

    # Return the path to the final processed genome file
    return fixed_non_te_genome


# Function to call summarise_rm_out.py to summarise TE composition
def run_summarise_rm_out(rmasker_tbl, libindex, rmasker_out, final_out):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/summarise_rm_out.py')
    try:
        # The command to run the external script
        summarise_rm_out_command = [
            'python3', script_path, 
            rmasker_tbl,
            libindex,
            rmasker_out,
            final_out
        ]
        
        # Run the command and capture the output
        log_path = os.path.join(final_out, "TEgenomeSimulator.log")
        run_subprocess_and_tee_output(summarise_rm_out_command, log_path)

        print(f"\nSummarising TE composition successfully. Output logged to {final_out}/TEgenomeSimulator.log", flush=True)
    
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred while running summarise_rm_out.py: {e}", flush=True)
        sys.exit(1)


# Function to call the prep_sim_TE_lib.py script to generate the TE library table
def run_prep_sim_TE_lib(prefix, repeat, maxcp, mincp, maxidn, minidn, maxsd, minsd, intact, seed, final_out):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/prep_sim_TE_lib.py')
    try:
        # The command to run the external script
        prep_telib_command = [
            'python3', script_path, 
            '-p', prefix,
            '-r', repeat, 
            '-m', str(maxcp),
            '-n', str(mincp),
            '--maxidn', str(maxidn),
            '--minidn', str(minidn),
            '--maxsd', str(maxsd),
            '--minsd', str(minsd),
            '-i', str(intact),
            '-s', str(seed),
            '-o', final_out
        ]
        
        # Run the command and capture the output
        log_path = os.path.join(final_out, "TEgenomeSimulator.log")
        run_subprocess_and_tee_output(prep_telib_command, log_path)

        print(f"\nTE library table generated successfully. Output logged to {final_out}/TEgenomeSimulator.log", flush=True)
    
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred while running prep_sim_TE_lib.py: {e}", flush=True)
        sys.exit(1)

def prep_sim_TE_lib_mode2(prefix, rmout_summary, seed, final_out):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/prep_sim_TE_lib_mode2.py')
    try:
        # The command to run the external script
        prep_telib_m2_command = [
            'python3', script_path, 
            prefix,
            rmout_summary, 
            str(seed), 
            final_out
        ]
        
        # Run the command and capture the output
        log_path = os.path.join(final_out, "TEgenomeSimulator.log")
        run_subprocess_and_tee_output(prep_telib_m2_command, log_path)

        print(f"\nTE library table for mode 2 generated successfully. Output logged to {final_out}/TEgenomeSimulator.log", flush=True)
    
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred while running prep_sim_TE_lib_mode2.py: {e}", flush=True)
        sys.exit(1)

# Function to call the prep_yml_config.py script for Random Genome Mode
def run_prep_config_random(prefix, chridx, repeat, te_table, seed, final_out):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/prep_yml_config.py')
    try:
        # Construct the command to run the prep_yml_config.py script using chromosome index file
        prep_yml_command = [
            'python3', script_path, 
            '-p', prefix, 
            '-c', chridx, 
            '-r', repeat, 
            '-t', te_table, 
            '-s', str(seed), 
            '-o', final_out
        ]

        # Run the command and append the output to the log file
        log_path = os.path.join(final_out, "TEgenomeSimulator.log")
        run_subprocess_and_tee_output(prep_yml_command, log_path)

        print(f"\nConfig file generated successfully. Output logged to {final_out}/TEgenomeSimulator.log", flush=True)
        
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred while running prep_yml_config.py: {e}", flush=True)
        sys.exit(1)

# Function to call the prep_yml_config.py script for Custom Genome Mode
def run_prep_config_custom(prefix, genome, repeat, te_table, seed, final_out):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/prep_yml_config.py')
    try:
        # Construct the command to run the prep_yml_config.py script using genome fasta file
        prep_yml_command = [
            'python3', script_path, 
            '-p', prefix, 
            '-g', genome, 
            '-r', repeat, 
            '-t', te_table, 
            '-s', str(seed), 
            '-o', final_out
        ]
            
        # Run the command and append the output to the log file
        log_path = os.path.join(final_out, "TEgenomeSimulator.log")
        run_subprocess_and_tee_output(prep_yml_command, log_path)    
        
        print(f"\nConfig file generated successfully. Output logged to {final_out}/TEgenomeSimulator.log", flush=True)
        
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred while running prep_yml_config.py: {e}", flush=True)
        sys.exit(1)

# Function to call TE_sim_random_insertion.py for non-overlape random TE insertion
def run_TE_sim_random_insertion(mode, prefix, frag_mode, alpha, beta, final_out):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/TE_sim_random_insertion.py')
    try:
        # Construct the command to run the TE_sim_random_insertion.py script
        prep_sim_command = [
            'python3', script_path, 
            '-M', str(mode),
            '-p', prefix, 
            '-f', str(frag_mode),
            '-a', str(alpha),
            '-b', str(beta), 
            '-o', final_out
        ]

        # Run the command and capture the output
        log_path = os.path.join(final_out, "TEgenomeSimulator.log")
        run_subprocess_and_tee_output(prep_sim_command, log_path)

        print(f"\nGenome with non-overlap random TE insertions was generated successfully. Output logged to {final_out}/TEgenomeSimulator.log", flush=True)
    
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred while running prep_sim_TE_lib.py: {e}", flush=True)
        sys.exit(1)

# Function to call for nested TE insertion
def run_TE_sim_nested_insertion(mode, prefix, frag_mode, alpha, beta, final_out):
    script_path = os.path.join(os.path.dirname(__file__), 'utils/TE_sim_nested_insertion.py')
    try:
        # Construct the command to run the TE_sim_nested_insertion.py script
        prep_nest_command = [
            'python3', script_path, 
            '-M', str(mode),
            '-p', prefix, 
            '-f', str(frag_mode),
            '-a', str(alpha),
            '-b', str(beta),
            '-o', final_out
        ]

        # Run the command and capture the output
        log_path = os.path.join(final_out, "TEgenomeSimulator.log")
        run_subprocess_and_tee_output(prep_nest_command, log_path)

        print(f"\nGenome with non-overlap random and nested TE insertions was generated successfully. Output logged to {final_out}/TEgenomeSimulator.log", flush=True)
    
    except subprocess.CalledProcessError as e:
        print(f"\nError occurred while running prep_sim_TE_lib.py: {e}", flush=True)
        sys.exit(1)


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="main arguments of TEgenomeSimulator to simulate TE mutation and insertion into genome.")

    # Define arguments
    parser.add_argument('-M', '--mode', type=mode_check, choices=[0, 1, 2], help="Mode for genome simulation (choose from 0, 1 or 2).", required=True)
    parser.add_argument('-k', '--to_mask', action='store_true', help="Mask and remove TE from provided genome when enabled")
    parser.add_argument('-p', '--prefix', type=str, help="Prefix for output files.", required=True)
    parser.add_argument('-r', '--repeat', type=str, help="TE family fasta file for random insertion simulation.", required=True)
    parser.add_argument('-r2', '--repeat2', type=str, help="TE family fasta file for masking genome (required when using --to_mask).")
    parser.add_argument('-m', '--maxcp', type=int, default=10, help="Maximum copies of TE family (default is 10).")
    parser.add_argument('-n', '--mincp', type=int, default=1, help="Minimum copies of TE family (default is 1).")
    parser.add_argument('--maxidn', type=int, default=95, help="The upper bound of mean sequence identity to be sampled for each TE family (default is 95; i.e. 95 percent).")
    parser.add_argument('--minidn', type=int, default=80, help="The lower bound of mean sequence identity to be sampled for each TE family (default is 80; i.e. 80 percent).")
    parser.add_argument('--maxsd', type=int, default=20, help="The upper bound of standard deviation of mean identity to be sampled for each TE family (default is 20).")
    parser.add_argument('--minsd', type=int, default=1, help="The lower bound of standard deviation of mean identity to be sampled for each TE family (default is 1).")
    parser.add_argument('-c', '--chridx', type=str, help="Chromosome index file if mode 0 is selected.")
    parser.add_argument('-g', '--genome', type=str, help="Genome fasta file if mode 1 or 2 is selected.")
    parser.add_argument('-f', '--frag_mode', type=int, default=0, help="Mode for TE fragmentation (0, 1, or 2). If frag_mode is 2, uses an integrity list read from --te_table, otherwise beta distribution.")
    parser.add_argument('-a', '--alpha', type=float, default=0.5, help="Alpha value for the beta distribution used for fragmentation simulation (default is 0.5).")
    parser.add_argument('-b', '--beta', type=float, default=0.7, help="Beta value for the beta distribution used for fragmentation simulation (default is 0.7).")
    parser.add_argument('-i', '--intact', type=float, default=0.001, help="Maximum probability of inserting intact TEs per family (default is 0.001; i.e. 0.1 percent).")
    parser.add_argument('-s', '--seed', type=int, default=1, help="Random seed (default is 1).")
    parser.add_argument('-t', '--threads', type=int, default=1, help="Threads for running RepeatMasker (default is 1)")
    parser.add_argument('-o', '--outdir', type=str, help="Output directory.", required=True)
    parser.add_argument('--te_table', type=str, help="Path to an existing TE mutagenesis settings table (modes 0 and 1 only).")
    parser.add_argument('--te_scan_only', action='store_true', help="Mode 2 only: stop after TE composition profiling and TE table generation.")
    
    # Parse arguments
    args = parser.parse_args()

    # Validate --te_table and --te_scan_only options
    if args.te_table:
        if args.mode not in [0, 1]:
            print("Error: --te_table can only be used with mode 0 or mode 1.", flush=True)
            sys.exit(1)
        args.te_table = Path(args.te_table).resolve()
        if not args.te_table.is_file():
            print(f"Error: The specified TE table file does not exist: {args.te_table}", flush=True)
            sys.exit(1)

    if args.te_scan_only and args.mode != 2:
        print("Error: --te_scan_only can only be used with mode 2.", flush=True)
        sys.exit(1)

    # Convert user-provided output path to absolute path
    outdir = Path(args.outdir).resolve()
    
    # Specify final output dir for each project
    final_out = outdir / f"TEgenomeSimulator_{args.prefix}_result"
    Path(final_out).mkdir(parents=True, exist_ok=True)
    
    # Set up path to log file
    log_path = os.path.join(final_out, "TEgenomeSimulator.log")

    # Create and initialize the log file early
    with open(log_path, "w") as log_file:
        print(f"[{datetime.now()}] TEgenomeSimulator started.")
        print(f"[{datetime.now()}] Arguments: {vars(args)}")

    # Redirect stdout and stderr to TeeLogger after initialization
    sys.stdout = TeeLogger(log_path)
    sys.stderr = TeeErrorLogger(log_path)

    print(f"Output Directory: {final_out}", flush=True)

    # Mode-based validation
    if args.mode == 0:
        # Mode 0 requires --chridx to be specified
        if not args.chridx:
            print("Error: When --mode 0 is selected, --chridx must also be specified.", flush=True)
            sys.exit(1)
    elif args.mode == 1:
        # Mode 1 requires --genome to be specified
        if not args.genome:
            print("Error: When --mode 1 is selected, --genome must also be specified.", flush=True)
            sys.exit(1)
        # If --to_mask is specified, then --repeat2 must also be provided
        if args.to_mask and not args.repeat2:
            print("Error: When --mode 1 is selected and --to_mask is specified, --repeat2 must also be provided.", flush=True)
            sys.exit(1)
    elif args.mode == 2:
        # Mode 2 requires --genome and --repeat2 to be specified
        if not args.genome or not args.repeat2:
            print("Error: When --mode 2 is selected, --genome and --repeat2 must also be specified.", flush=True)
            sys.exit(1)


    # Output parsed arguments (for demonstration)
    print(f"Mode: {args.mode}", flush=True)
    if args.mode == 0:
        print(f"Running Random Synthesized mode.", flush=True)
    elif args.mode == 1:
        print(f"Running Custom Genome mode.", flush=True)
        if args.to_mask:
            print(f"To mask genome? {args.to_mask} enabled. Yes.", flush=True)
        else:
            print(f"To mask genome? No.")
    elif args.mode == 2:
        print(f"Running TE Composition Approximation mode. TEgenomeSimulator will mask the genome under this mode.", flush=True)
        
    print(f"Prefix: {args.prefix}", flush=True)
    print(f"Repeat: {args.repeat}", flush=True)
    if args.mode == 1 or args.mode == 2:
        print(f"Repeat2: {args.repeat2}", flush=True)
    print(f"Chromosome Index: {args.chridx}", flush=True)
    print(f"Genome File: {args.genome}", flush=True)
    

    if args.mode == 0 or args.mode == 1:
        if args.te_table:
            print(f"Using user-provided TE mutagenesis settings table: {args.te_table}", flush=True)
            if args.frag_mode in [0, 1]:
               print(f"frag_mode={args.frag_mode}. Expecting 11 columns in the table and beta distribution will be used for integrity simulation. See README for table format.", flush=True)
            if args.frag_mode == 2:
               print(f"frag_mode={args.frag_mode}. Expecting 12 columns in the table and will use the integrity list for integrity simulation. See README for table format.", flush=True)
        else:  
            print(f"Max Copies: {args.maxcp}", flush=True)
            print(f"Min Copies: {args.mincp}", flush=True)        
            print(f"Upper bound of mean identity: {args.maxidn}", flush=True)
            print(f"Lower bound of mean identity: {args.minidn}", flush=True)
            print(f"Upper bound of sd of mean identity: {args.maxsd}", flush=True)
            print(f"Lower bound of sd of mean ideneity: {args.minsd}", flush=True)
            print(f"Max chance of intact insertion: {args.intact}", flush=True)
            print(f"Alpha: {args.alpha}", flush=True)
            print(f"Beta: {args.beta}", flush=True)
    else:
        print(f"TE copy number, mean sequence identity, sd, and max chance of intact insertion evaluated from repeatmasker output.", flush=True)

    print(f"Seed: {args.seed}", flush=True)

    # ===================== Mode-specific steps =====================

    # Call the prep_sim_TE_lib.py script to generate the TE library table
    if args.mode in [0, 1]:
        if args.te_table:
            te_table = args.te_table            
        else:
            print("Generating TE mutagenesis settings table...", flush=True)
            run_prep_sim_TE_lib(args.prefix, args.repeat, args.maxcp, args.mincp, args.maxidn, args.minidn, args.maxsd, args.minsd, args.intact, args.seed, final_out)
            te_table = os.path.join(final_out, "TElib_sim_list.table")

    # Call the prep_yml_config.py script to generate the config file using the TE library table generated from previous step or supplied by users    
    if args.mode == 0:
        print("mode=0, running prep_yml_config.py for Random Genome Mode.", flush=True)
        run_prep_config_random(args.prefix, args.chridx, args.repeat, te_table, args.seed, final_out)

    elif args.mode == 1:
        if not args.to_mask: 
            print("mode=1, and --to_mask not specified, assumming user-provided genome is TE-depleted.", 
                  "Running prep_yml_config.py for Custom Genome Mode.", flush=True) 
            run_prep_config_custom(args.prefix, args.genome, args.repeat, te_table, args.seed, final_out)
        else: 
            # when --to_mask is enabled
            print("mode=1, and --to_mask is specified, assumming user-provided genome is not TE-depleted.", 
                  "Running mask_TE.py before prep_yml_config.py for Custom Genome Mode.", flush=True)            
            fixed_non_te_genome = to_mask(args, final_out)
            run_prep_config_custom(args.prefix, fixed_non_te_genome, args.repeat, te_table, args.seed, final_out)

    elif args.mode == 2:
        print("mode=2, summarising TE composition and masking TE in user-provided genome.", 
              "Running mask_TE.py and summarise_rm_out.py before prep_yml_config.py for TE Composition Approximation Mode.", flush=True)
        # --to_mask is enabled under this mode
        fixed_non_te_genome = to_mask(args, final_out)
        # Specify the required input files for run_summarise_rm_out()
        repeat2_name = f"{os.path.basename(args.repeat2)}"
        unmasked_genome_name = f"{os.path.basename(args.genome)}"
        libindex = os.path.join(final_out, f"{repeat2_name}.stitched.fai")
        rmasker_tbl = os.path.join(final_out, "repeatmasker", f"{unmasked_genome_name}.tbl")
        rmasker_out = os.path.join(final_out, "repeatmasker", f"{unmasked_genome_name}.out")
        # Create the summary files
        run_summarise_rm_out(rmasker_tbl, libindex, rmasker_out, final_out)
        rmout_summary = os.path.join(final_out, "summarise_repeatmasker_out_family.txt")
        prep_sim_TE_lib_mode2(args.prefix, rmout_summary, args.seed, final_out)
        # Use the summarise file as te_table
        te_table = os.path.join(final_out, "TElib_sim_list_mode2.table")
        # Check whether to proceed to full mode 2 simulation
        if args.te_scan_only:
            print(f"Mode 2: te_scan_only was enabled.")
            print(f"TE composition profiling complete. TE mutagenesis settings table generated at {te_table}. Stopping as requested.", flush=True)
            sys.exit(0)
        # Prepare config .yaml file
        sim_repeat = os.path.join(final_out, f"{repeat2_name}.stitched")
        run_prep_config_custom(args.prefix, fixed_non_te_genome, sim_repeat, te_table, args.seed, final_out)
   
    # ===================== Final simulation steps =====================
    # Non-overlape random TE insertion
    run_TE_sim_random_insertion(args.mode, args.prefix, args.frag_mode, args.alpha, args.beta, final_out)

    # Nested TE insertion
    run_TE_sim_nested_insertion(args.mode, args.prefix, args.frag_mode, args.alpha, args.beta, final_out)

if __name__ == "__main__":
    main()


