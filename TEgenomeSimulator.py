import os
import sys
import argparse
import subprocess
from pathlib import Path
#import random
#import re
#from Bio import SeqIO

# function to check if the mode is 1 (user-provided genome) or 0 (random genome)
def mode_check(value):
    if value not in ['0', '1']:
        raise argparse.ArgumentTypeError("Mode must be either '0' or '1'")
    return int(value)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="main arguments of TEgenomeSimulator to simulate TE mutation and insertion into genome.")

    # Define arguments
    parser.add_argument('-M', '--mode', type=mode_check, choices=[0, 1], help="Mode for genome simulation (either 0 or 1).", required=True)
    parser.add_argument('-p', '--prefix', type=str, help="Prefix for output files.", required=True)
    parser.add_argument('-r', '--repeat', type=str, help="TE family fasta file.", required=True)
    parser.add_argument('-m', '--maxcp', type=int, help="Maximum copies of TE family.", required=True)
    parser.add_argument('-n', '--mincp', type=int, help="Minimum copies of TE family.", required=True)
    parser.add_argument('-c', '--chridx', type=str, help="Chromosome index file if mode 0 is selected.")
    parser.add_argument('-g', '--genome', type=str, help="Genome fasta file if mode 1 is selected.")
    parser.add_argument('-s', '--seed', type=int, default=1, help="Random seed (default is 1).")
    parser.add_argument('-o', '--outdir', type=str, help="Output directory.", required=True)
    
    # Parse arguments
    args = parser.parse_args()

    # Mode-based validation
    if args.mode == 0:
        # Mode 0 requires --chridx to be specified
        if not args.chridx:
            print("Error: When --mode 0 is selected, --chridx must also be specified.")
            sys.exit(1)
    elif args.mode == 1:
        # Mode 1 requires --genome to be specified
        if not args.genome:
            print("Error: When --mode 1 is selected, --genome must also be specified.")
            sys.exit(1)

    # Output parsed arguments (for demonstration)
    print(f"Mode: {args.mode}")
    print(f"Prefix: {args.prefix}")
    print(f"Repeat: {args.repeat}")
    print(f"Max Copies: {args.maxcp}")
    print(f"Min Copies: {args.mincp}")
    print(f"Chromosome Index: {args.chridx}")
    print(f"Genome File: {args.genome}")
    print(f"Seed: {args.seed}")
    print(f"Output Directory: {args.outdir}")
    
    # Specify output dir for each project
    final_out = args.outdir + "/TEgenomeSimulator_" + args.prefix + "_result"
    Path(final_out).mkdir(parents=True, exist_ok=True)

    # Call the prep_sim_TE_lib.py script to generate the TE library table
    try:
        # The command to run the external script
        prep_telib_command = [
            'python3', './scripts/prep_sim_TE_lib.py', 
            '-p', args.prefix,
            '-r', args.repeat, 
            '-m', str(args.maxcp),
            '-n', str(args.mincp),
            '-o', args.outdir
        ]
        
        # Run the command and capture the output
        
        with open(f"{final_out}/TEgenomeSimulator.log", "w") as log_file:
            subprocess.run(prep_telib_command, check=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        print(f"TE library table generated successfully. Output logged to {final_out}/TEgenomeSimulator.log")
    
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running prep_sim_TE_lib.py: {e}")
        sys.exit(1)

    # Call the prep_yml_config.py script to generate the config file using the TE library table generated from previous step
    te_table = os.path.join(final_out, "TElib_sim_list.table")

    if args.mode == 0:
        print("mode=0, running prep_yml_config.py for Random Genome Mode.")
        try:
            # Construct the command to run the prep_yml_config.py script using chromosome index file
            prep_yml_command = [
                'python3', './scripts/prep_yml_config.py', 
                '-p', args.prefix, 
                '-c', args.chridx, 
                '-r', args.repeat, 
                '-t', te_table, 
                '-s', str(args.seed), 
                '-o', args.outdir
            ]

            # Run the command and append the output to the log file
            with open(f"{final_out}/TEgenomeSimulator.log", "a") as log_file:
                subprocess.run(prep_yml_command, check=True, stdout=log_file, stderr=subprocess.STDOUT)
            
            print(f"Config file generated successfully. Output logged to {final_out}/TEgenomeSimulator.log")
        
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running prep_yml_config.py: {e}")
            sys.exit(1)
    elif args.mode == 1:
        print("mode=1, running prep_yml_config.py for Custom Genome Mode.")
        try:
            # Construct the command to run the prep_yml_config.py script using genome fasta file
            prep_yml_command = [
                'python3', './scripts/prep_yml_config.py', 
                '-p', args.prefix, 
                '-g', args.genome, 
                '-r', args.repeat, 
                '-t', te_table, 
                '-s', str(args.seed), 
                '-o', args.outdir
            ]
            
            # Run the command and append the output to the log file
            with open(f"{final_out}/TEgenomeSimulator.log", "a") as log_file:
                subprocess.run(prep_yml_command, check=True, stdout=log_file, stderr=subprocess.STDOUT)
            
            print(f"Config file generated successfully. Output logged to {final_out}/TEgenomeSimulator.log")
        
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running prep_yml_config.py: {e}")
            sys.exit(1)
    
    # Non-overlape random TE insertion
    try:
        # Construct the command to run the TE_sim_random_insertion.py script
        prep_sim_command = [
            'python3', './scripts/TE_sim_random_insertion.py', 
            '-M', str(args.mode),
            '-p', args.prefix, 
            '-o', args.outdir
        ]

        # Run the command and capture the output
        with open(f"{final_out}/TEgenomeSimulator.log", "w") as log_file:
            subprocess.run(prep_sim_command, check=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        print(f"Genome with non-overlap random TE insertions was generated successfully. Output logged to {final_out}/TEgenomeSimulator.log")
    
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running prep_sim_TE_lib.py: {e}")
        sys.exit(1)

    # Nested TE insertion
    try:
        # Construct the command to run the TE_sim_nested_insertion.py script
        prep_nest_command = [
            'python3', './scripts/TE_sim_nested_insertion.py', 
            '-M', str(args.mode),
            '-p', args.prefix, 
            '-o', args.outdir
        ]

        # Run the command and capture the output
        with open(f"{final_out}/TEgenomeSimulator.log", "w") as log_file:
            subprocess.run(prep_nest_command, check=True, stdout=log_file, stderr=subprocess.STDOUT)
        
        print(f"Genome with non-overlap random and nested TE insertions was generated successfully. Output logged to {final_out}/TEgenomeSimulator.log")
    
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running prep_sim_TE_lib.py: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()


