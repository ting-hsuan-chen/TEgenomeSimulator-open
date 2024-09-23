#!/bin/bash

ml conda
conda activate TEgenomeSimulator

# 2024-04-17 rerun simulation after the mis-calculation of sd has been fixed
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=1G
TIME="06:00:00"
JOB="TEgenSim_tair10"
sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ./log/${JOB}.out
#SBATCH -e ./log/${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

python $TEgenomeSimulator/custom_genome_seq_TEs_tair10.py
python $TEgenomeSimulator/custom_genome_nest_TEs_tair10.py

EOF
#Submitted batch job 5256967
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:03:41
CPU Efficiency: 91.70% of 00:04:01 core-walltime
Job Wall-clock time: 00:04:01
Memory Utilized: 444.69 MB
Memory Efficiency: 43.43% of 1.00 GB


# 2024-03-27
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=1G
TIME="06:00:00"
JOB="TEgenSim_tair10"
sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ./log/${JOB}.out
#SBATCH -e ./log/${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

python $TEgenomeSimulator/custom_genome_seq_TEs_tair10.py
python $TEgenomeSimulator/custom_genome_nest_TEs_tair10.py

EOF
#Submitted batch job 5019564
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:04:44
#CPU Efficiency: 94.35% of 00:05:01 core-walltime
#Job Wall-clock time: 00:05:01
#Memory Utilized: 426.00 MB
#Memory Efficiency: 41.60% of 1.00 GB
