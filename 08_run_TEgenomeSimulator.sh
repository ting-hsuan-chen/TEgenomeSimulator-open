#!/bin/bash

ml conda
conda activate TEgenomeSimulator
# 21 Feb 2024 after debugging "--" in the simulated genome
# Try TE loci: 5-10 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=1G
TIME="06:00:00"
JOB="TEgenSim_5_10"
sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ./log/${JOB}.out
#SBATCH -e ./log/${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_10.py
python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_10.py

EOF
Submitted batch job 4438066
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:00:27
CPU Efficiency: 52.94% of 00:00:51 core-walltime
Job Wall-clock time: 00:00:51
Memory Utilized: 842.30 MB
Memory Efficiency: 82.26% of 1.00 GB

# Try TE loci: 5-100 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=1G
TIME="06:00:00"
JOB="TEgenSim_5_100"
sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ./log/${JOB}.out
#SBATCH -e ./log/${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_100.py
python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_100.py

EOF
Submitted batch job 4438072
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:03:13
CPU Efficiency: 88.53% of 00:03:38 core-walltime
Job Wall-clock time: 00:03:38
Memory Utilized: 879.99 MB
Memory Efficiency: 85.94% of 1.00 GB

# Try TE loci: 5-500 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=5G
TIME="3-00:00:00"
JOB="TEgenSim_5_500"
sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ./log/${JOB}.out
#SBATCH -e ./log/${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_500.py
python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_500.py

EOF
Submitted batch job 4438074
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:54:48
CPU Efficiency: 96.99% of 00:56:30 core-walltime
Job Wall-clock time: 00:56:30
Memory Utilized: 2.81 GB
Memory Efficiency: 56.13% of 5.00 GB


# Try TE loci: 5-1000 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=10G
TIME="6-00:00:00"
JOB="TEgenSim_5_1000"
sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ./log/${JOB}.out
#SBATCH -e ./log/${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_1000.py
python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_1000.py

EOF
Submitted batch job 4438075
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 04:38:20
CPU Efficiency: 99.77% of 04:38:58 core-walltime
Job Wall-clock time: 04:38:58
Memory Utilized: 5.10 GB
Memory Efficiency: 50.98% of 10.00 GB

















# Feb 2024 after debugging the coordinate problem
# Try TE loci: 5-10 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=1G
TIME="06:00:00"
JOB="TEgenSim_5_10"
sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ./log/${JOB}.out
#SBATCH -e ./log/${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_10.py
python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_10.py

EOF
Submitted batch job 4387771
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:00:19
CPU Efficiency: 82.61% of 00:00:23 core-walltime
Job Wall-clock time: 00:00:23
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 1.00 GB (1.00 GB/node)

# Try TE loci: 5-100 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=1G
TIME="06:00:00"
JOB="TEgenSim_5_100"
sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ./log/${JOB}.out
#SBATCH -e ./log/${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_100.py
python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_100.py

EOF
Submitted batch job 4387932
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:02:31
CPU Efficiency: 94.97% of 00:02:39 core-walltime
Job Wall-clock time: 00:02:39
Memory Utilized: 1016.80 MB
Memory Efficiency: 99.30% of 1.00 GB


# Try TE loci: 5-500 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=5G
TIME="3-00:00:00"
JOB="TEgenSim_5_500"
sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ./log/${JOB}.out
#SBATCH -e ./log/${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_500.py
python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_500.py

EOF
Submitted batch job 4387975
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:57:20
CPU Efficiency: 98.85% of 00:58:00 core-walltime
Job Wall-clock time: 00:58:00
Memory Utilized: 2.88 GB
Memory Efficiency: 57.68% of 5.00 GB



# Try TE loci: 5-1000 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=10G
TIME="6-00:00:00"
JOB="TEgenSim_5_1000"
sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ./log/${JOB}.out
#SBATCH -e ./log/${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_1000.py
python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_1000.py

EOF
Submitted batch job 4387977
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 04:23:34
CPU Efficiency: 99.79% of 04:24:07 core-walltime
Job Wall-clock time: 04:24:07
Memory Utilized: 5.05 GB
Memory Efficiency: 50.48% of 10.00 GB







# Jan 2024
# Try TE loci: 5-10 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=10G
TIME="06:00:00"
JOB="TEgenSim_5_10"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J ${JOB} -o ./log/${JOB}.out -e ./log/${JOB}.err --wrap "python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_10.py"

#Submitted batch job 4141298
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:00:28
#CPU Efficiency: 75.68% of 00:00:37 core-walltime
#Job Wall-clock time: 00:00:37
#Memory Utilized: 842.50 MB
#Memory Efficiency: 8.23% of 10.00 GB

THREADS=1
RAM=1G
TIME="06:00:00"
JOB="TEgenSim_nest_5_10"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J ${JOB} -o ./log/${JOB}.out -e ./log/${JOB}.err --wrap "python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_10_v2.py"

#Submitted batch job 4141323
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:00:04
#CPU Efficiency: 40.00% of 00:00:10 core-walltime
#Job Wall-clock time: 00:00:10
#Memory Utilized: 0.00 MB (estimated maximum)
#Memory Efficiency: 0.00% of 1.00 GB (1.00 GB/node)

### 2024-02-16 ###
### test if v2 nested script fix the bug
THREADS=1
RAM=1G
TIME="06:00:00"
JOB="TEgenSim_nest_5_10"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J ${JOB} -o ./log/${JOB}.out -e ./log/${JOB}.err --wrap "python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_10_v2.py"
Submitted batch job 4362023
Job ID: 4362023
Cluster: powerplant
User/Group: cflthc/cflthc
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:00:04
CPU Efficiency: 50.00% of 00:00:08 core-walltime
Job Wall-clock time: 00:00:08
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 1.00 GB (1.00 GB/node)


# Try TE loci: 5-100 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=10G
TIME="1-00:00:00"
JOB="TEgenSim_5_100"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J ${JOB} -o ./log/${JOB}.out -e ./log/${JOB}.err --wrap "python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_100.py"
# Submitted batch job 4140987
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:16:26
#CPU Efficiency: 99.80% of 00:16:28 core-walltime
#Job Wall-clock time: 00:16:28
#Memory Utilized: 1.07 GB
#Memory Efficiency: 10.72% of 10.00 GB

THREADS=1
RAM=1G
TIME="06:00:00"
JOB="TEgenSim_nest_5_100"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J ${JOB} -o ./log/${JOB}.out -e ./log/${JOB}.err --wrap "python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_100.py"
# Submitted batch job 3838093 # Completed; Job Wall-clock time: 00:23:55; Memory Utilized: 907.72 MB
#Submitted batch job 4148180
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:05:37
#CPU Efficiency: 98.54% of 00:05:42 core-walltime
#Job Wall-clock time: 00:05:42
#Memory Utilized: 768.39 MB
#Memory Efficiency: 75.04% of 1.00 GB

# Try TE loci: 5-500 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=10G
TIME="7-00:00:00"
JOB="TEgenSim_5_500"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J ${JOB} -o ./log/${JOB}.out -e ./log/${JOB}.err --wrap "python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_500.py"
# Submitted batch job 3836975 #Completed; Job Wall-clock time: 12:38:39; Memory Utilized: 1.86 GB
# Submitted batch job 4140946
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 11:01:57
#CPU Efficiency: 99.96% of 11:02:12 core-walltime
#Job Wall-clock time: 11:02:12
#Memory Utilized: 1.92 GB
#Memory Efficiency: 19.24% of 10.00 GB

THREADS=1
RAM=1G
TIME="3-00:00:00"
JOB="TEgenSim_nest_5_500"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J ${JOB} -o ./log/${JOB}.out -e ./log/${JOB}.err --wrap "python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_500.py"
#Submitted batch job 3868118; # Completed; Job Wall-clock time: 07:52:09; Memory Utilized: 2.93 GB
#Submitted batch job 4148191
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 01:18:52
#CPU Efficiency: 99.62% of 01:19:10 core-walltime
#Job Wall-clock time: 01:19:10
#Memory Utilized: 2.96 GB
#Memory Efficiency: 295.53% of 1.00 GB

# Try TE loci: 5-1000 per family
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator
THREADS=1
RAM=15G
TIME="7-00:00:00"
JOB="TEgenSim_5_1000"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J ${JOB} -o ./log/${JOB}.out -e ./log/${JOB}.err --wrap "python $TEgenomeSimulator/custom_genome_seq_TEs_DH_5_1000.py"
# Submitted batch job 3837508 #Completed; Job Wall-clock time: 1-22:06:04; Memory Utilized: 3.38 GB
#Submitted batch job 4140956
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 1-10:39:16
#CPU Efficiency: 99.99% of 1-10:39:26 core-walltime
#Job Wall-clock time: 1-10:39:26
#Memory Utilized: 3.31 GB
#Memory Efficiency: 22.05% of 15.00 GB

THREADS=1
RAM=1G
TIME="3-00:00:00"
JOB="TEgenSim_nest_5_1000"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J ${JOB} -o ./log/${JOB}.out -e ./log/${JOB}.err --wrap "python $TEgenomeSimulator/custom_genome_nest_TEs_DH_5_1000.py"
# Submitted batch job 3868222; # Completed; Job Wall-clock time: 1-07:20:50; Memory Utilized: 5.59 GB
#Submitted batch job 4155294
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 05:27:48
#CPU Efficiency: 99.85% of 05:28:17 core-walltime
#Job Wall-clock time: 05:28:17
#Memory Utilized: 5.17 GB
#Memory Efficiency: 517.48% of 1.00 GB