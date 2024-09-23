#!/bin/bash
#Ting-Hsuan Chen
#2024-04-04
#Show GC content

module load emboss/6.6.0
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/GC_analysis
genomename="tair10_sim"
genome=$WRK/10.01_TEgenomeSimulator_output/sim_tair10_genome/tair10_genome_sequence_out_nest.fasta

JOB="emboss"
RAM=1G
THREADS=1
TIME="06:00:00"
mkdir -p $OUT/tair10_sim
cd $OUT/tair10_sim
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "infoseq -auto -only -name -length -pgc -outfile tair10_sim_GC.txt $genome"
#Submitted batch job 5258775
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:00:03
CPU Efficiency: 33.33% of 00:00:09 core-walltime
Job Wall-clock time: 00:00:09
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 1.00 GB (1.00 GB/node)


### The original tair10 genome
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/GC_analysis
genomename="tair10_ori"
genome=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_ref/AthaGenome.fa

JOB="emboss"
RAM=1G
THREADS=1
TIME="06:00:00"
mkdir -p $OUT/tair10_ori
cd $OUT/tair10_ori
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "infoseq -auto -only -name -length -pgc -outfile tair10_ori_GC.txt $genome"
#Submitted batch job 5197624
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:00:02
CPU Efficiency: 50.00% of 00:00:04 core-walltime
Job Wall-clock time: 00:00:04
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 1.00 GB (1.00 GB/node)

