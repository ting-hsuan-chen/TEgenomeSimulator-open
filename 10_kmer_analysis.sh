#!/bin/bash
#Ting-Hsuan Chen
#2024-03-05
#Show kmer profile

module load jellyfish/2.2.10
loci_range="5_100"
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/kmer_analysis
genomename="donghong_"$loci_range
genome=$WRK/10.01_TEgenomeSimulator_output/"sim_donghong_"$loci_range"_genome"/donghong_$loci_range"_genome_sequence_out_nest.fasta"

JOB="jellyfish_bc"
RAM=16G
THREADS=20
TIME="06:00:00"
mkdir -p $OUT
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish bc -m 21 -s 100G -t 20 -o $genomename.bc $genome"
#Submitted batch job 4483384
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:11:21
CPU Efficiency: 8.41% of 02:15:00 core-walltime
Job Wall-clock time: 00:06:45
Memory Utilized: 260.77 GB
Memory Efficiency: 1629.83% of 16.00 GB

JOB="jellyfish_count"
RAM=280G
THREADS=16
TIME="06:00:00"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish count -m 21 -s 100G -t 16 --bc $genomename.bc $genome"
Submitted batch job 4483496
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 16
CPU Utilized: 00:48:28
CPU Efficiency: 9.04% of 08:56:00 core-walltime
Job Wall-clock time: 00:33:30
Memory Utilized: 588.46 GB
Memory Efficiency: 210.16% of 280.00 GB



JOB="jellyfish_hist"
RAM=30G
THREADS=16
TIME="06:00:00"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish histo -t 16 mer_counts.jf"
Submitted batch job 4691444
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 16
CPU Utilized: 00:00:01
CPU Efficiency: 1.04% of 00:01:36 core-walltime
Job Wall-clock time: 00:00:06
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 30.00 GB (30.00 GB/node)


jellyfish bc -m 25 -s 100G -t 16 -o homo_sapiens.bc homo_sapiens.fa
jellyfish count -m 25 -s 3G -t 16 --bc homo_sapiens.bc homo_sapiens.fa


jellyfish histo -t 10 WDR/reads.histo
jellyfish histo 20and1.jf









### The original donghong genome

WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/kmer_analysis
genomename="donghong_original"
genome=$WRK/09.02_EDTA/ref/reformated/Donghong.chromosomes.only.fa

JOB="kmer"
RAM=16G
THREADS=20
TIME="06:00:00"
mkdir -p $OUT/$genomename
cd $OUT/$genomename
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish count -m 21 -s 100G --bf-size 100G -t 20 $genome"
Submitted batch job 4483294



### Didn't find this useful
### Deleted kmer output 2024-08-01 to release space in pP
# Clean up related scratch folder in workspace
DATE=20240801
cd /workspace/cflthc/scratch/2022_Actinidia_TE/10.01_TEgenomeSimulator_output
sbatch -t 1-00:00:00 --mem 10G -J rm_kmer_analysis_${DATE} -o rm_kmer_analysis_${DATE}.out -e rm_kmer_analysis_${DATE}.err --wrap "rm -rf kmer_analysis"
# Submitted batch job 6858192
