#!/bin/bash
# Navigate to the cloned repository of TEgenomeSimulator
ml apptainer
WRK=$(pwd)
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/pub/input
OUT=$WRK/pub/output
mkdir -p $IN $OUT

prefix=rice
genome=$IN/GCF_034140825.1_ASM3414082v1_genomic.fna
repeat=$IN/rice6.9.5.liban

JOB=TEGS_rice
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME --nodelist $NODE -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
# Submitted batch job 8791110
# Job ID: 8791110
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 03:34:36
# CPU Efficiency: 24.70% of 14:29:00 core-walltime
# Job Wall-clock time: 00:43:27
# Memory Utilized: 5.43 GB
# Memory Efficiency: 10.87% of 50.00 GB

cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

##############################################
# Weird high proportion of intact TE         #
# Updated TEgenomeSimulator.sif and rerun    #
##############################################
ml apptainer
WRK=$(pwd)
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/pub/input
OUT=../TEgenomeSimulator_pub/output
mkdir -p $IN $OUT

prefix=rice_run2
genome=$IN/GCF_034140825.1_ASM3414082v1_genomic.fna
repeat=$IN/rice6.9.5.liban

JOB=TEGS_rice_run2
TIME="1-00:00:00"
THREADS=20
MEM=50G
#NODE="aklppg31"
cd $OUT
OUT=$(pwd)
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
# Submitted batch job 8825560


cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'