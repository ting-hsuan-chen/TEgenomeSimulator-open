#!/bin/bash
# Navigate to the cloned repository of TEgenomeSimulator
ml apptainer
WRK=$(pwd)
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/pub/input
OUT=$WRK/pub/output
mkdir -p $IN $OUT

prefix=danio
genome=$IN/GCF_049306965.1_GRCz12tu_genomic.fna
repeat=$IN/danRer10-families.fa

JOB=TEGS_danio
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME --nodelist $NODE -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
# Submitted batch job 8790057
# Job ID: 8790057
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 1-21:33:59
# CPU Efficiency: 18.51% of 10-06:12:20 core-walltime
# Job Wall-clock time: 12:18:37
# Memory Utilized: 32.86 GB
# Memory Efficiency: 65.71% of 50.00 GB

cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
