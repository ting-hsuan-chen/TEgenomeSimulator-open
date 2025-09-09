#!/bin/bash
# Navigate to the cloned repository of TEgenomeSimulator
ml apptainer
WRK=$(pwd)
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/pub/input
OUT=$WRK/pub/output
mkdir -p $IN $OUT

prefix=maize
genome=$IN/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna
repeat=$IN/maizeTE11122019

JOB=TEGS_maize
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME --nodelist $NODE -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
# Submitted batch job 8790055
# Job ID: 8790055
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 1-11:53:04
# CPU Efficiency: 13.66% of 10-22:36:40 core-walltime
# Job Wall-clock time: 13:07:50
# Memory Utilized: 32.30 GB
# Memory Efficiency: 64.60% of 50.00 GB

cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
