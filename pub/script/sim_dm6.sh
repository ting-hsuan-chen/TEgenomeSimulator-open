#!/bin/bash
# Navigate to the cloned repository of TEgenomeSimulator
ml apptainer
WRK=$(pwd)
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/pub/input
OUT=$WRK/pub/output
mkdir -p $IN $OUT

# Take pseudo-chromosome only
cd $IN
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna |
grep -A1 -E '>NC_004354.4|>NT_033779.5|>NT_033778.4|>NT_037436.4|>NT_033777.3|>NC_004353.4|>NC_024512.1' > GCF_000001215.4_Release_6_chromosome_only.fna

prefix=dm6
genome=$IN/GCF_000001215.4_Release_6_chromosome_only.fna
repeat=$IN/Dmel-families.fa

JOB=TEGS_dm6
TIME="1-00:00:00"
THREADS=20
MEM=10G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME --nodelist $NODE -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
# Submitted batch job 8791194 # failed at fix_empty_seq.py. Try to keep only pseudo-chromosomes and exclude mitochondrion.
# Submitted batch job 8792680 # failed with error: OSError: [Errno 30] Read-only file system: '/workspace/cflthc/script/TEgenomeSimulator/pub/output/TEgenomeSimulator_dm6_result'
# Remove the result folder and rerun
# Submitted batch job 8795239
# Job ID: 8795239
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 00:15:02
# CPU Efficiency: 16.05% of 01:33:40 core-walltime
# Job Wall-clock time: 00:04:41
# Memory Utilized: 1.56 GB
# Memory Efficiency: 15.58% of 10.00 GB




cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
