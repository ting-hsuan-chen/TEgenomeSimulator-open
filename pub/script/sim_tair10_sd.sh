#!/bin/bash

##########################
## scenario I 
## high idn 90-100 
## high sd 15-20
##########################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_sd_scenario1
minidn=90
maxidn=100
minsd=15
maxsd=20
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_scenario1
TIME="1-00:00:00"
THREADS=1
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -o $OUT"
# Submitted batch job 8503102
# Job ID: 8503102
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 02:21:22
# CPU Efficiency: 99.58% of 02:21:58 core-walltime
# Job Wall-clock time: 02:21:58
# Memory Utilized: 2.75 GB
# Memory Efficiency: 5.51% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########################
## scenario II
## low idn 70-80 
## high sd 15-20
##########################
ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_sd_scenario2
minidn=70
maxidn=80
minsd=15
maxsd=20
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_scenario2
TIME="1-00:00:00"
THREADS=1
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -o $OUT"
# Submitted batch job 8503103
# Job ID: 8503103
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 02:26:44
# CPU Efficiency: 99.63% of 02:27:17 core-walltime
# Job Wall-clock time: 02:27:17
# Memory Utilized: 2.60 GB
# Memory Efficiency: 5.20% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########################
## scenario III
## low idn 70-80 
## low sd 1-5
##########################
ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_sd_scenario3
minidn=70
maxidn=80
minsd=1
maxsd=5
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_scenario3
TIME="1-00:00:00"
THREADS=1
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -o $OUT"
# Submitted batch job 8503104
Job ID: 8503104
Cluster: powerplant
User/Group: cflthc/cflthc
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 02:26:41
CPU Efficiency: 99.60% of 02:27:16 core-walltime
Job Wall-clock time: 02:27:16
Memory Utilized: 2.60 GB
Memory Efficiency: 5.20% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########################
## scenario IV
## high idn 90-100 
## low sd 1-5
##########################
ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_sd_scenario4
minidn=90
maxidn=100
minsd=1
maxsd=5
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_scenario4
TIME="1-00:00:00"
THREADS=1
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -o $OUT"
# Submitted batch job 8503105
# Job ID: 8503105
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 02:22:51
# CPU Efficiency: 99.59% of 02:23:26 core-walltime
# Job Wall-clock time: 02:23:26
# Memory Utilized: 2.59 GB
# Memory Efficiency: 5.18% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'