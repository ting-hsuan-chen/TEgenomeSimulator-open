#!/bin/bash

##########################
### alpha 0.5 beta 0.5 ###
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

prefix=tair10_itg_1
min=5
max=100
alpha=0.5
beta=0.5
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_itg_1
TIME="00:10:00"
THREADS=1
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -a $alpha -b $beta -o $OUT"
# Submitted batch job 8515120
# Job ID: 8515120
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:02:22
# CPU Efficiency: 87.65% of 00:02:42 core-walltime
# Job Wall-clock time: 00:02:42
# Memory Utilized: 568.04 MB
# Memory Efficiency: 11.09% of 5.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

###########################
### alpha 0.5 beta 0.75 ###
###########################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_itg_2
min=5
max=100
alpha=0.5
beta=0.75
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_itg_2
TIME="00:10:00"
THREADS=1
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -a $alpha -b $beta -o $OUT"
# Submitted batch job 8515121
# Job ID: 8515121
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:02:10
# CPU Efficiency: 87.25% of 00:02:29 core-walltime
# Job Wall-clock time: 00:02:29
# Memory Utilized: 498.22 MB
# Memory Efficiency: 9.73% of 5.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

########################
### alpha 0.5 beta 1 ###
########################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_itg_3
min=5
max=100
alpha=0.5
beta=1
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_itg_3
TIME="00:10:00"
THREADS=1
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -a $alpha -b $beta -o $OUT"
# Submitted batch job 8515122
# Job ID: 8515122
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:02:04
# CPU Efficiency: 86.71% of 00:02:23 core-walltime
# Job Wall-clock time: 00:02:23
# Memory Utilized: 514.92 MB
# Memory Efficiency: 10.06% of 5.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

###########################
### alpha 0.75 beta 0.5 ###
###########################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_itg_4
min=5
max=100
alpha=0.75
beta=0.5
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_itg_4
TIME="00:10:00"
THREADS=1
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -a $alpha -b $beta -o $OUT"
# Submitted batch job 8515123
# Job ID: 8515123
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:02:37
# CPU Efficiency: 87.22% of 00:03:00 core-walltime
# Job Wall-clock time: 00:03:00
# Memory Utilized: 551.24 MB
# Memory Efficiency: 10.77% of 5.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


########################
### alpha 1 beta 0.5 ###
########################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_itg_5
min=5
max=100
alpha=1
beta=0.5
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_itg_5
TIME="00:10:00"
THREADS=1
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -a $alpha -b $beta -o $OUT"
# Submitted batch job 8515124
# Job ID: 8515124
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:02:53
# CPU Efficiency: 89.18% of 00:03:14 core-walltime
# Job Wall-clock time: 00:03:14
# Memory Utilized: 649.02 MB
# Memory Efficiency: 12.68% of 5.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'