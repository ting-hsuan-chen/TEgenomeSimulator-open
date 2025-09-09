#!/bin/bash

######################
## copy number 1-10 ##
######################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_cn_1_10
min=1
max=10
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"
# Submitted batch job 8498005
# Job ID: 8498005
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 00:00:21
# CPU Efficiency: 3.50% of 00:10:00 core-walltime
# Job Wall-clock time: 00:00:30
# Memory Utilized: 218.34 MB
# Memory Efficiency: 0.43% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


#######################
## copy number 5-100 ##
#######################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_cn_5_100
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"
# Submitted batch job 8498013
# Job ID: 8498013
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 00:02:04
# CPU Efficiency: 4.59% of 00:45:00 core-walltime
# Job Wall-clock time: 00:02:15
# Memory Utilized: 709.72 MB
# Memory Efficiency: 1.39% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#######################
## copy number 5-500 ##
#######################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_cn_5_500
min=5
max=500
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"
# Submitted batch job 8498015
# Job ID: 8498015
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 00:29:35
# CPU Efficiency: 4.89% of 10:04:40 core-walltime
# Job Wall-clock time: 00:30:14
# Memory Utilized: 1.53 GB
# Memory Efficiency: 3.07% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

########################
## copy number 5-1000 ##
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

prefix=tair10_cn_5_1000
min=5
max=1000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"
# Submitted batch job 8498016
# Job ID: 8498016
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 02:07:45
# CPU Efficiency: 4.95% of 1-18:59:40 core-walltime
# Job Wall-clock time: 02:08:59
# Memory Utilized: 2.60 GB
# Memory Efficiency: 5.21% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

########################
## copy number 5-2000 ##
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

prefix=tair10_cn_5_2000
min=5
max=2000
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $OUT -t $THREADS"
# Submitted batch job 8498017
# Job ID: 8498017
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 09:48:38
# CPU Efficiency: 4.97% of 8-05:29:40 core-walltime
# Job Wall-clock time: 09:52:29
# Memory Utilized: 6.04 GB
# Memory Efficiency: 12.08% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'