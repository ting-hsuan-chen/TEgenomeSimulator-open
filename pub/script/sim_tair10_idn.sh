#!/bin/bash

##########################
## mean identity 90-100 ##
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

prefix=tair10_idn_90_100
minidn=90
maxidn=100
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"
# Submitted batch job 8498051
# Job ID: 8498051
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 00:01:49
# CPU Efficiency: 3.66% of 00:49:40 core-walltime
# Job Wall-clock time: 00:02:29
# Memory Utilized: 526.64 MB
# Memory Efficiency: 1.03% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########################
## mean identity 80-90 ##
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

prefix=tair10_idn_80_90
minidn=80
maxidn=90
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=1
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT"
# Submitted batch job 8503101


# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########################
## mean identity 70-80 ##
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

prefix=tair10_idn_70_80
minidn=70
maxidn=80
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=1
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT"
# Submitted batch job 8499044


# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########################
## mean identity 85-100 ##
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

prefix=tair10_idn_85_100
minidn=85
maxidn=100
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"
# Submitted batch job 8498214
# Job ID: 8498214
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 00:01:48
# CPU Efficiency: 3.58% of 00:50:20 core-walltime
# Job Wall-clock time: 00:02:31
# Memory Utilized: 672.22 MB
# Memory Efficiency: 1.31% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 80-95 ##
#########################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_80_95
minidn=80
maxidn=95
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"
# Submitted batch job 8498342
# Job ID: 8498342
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 00:02:09
# CPU Efficiency: 4.27% of 00:50:20 core-walltime
# Job Wall-clock time: 00:02:31
# Memory Utilized: 714.64 MB
# Memory Efficiency: 1.40% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 75-90 ##
#########################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_75_90
minidn=75
maxidn=90
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"
# Submitted batch job 8498468
# Job ID: 8498468
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 00:01:56
# CPU Efficiency: 5.04% of 00:38:20 core-walltime
# Job Wall-clock time: 00:01:55
# Memory Utilized: 550.22 MB
# Memory Efficiency: 1.07% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 70-85 ##
#########################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_70_85
minidn=70
maxidn=85
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"
# Submitted batch job 8498482
# ob ID: 8498482
# luster: powerplant
# ser/Group: cflthc/cflthc
# tate: COMPLETED (exit code 0)
# odes: 1
# ores per node: 20
# PU Utilized: 00:02:14
# PU Efficiency: 4.59% of 00:48:40 core-walltime
# ob Wall-clock time: 00:02:26
# emory Utilized: 602.37 MB
# emory Efficiency: 1.18% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########################
## mean identity 80-100 ##
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

prefix=tair10_idn_80_100
minidn=80
maxidn=100
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"
# Submitted batch job 8498609
# Job ID: 8498609
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 00:01:51
# CPU Efficiency: 5.05% of 00:36:40 core-walltime
# Job Wall-clock time: 00:01:50
# Memory Utilized: 552.08 MB
# Memory Efficiency: 1.08% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 75-95 ##
#########################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_75_95
minidn=75
maxidn=95
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT -t $THREADS"
# Submitted batch job 8498610
# Job ID: 8498610
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 00:02:10
# CPU Efficiency: 4.55% of 00:47:40 core-walltime
# Job Wall-clock time: 00:02:23
# Memory Utilized: 604.76 MB
# Memory Efficiency: 1.18% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 70-90 ##
#########################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_70_90
minidn=70
maxidn=90
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=1
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT"
# Submitted batch job 8498735
# Job ID: 8498735
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:01:46
# CPU Efficiency: 92.98% of 00:01:54 core-walltime
# Job Wall-clock time: 00:01:54
# Memory Utilized: 728.86 MB
# Memory Efficiency: 1.42% of 50.00 GB


# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#########################
## mean identity 70-95 ##
#########################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_idn_70_95
minidn=70
maxidn=95
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=1
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT"
# Submitted batch job 8498861
# Job ID: 8498861
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:01:38
# CPU Efficiency: 95.15% of 00:01:43 core-walltime
# Job Wall-clock time: 00:01:43
# Memory Utilized: 702.69 MB
# Memory Efficiency: 1.37% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

##########################
## mean identity 70-100 ##
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

prefix=tair10_idn_70_100
minidn=70
maxidn=100
min=5
max=100
genome=$OUT/TEgenomeSimulator_tair10_run5_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed
repeat=$OUT/TEgenomeSimulator_tair10_run5_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched

JOB=TEgenomeSimulator_idn_$minidn'_'$maxidn
TIME="1-00:00:00"
THREADS=1
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min --maxidn $maxidn --minidn $minidn -o $OUT"
# Submitted batch job 8498864
# Job ID: 8498864
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:01:52
# CPU Efficiency: 94.12% of 00:01:59 core-walltime
# Job Wall-clock time: 00:01:59
# Memory Utilized: 694.86 MB
# Memory Efficiency: 1.36% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'