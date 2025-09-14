#!/bin/bash
ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
tetools=/workspace/cflthc/bin/TETools/dfam-tetools.sif
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

# tair10 genome was downloaded from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/
# arabidopsis TE library was downloaded from https://github.com/oushujun/EDTA/blob/master/database/athrep.updated.nonredun.fasta
#cd $IN
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/GCF_000001735.4_TAIR10.1_genomic.fna ./
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/athrep.updated.nonredun.fasta ./

############################################
##remove centromeric and satellite repeat ##
############################################
ml seqkit/0.7.0
cd $IN
samtools faidx athrep.updated.nonredun.fasta
grep -E "centromeric|Satellite" athrep.updated.nonredun.fasta.fai | awk -F'\t' '{print $1}' > repeat_blocklist.txt
seqkit grep -v -f repeat_blocklist.txt athrep.updated.nonredun.fasta > athrep.updated.nonredun.noCenSatelli.fasta

#############
prefix=tair10
genome=$IN/GCF_000001735.4_TAIR10.1_genomic.fna
repeat=$IN/athrep.updated.nonredun.fasta

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=16
MEM=10G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 2 -S $tetools -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
# Submitted batch job 8422360
# Job ID: 8422360
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 16
# CPU Utilized: 00:20:49
# CPU Efficiency: 10.99% of 03:09:20 core-walltime
# Job Wall-clock time: 00:11:50
# Memory Utilized: 1.31 GB
# Memory Efficiency: 13.13% of 10.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


###########################################################################
## rerun after applying beta distribution for fragmentation distribution ##
###########################################################################

##########
## run2 ##
##########

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
tetools=/workspace/cflthc/bin/TETools/dfam-tetools.sif
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

# tair10 genome was downloaded from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/
# arabidopsis TE library was downloaded from https://github.com/oushujun/EDTA/blob/master/database/athrep.updated.nonredun.fasta
#cd $IN
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/GCF_000001735.4_TAIR10.1_genomic.fna ./
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/athrep.updated.nonredun.fasta ./

prefix=tair10_run2
genome=$IN/GCF_000001735.4_TAIR10.1_genomic.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
alpha=0.5
beta=0.7
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME --nodelist $NODE -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 2 -S $tetools -p $prefix -g $genome -r $repeat -r2 $repeat -a $alpha -b $beta -o $OUT -t $THREADS"
#Submitted batch job 8431692
#Job ID: 8431692
#Cluster: powerplant
#User/Group: cflthc/cflthc
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 20
#CPU Utilized: 00:18:35
#CPU Efficiency: 14.26% of 02:10:20 core-walltime
#Job Wall-clock time: 00:06:31
#Memory Utilized: 1.25 GB
#Memory Efficiency: 2.50% of 50.00 GB

JOB=TEgenomeSimulator
TIME="01:00:00"
THREADS=20
MEM=5G
NODE="aklppg31"
alpha=0.5
beta=0.7
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 2 -S $tetools -p $prefix -g $genome -r $repeat -r2 $repeat -a $alpha -b $beta -o $OUT -t $THREADS"
# Submitted batch job 8515852 # this run used the old fragmentation method (beta distribution)
# Job ID: 8515852
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 20
# CPU Utilized: 00:18:05
# CPU Efficiency: 14.16% of 02:07:40 core-walltime
# Job Wall-clock time: 00:06:23
# Memory Utilized: 1.26 GB
# Memory Efficiency: 25.12% of 5.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

##########
## run3 ##
##########

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
tetools=/workspace/cflthc/bin/TETools/dfam-tetools.sif
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

# tair10 genome was downloaded from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/
# arabidopsis TE library was downloaded from https://github.com/oushujun/EDTA/blob/master/database/athrep.updated.nonredun.fasta
#cd $IN
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/GCF_000001735.4_TAIR10.1_genomic.fna ./
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/athrep.updated.nonredun.fasta ./

prefix=tair10_run3
genome=$IN/GCF_000001735.4_TAIR10.1_genomic.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
alpha=0.7
beta=0.5
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME --nodelist $NODE -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 1 --to_mask -S $tetools -p $prefix -g $genome -r $repeat -r2 $repeat -a $alpha -b $beta -o $OUT -t $THREADS"
#Submitted batch job 8511580
#Job ID: 8511580
#Cluster: powerplant
#User/Group: cflthc/cflthc
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 20
#CPU Utilized: 00:15:07
#CPU Efficiency: 24.12% of 01:02:40 core-walltime
#Job Wall-clock time: 00:03:08
#Memory Utilized: 938.53 MB
#Memory Efficiency: 1.83% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

##########
## run4 ##
##########

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
tetools=/workspace/cflthc/bin/TETools/dfam-tetools.sif
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

# tair10 genome was downloaded from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/
# arabidopsis TE library was downloaded from https://github.com/oushujun/EDTA/blob/master/database/athrep.updated.nonredun.fasta
#cd $IN
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/GCF_000001735.4_TAIR10.1_genomic.fna ./
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/athrep.updated.nonredun.fasta ./

prefix=tair10_run4
genome=$IN/GCF_000001735.4_TAIR10.1_genomic.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
alpha=0.7
beta=0.5
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME --nodelist $NODE -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 2 -S $tetools -p $prefix -g $genome -r $repeat -r2 $repeat -a $alpha -b $beta -o $OUT -t $THREADS"
#Submitted batch job 8433592
Job ID: 8433592
Cluster: powerplant
User/Group: cflthc/cflthc
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:16:54
CPU Efficiency: 15.55% of 01:48:40 core-walltime
Job Wall-clock time: 00:05:26
Memory Utilized: 1.01 GB
Memory Efficiency: 2.01% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

##########
## run5 ##
##########

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
tetools=/workspace/cflthc/bin/TETools/dfam-tetools.sif
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

# tair10 genome was downloaded from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/
# arabidopsis TE library was downloaded from https://github.com/oushujun/EDTA/blob/master/database/athrep.updated.nonredun.fasta
#cd $IN
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/GCF_000001735.4_TAIR10.1_genomic.fna ./
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/athrep.updated.nonredun.fasta ./

prefix=tair10_run5
genome=$IN/GCF_000001735.4_TAIR10.1_genomic.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

JOB=TEgenomeSimulator
TIME="1-00:00:00"
THREADS=20
MEM=50G
NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME --nodelist $NODE -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 2 -S $tetools -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
#Submitted batch job 8433620
#Job ID: 8433620
#Cluster: powerplant
#User/Group: cflthc/cflthc
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 20
#CPU Utilized: 00:17:05
#CPU Efficiency: 16.02% of 01:46:40 core-walltime
#Job Wall-clock time: 00:05:20
#Memory Utilized: 983.29 MB
#Memory Efficiency: 1.92% of 50.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

##########
## run6 ##
##########
# had fixed wrong labeling of integrity in gff3
# the actual fragmentation step was not affected

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
tetools=/workspace/cflthc/bin/TETools/dfam-tetools.sif
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

# tair10 genome was downloaded from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/
# arabidopsis TE library was downloaded from https://github.com/oushujun/EDTA/blob/master/database/athrep.updated.nonredun.fasta
#cd $IN
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/GCF_000001735.4_TAIR10.1_genomic.fna ./
#cp -v /workspace/cflthc/script/TEgenomeSimulator_test/input/athrep.updated.nonredun.fasta ./

prefix=tair10_run6
genome=$IN/GCF_000001735.4_TAIR10.1_genomic.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

JOB=TEgenomeSimulator
TIME="01:00:00"
THREADS=20
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 2 -S $tetools -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
#Submitted batch job 8513117
#Job ID: 8513117
#Cluster: powerplant
#User/Group: cflthc/cflthc
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 20
#CPU Utilized: 00:19:12
#CPU Efficiency: 10.85% of 02:57:00 core-walltime
#Job Wall-clock time: 00:08:51
#Memory Utilized: 932.64 MB
#Memory Efficiency: 18.22% of 5.00 GB

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

##########
## run7 ##
##########
# had fixed integrity > 1 
# test if the program runs correctly

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
tetools=/workspace/cflthc/bin/TETools/dfam-tetools.sif
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_run7
genome=$IN/GCF_000001735.4_TAIR10.1_genomic.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

JOB=TEgenomeSimulator
TIME="01:00:00"
THREADS=20
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 2 -S $tetools -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
#Submitted batch job 8559090 #completed

# indexing
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'


##########
## run8 ##
##########
# Test TEGenomeSimulator.sif 

#ml singularity/3
ml apptainer/1.1

TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator_v1.0.0_apptainer/TEgenomeSimulator_v1.0.0.sif
#TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/tair10
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

prefix=tair10_run8
genome=$IN/GCF_000001735.4_TAIR10.1_genomic.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

JOB=TEgenomeSimulator
TIME="01:00:00"
THREADS=20
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
#Submitted batch job 8560600 #completed
#Job ID: 8560600
Cluster: powerplant
User/Group: cflthc/cflthc
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:20:09
CPU Efficiency: 12.41% of 02:42:20 core-walltime
Job Wall-clock time: 00:08:07
Memory Utilized: 1.32 GB
Memory Efficiency: 26.42% of 5.00 GB

# indexing
module load samtools
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'