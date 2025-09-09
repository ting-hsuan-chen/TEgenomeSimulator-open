#!/bin/bash

######################
## copy number 1-10 ##
######################

ml conda
ml singularity/3
ml samtools/1.20
conda activate TEgenomeSimulator
TEgenomeSimulator=/workspace/cflthc/script/TEgenomeSimulator/TEgenomeSimulator/TEgenomeSimulator.py
WRK=/workspace/cflthc/scratch/TEGE/TEgenomeSimulator/mode0
IN=$WRK/input
OUT=$WRK/output
mkdir -p $IN $OUT

cd $IN
echo "chr1,10000,35" > random_genome_chr_index_10k.csv
echo "chr1,100000,35" > random_genome_chr_index_100k.csv
echo "chr1,1000000,35" > random_genome_chr_index_1000k.csv
echo "chr1,10000000,35" > random_genome_chr_index_10000k.csv
echo "chr1,100000000,35" > random_genome_chr_index_100000k.csv

echo "chr1,100000000,35" > random_genome_chr_index_1000000k.csv
echo "chr2,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr3,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr4,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr5,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr6,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr7,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr8,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr9,100000000,35" >> random_genome_chr_index_1000000k.csv
echo "chr10,100000000,35" >> random_genome_chr_index_1000000k.csv

cp /workspace/cflthc/script/TEgenomeSimulator/test/input/combined_curated_TE_lib_ATOSZM_selected.fasta.gz ./
gzip -d combined_curated_TE_lib_ATOSZM_selected.fasta.gz
repeat=$IN/combined_curated_TE_lib_ATOSZM_selected.fasta

min=1
max=10

for chr in 10 100 1000 10000 100000
do 
indcsv=$IN/random_genome_chr_index_${chr}k.csv
prefix=m0_chr${chr}k
JOB=TEgenomeSimulator_m0_chr${chr}k
TIME="1-00:00:00"
THREADS=1
MEM=30G
NODE="aklppg31"
cd $OUT
#sbatch --cpus-per-task=$THREADS --nodelist $NODE --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 0 -c $indcsv -p $prefix -r $repeat -m $max -n $min -o $OUT -t $THREADS"
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 0 -c $indcsv -p $prefix -r $repeat -m $max -n $min -o $OUT -t $THREADS"
done
#Submitted batch job 8757347
#Submitted batch job 8757348
#Submitted batch job 8757349
#Submitted batch job 8757350
#Submitted batch job 8757351

#Submitted batch job 8812916
#Submitted batch job 8812917

chr=1000000
indcsv=$IN/random_genome_chr_index_${chr}k.csv
prefix=m0_chr${chr}k
JOB=TEgenomeSimulator_m0_chr${chr}k
TIME="1-00:00:00"
THREADS=1
MEM=100G
#NODE="aklppg31"
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 0 -c $indcsv -p $prefix -r $repeat -m $max -n $min -o $OUT -t $THREADS"
#Submitted batch job 8757376
#Job ID: 8757376
#Cluster: powerplant
#User/Group: cflthc/cflthc
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:03:52
#CPU Efficiency: 84.06% of 00:04:36 core-walltime
#Job Wall-clock time: 00:04:36
#Memory Utilized: 12.62 GB
#Memory Efficiency: 12.62% of 100.00 GB

# indexing
for chr in 10 100 1000 10000 100000
do 
prefix=m0_chr${chr}k
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'
done

chr=1000000
prefix=m0_chr${chr}k
cd $OUT/TEgenomeSimulator_$prefix'_result'
samtools faidx $prefix'_genome_sequence_out_final.fasta'
samtools faidx $prefix'_repeat_sequence_out_final.fasta'

#Job ID: 8757347
#Cluster: powerplant
#User/Group: cflthc/cflthc
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:00:07
#CPU Efficiency: 22.58% of 00:00:31 core-walltime
#Job Wall-clock time: 00:00:31
#Memory Utilized: 163.78 MB
#Memory Efficiency: 0.53% of 30.00 GB

#Job ID: 8757348
#Cluster: powerplant
#User/Group: cflthc/cflthc
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:00:07
#CPU Efficiency: 22.58% of 00:00:31 core-walltime
#Job Wall-clock time: 00:00:31
#Memory Utilized: 162.15 MB
#Memory Efficiency: 0.53% of 30.00 GB

#Job ID: 8757349
#Cluster: powerplant
#User/Group: cflthc/cflthc
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:00:08
#CPU Efficiency: 25.81% of 00:00:31 core-walltime
#Job Wall-clock time: 00:00:31
#Memory Utilized: 168.02 MB
#Memory Efficiency: 0.55% of 30.00 GB

#Job ID: 8757350
#Cluster: powerplant
#User/Group: cflthc/cflthc
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:00:10
#CPU Efficiency: 29.41% of 00:00:34 core-walltime
#Job Wall-clock time: 00:00:34
#Memory Utilized: 195.22 MB
#Memory Efficiency: 0.64% of 30.00 GB

#Job ID: 8757351
#Cluster: powerplant
#User/Group: cflthc/cflthc
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:00:50
#CPU Efficiency: 67.57% of 00:01:14 core-walltime
#Job Wall-clock time: 00:01:14
#Memory Utilized: 3.83 GB
#Memory Efficiency: 12.76% of 30.00 GB