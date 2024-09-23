#!/bin/bash

#Ting-Hsuan Chen
#2024-03-19
ml singularity/3
WRK=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
REF=$WRK/model_species_ref
OUT=$WRK/model_species_RM
TOOL=/workspace/cflthc/bin/TETools
JOB=RepeatMasker

mkdir -p $OUT
cd $OUT

#####################
## RepeatMasker #####
## Mask repetitive ##
## regions in 'X's ##
#####################

# Arabidopsis 
THREADS=6
RAM=1G
TIME="06:00:00"
mkdir -p $OUT/tair10
cd $OUT/tair10
genome=AthaGenome.fa
telib=$REF/athaTEref_ru1.classified.fasta
ln -s $REF/$genome $genome

sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $OUT/$JOB.out -e $OUT/$JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -a -x -q -no_is -norna -nolow -div 40 -lib $telib -cutoff 225 $genome"
# Submitted batch job 4867846
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 6
#CPU Utilized: 00:15:46
#CPU Efficiency: 34.58% of 00:45:36 core-walltime
#Job Wall-clock time: 00:07:36
#Memory Utilized: 445.26 MB
#Memory Efficiency: 43.48% of 1.00 GB

# convert multi-line fasta to single-line fasta
MASKED='AthaGenome.fa.masked'
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $MASKED | tail -n +2 > temp.fa && yes | mv temp.fa $MASKED

# Remove TE nucleotides masked as 'X's
OUTPUT='AthaGenome.fa.nonTE'
sed 's/X//g' $MASKED > $OUTPUT

# create index
ml samtools/1.16
samtools faidx $OUTPUT


#------------------
# Drosophila genome
mkdir -p $OUT/dm6
cd $OUT/dm6
genome=dm6.fa
telib=$REF/Dfam_curatedonly.fasta
ln -s $REF/$genome $genome

sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $OUT/$JOB.out -e $OUT/$JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -a -x -q -no_is -norna -nolow -div 40 -lib $telib -cutoff 225 $genome"
# Submitted batch job 4855074
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 12
#CPU Utilized: 05:17:42
#CPU Efficiency: 90.67% of 05:50:24 core-walltime
#Job Wall-clock time: 00:29:12
#Memory Utilized: 661.66 MB
#Memory Efficiency: 5.38% of 12.00 GB

# convert multi-line fasta to single-line fasta
#MASKED='dm6.fa.masked'
#awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $MASKED | tail -n +2 > temp.fa && yes | mv temp.fa $MASKED

#------------------
# Drosophila genome chromosome only
THREADS=12
RAM=2G
TIME="06:00:00"
mkdir -p $OUT/dm6
cd $OUT/dm6
genome=dm6_chromosome_only.fa
telib=$REF/Dfam_curatedonly.fasta
ln -s $REF/$genome $genome

sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $OUT/$JOB.out -e $OUT/$JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -a -x -q -no_is -norna -nolow -div 40 -lib $telib -cutoff 225 $genome"
# Submitted batch job 4919523
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 12
#CPU Utilized: 04:40:58
#CPU Efficiency: 93.34% of 05:01:00 core-walltime
#Job Wall-clock time: 00:25:05
#Memory Utilized: 1.13 GB
#Memory Efficiency: 56.45% of 2.00 GB

# convert multi-line fasta to single-line fasta
MASKED='dm6_chromosome_only.fa.masked'
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $MASKED | tail -n +2 > temp.fa && yes | mv temp.fa $MASKED



#------------------
# Human genome
THREADS=24
RAM=20G
TIME="3-00:00:00"
mkdir -p $OUT/hg38
cd $OUT/hg38
genome=hg38.simple.fa
telib=$REF/Dfam_curatedonly.fasta
ln -s $REF/$genome $genome

sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $OUT/$JOB.out -e $OUT/$JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -a -x -q -no_is -norna -nolow -div 40 -lib $telib -cutoff 225 $genome"
# Submitted batch job 4855408
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 24
#CPU Utilized: 11-07:39:33
#CPU Efficiency: 71.85% of 15-18:04:00 core-walltime
#Job Wall-clock time: 15:45:10
#Memory Utilized: 54.69 GB
#Memory Efficiency: 273.45% of 20.00 GB

# convert multi-line fasta to single-line fasta
MASKED='hg38.simple.fa.masked'
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $MASKED | tail -n +2 > temp.fa && yes | mv temp.fa $MASKED

# Remove TE nucleotides masked as 'X's
OUTPUT='hg38.simple.fa.nonTE'
sed 's/X//g' $MASKED > $OUTPUT

# create index
ml samtools/1.16
samtools faidx $OUTPUT