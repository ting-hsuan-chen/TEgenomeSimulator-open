#!/bin/bash

#Ting-Hsuan Chen
#2023-12-21

WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/09.10_BM_modules/curated_TElib
cd $OUT

ATlib=$OUT/AT2/athrep.updated.nonredun.LTRstitched.fasta
OSlib=$OUT/OS/rice6.9.5.liban.LTRstitched.fasta
ZMlib=$OUT/ZM/maizeTE11122019

cat $ATlib $OSlib $AMlib > combined_curated_TE_lib_ATOSZM.fasta

# Remove redudandt TE family seq
# Use the 80-80-80 rules proposed by Wicker et al. 2007 and adopted by RM2
# -c  sequence identity threshold, default 0.9
# -n  word_length, default 5, for thresholds 0.80 ~ 0.85
# -l  length of throw_away_sequences, default 10
# -aS alignment coverage for the shorter sequence, default 0.0; if set to 0.9, the alignment must covers 90% of the sequence
MEM=4000 # i.e. 4G
THREADS=4
TIME="01:00:00"
JOB="cdhit"
module load conda
conda activate cflthc_cdhit
sbatch -t $TIME --mem $MEM --cpus-per-task=$THREADS -J $JOB -o $JOB.out -e $JOB.err --wrap "cd-hit-est -i combined_curated_TE_lib_ATOSZM.fasta -o combined_curated_TE_lib_ATOSZM_cdhit -c 0.8 -n 5 -l 80 -aS 0.8 -M $MEM -T $THREADS"
#Submitted batch job 3731863

mv combined_curated_TE_lib_ATOSZM_cdhit combined_curated_TE_lib_ATOSZM_cdhit.fasta

conda deactivate
module unload conda
# grep ">" -c combined_curated_TE_lib_ATOSZM.fasta
# 2411
# grep ">" -c combined_curated_TE_lib_ATOSZM_cdhit.fasta
# 1918