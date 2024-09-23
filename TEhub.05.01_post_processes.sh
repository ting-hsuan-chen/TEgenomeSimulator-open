#!/bin/bash

# Move result to a different directory
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator/result
OUT=/workspace/cflthc/scratch/2022_Actinidia_TE/10.01_TEgenomeSimulator_output

mv sim_tair10_genome $OUT/sim_tair10_genome

# Index the result fasta files
ml samtools/1.16
OUT=/workspace/cflthc/scratch/2022_Actinidia_TE/10.01_TEgenomeSimulator_output

cd $OUT/sim_tair10_genome
samtools faidx tair10_genome_sequence_out_nest.fasta
samtools faidx tair10_repeat_sequence_out_nest.fasta

