#!/bin/bash

# Move result to a different directory
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
cd $TEgenomeSimulator/result
OUT=/workspace/cflthc/scratch/2022_Actinidia_TE/10.01_TEgenomeSimulator_output
#mv sim_random_genome $OUT/sim_random_genome
#mv sim_custom_genome $OUT/sim_custom_genome
mv sim_donghong_5_10_genome $OUT/sim_donghong_5_10_genome
mv sim_donghong_5_100_genome $OUT/sim_donghong_5_100_genome
mv sim_donghong_5_500_genome $OUT/sim_donghong_5_500_genome
mv sim_donghong_5_1000_genome $OUT/sim_donghong_5_1000_genome


# Index the result fasta files
ml samtools/1.16
OUT=/workspace/cflthc/scratch/2022_Actinidia_TE/10.01_TEgenomeSimulator_output

range=5_10
cd $OUT/sim_donghong_$range'_genome'
samtools faidx donghong_$range'_genome_sequence_out_nest.fasta'
samtools faidx donghong_$range'_repeat_sequence_out_nest.fasta'

range=5_100
cd $OUT/sim_donghong_$range'_genome'
samtools faidx donghong_$range'_genome_sequence_out_nest.fasta'
samtools faidx donghong_$range'_repeat_sequence_out_nest.fasta'

range=5_500
cd $OUT/sim_donghong_$range'_genome'
samtools faidx donghong_$range'_genome_sequence_out_nest.fasta'
samtools faidx donghong_$range'_repeat_sequence_out_nest.fasta'

range=5_1000
cd $OUT/sim_donghong_$range'_genome'
samtools faidx donghong_$range'_genome_sequence_out_nest.fasta'
samtools faidx donghong_$range'_repeat_sequence_out_nest.fasta'

