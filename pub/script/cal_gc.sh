#!/bin/bash
# Navigate to the cloned repository of TEgenomeSimulator
WRK=$(pwd)
ml conda
conda activate hraaxt_emboss
cd $OUT/TEgenomeSimulator_tair10_ch1_m2_result 
gc=$(cat GCF_000001735.4_TAIR10.1_genomic_chr1.fna.masked.reformatted.nonTE.emptfixed | infoseq -auto -only -name -length -pgc stdin | awk 'NR>1{print $3}')
