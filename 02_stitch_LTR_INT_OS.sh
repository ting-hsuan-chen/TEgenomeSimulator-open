#!/bin/bash

#Ting-Hsuan Chen
#2023-12-18
ml conda
conda activate TEgenomeSimulator

WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/09.10_BM_modules/curated_TElib/OS
cd $OUT
#cp /workspace/cflthc/bin/EDTA2.1/EDTA/database/athrep.updated.nonredun.fasta ./

# Grep LTR retrotransposon
grep -A1 'LTR' rice6.9.5.liban > rice6.9.5.liban.LTRretro.fasta

# Create a header list
grep '>' rice6.9.5.liban.LTRretro.fasta > header.list

# Extract the header of LTR ('LTR') and INT ('I') and then use sed to extract the common TE id
grep -e 'LTR#' -e '_LTR_' header.list | sed -e 's/^>//g' -e 's/ .*//g' > header_LTR.list
grep -e 'INT#' -e '_INT_' header.list | sed -e 's/^>//g' -e 's/ .*//g' > header_INT.list
sed -e 's/LTR#.*//g' -e 's/_LTR_.*//g' -e 's/_$//g' -e 's/-$//g' header_LTR.list > header_LTR_trim.list
sed -e 's/INT#.*//g' -e 's/_INT_.*//g' -e 's/_$//g' -e 's/-$//g' header_INT.list > header_INT_trim.list

# Run 02_stitch_LTR_INT.py to find exact match, stitch LTR-INT-LTR and write fasta file
# and repleace the LTR-INT-separated entries with the stitched sequence
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/09_benchmarking/TEgenomeSimulator
original_te_fa='rice6.9.5.liban'
ltr_retro_fa='rice6.9.5.liban.LTRretro.fasta'
out_file_fa_name='rice6.9.5.liban.LTRstitched.fasta'
python $TEgenomeSimulator/02_stitch_LTR_INT.py $original_te_fa $ltr_retro_fa $out_file_fa_name

rm -f header* *matched.txt

