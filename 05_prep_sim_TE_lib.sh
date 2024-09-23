#!/bin/bash

#Ting-Hsuan Chen
#2024-01-10
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
OUT=$WRK/09.10_BM_modules/curated_TElib
cd $OUT

# combine TE families to be included for simulation
# Copia, Gypsy, LTRsolo, LTRunknown, LINE, SINE, DTA, DTC, DTH, DTM DTT, Helitron, MITE
cat combined_curated_TE_lib_ATOSZM_cdhit_LTRcopia.fasta combined_curated_TE_lib_ATOSZM_cdhit_LTRgypsy.fasta \
combined_curated_TE_lib_ATOSZM_cdhit_LTRsolo.fasta combined_curated_TE_lib_ATOSZM_cdhit_LTRunknown.fasta \
combined_curated_TE_lib_ATOSZM_cdhit_LINE.fasta combined_curated_TE_lib_ATOSZM_cdhit_SINE.fasta \
combined_curated_TE_lib_ATOSZM_cdhit_DNAhat.fasta combined_curated_TE_lib_ATOSZM_cdhit_DNAcacta.fasta \
combined_curated_TE_lib_ATOSZM_cdhit_DNAharbinger.fasta combined_curated_TE_lib_ATOSZM_cdhit_DNAmudr.fasta \
combined_curated_TE_lib_ATOSZM_cdhit_DNAmariner.fasta combined_curated_TE_lib_ATOSZM_cdhit_RChelitron.fasta \
combined_curated_TE_lib_ATOSZM_cdhit_MITEstow.fasta combined_curated_TE_lib_ATOSZM_cdhit_MITEtourist.fasta \
> combined_curated_TE_lib_ATOSZM_selected.fasta

ml samtools/1.16
samtools faidx combined_curated_TE_lib_ATOSZM_selected.fasta

# soft link TElib fasta to $TEgenomeSimulator dir
ln -s $OUT/combined_curated_TE_lib_ATOSZM_selected.fasta $TEgenomeSimulator/combined_curated_TE_lib_ATOSZM_selected.fasta

# prepare the TElib_sim_list.table
#loci_count: 10 - 1000
#idn: 80 - 95
#sd: 1-20
#indels: 5-20 (as a proportion to total SNP = substitution + indel)
#tsd: use prior knowledge
#length: from .fai
#frag: 50 - 90
#nest: 0 - 30 (only for Copia and Gypsy)

# create another table that has fewer TE loci: 5-10 loci per family
ml conda
conda activate TEgenomeSimulator
cd $TEgenomeSimulator
python $TEgenomeSimulator/05_prep_sim_TE_lib_5_10.py

# create another table that has fewer TE loci: 5-100 loci per family
ml conda
conda activate TEgenomeSimulator
cd $TEgenomeSimulator
python $TEgenomeSimulator/05_prep_sim_TE_lib_5_100.py

# create another table that has fewer TE loci: 5-500 loci per family
ml conda
conda activate TEgenomeSimulator
cd $TEgenomeSimulator
python $TEgenomeSimulator/05_prep_sim_TE_lib_5_500.py

# create another table that has fewer TE loci: 5-1000 loci per family
ml conda
conda activate TEgenomeSimulator
cd $TEgenomeSimulator
python $TEgenomeSimulator/05_prep_sim_TE_lib_5_1000.py