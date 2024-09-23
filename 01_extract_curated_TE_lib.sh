#!/bin/bash

#Ting-Hsuan Chen
#2023-11-28

# Arabidopsis thaliana TE libarary from EDTA's github
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/09.10_BM_modules/curated_TElib/AT2
mkdir -p $OUT
cd $OUT
cp /workspace/cflthc/bin/EDTA2.1/EDTA/database/athrep.updated.nonredun.fasta ./

# Rice TE library
# The TE library used for benchmarking in EDTA's and RM2's original papers
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/09.10_BM_modules/curated_TElib/OS
mkdir -p $OUT
cd $OUT
cp /workspace/cflthc/bin/EDTA2.1/EDTA/database/rice6.9.5.liban ./

# Maize TE library
# The TE library used for benchmarking in EDTA's original paper
# LTR retrotransposons have been stiched as LTR-INT-LTR
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/09.10_BM_modules/curated_TElib/ZM
mkdir -p $OUT
cd $OUT
cp /workspace/cflthc/bin/EDTA2.1/EDTA/database/maizeTE11122019 ./

# Tomato TE library from the Sol Genomics Network
# seems to have too many TE sequences
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/09.10_BM_modules/curated_TElib/SL
mkdir -p $OUT
cd $OUT
wget https://solgenomics.net/ftp/tomato_genome/repeats/repeats_master.v5.fasta.gz
gzip -d repeats_master.v5.fasta.gz


