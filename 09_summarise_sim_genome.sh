#!/bin/bash
#Ting-Hsuan Chen
# Example:
module load R/4.2.3-foss-2022b-PFR
SCRIPT=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output
cd $SCRIPT
Rscript -e 'library(rmarkdown); render("09_summarise_sim_genome_DH_5_10.Rmd", github_document())'
# Alternatively, open up R studio, nevigating to the Rmd file and click "knit"

cp *.pdf $OUT