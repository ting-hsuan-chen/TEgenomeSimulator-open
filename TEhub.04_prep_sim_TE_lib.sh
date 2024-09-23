#!/bin/bash

#Ting-Hsuan Chen
#2024-03-27
TEgenomeSimulator=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator

# create table for tair10
ml conda
conda activate TEgenomeSimulator
cd $TEgenomeSimulator
python $TEgenomeSimulator/TEhub.04_prep_sim_TE_lib_tair10.py

