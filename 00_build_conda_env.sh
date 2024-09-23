#!/bin/bash

#Ting-Hsuan Chen
#2023-12-18

module load conda
conda create --name TEgenomeSimulator python=3.9
conda activate TEgenomeSimulator
pip3 install biopython pyyaml numpy
pip install pandas