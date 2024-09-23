#!/bin/bash

#Ting-Hsuan Chen
#2024-03-18

OUT=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_ref
cd $OUT

##############
### genome ###
##############
# Arabidopsis genome
wget https://www.repeatmasker.org/~cgoubert/GARLIC24/inputs/AT/AthaGenome.fa
    # index
ml samtools/1.16
samtools faidx AthaGenome.fa


#AT=/workspace/cflthc/database/Arabidopsis/Arabidopsis_thaliana/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna
#ln -s $AT tair10.fa

# Drosophila genome
wget https://www.repeatmasker.org/~cgoubert/GARLIC24/inputs/dm6/dm6.fa
    # grep chromosomes only
    # some mysterious line would be capture if grep them all together
    # so need to seperate them into three files by using "grep" three times
    # before combining them together
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' dm6.fa |\
grep -A1 "chr2L$\|chr2R$\|chr3L$\|chr3R$\|chr4$" > dm6_temp.fa
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' dm6.fa |\
grep -A1 "chrX$" > dm6_temp.x.fa
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' dm6.fa |\
grep -A1 "chrY$" > dm6_temp.y.fa
cat dm6_temp.fa dm6_temp.x.fa dm6_temp.y.fa > dm6_chromosome_only.fa
rm -f dm6_temp.*

# Human genome
wget https://www.repeatmasker.org/~cgoubert/GARLIC24/inputs/hg38/hg38.simple.fa

##############
##### TE #####
##############

# Arabidopsis thaliana TE libarary from Clement github
wget https://www.repeatmasker.org/~cgoubert/GARLIC24/inputs/AT/athaTEref_ru1.classified.fasta
#ln -s /workspace/cflthc/bin/EDTA2.1/EDTA/database/athrep.updated.nonredun.fasta tair10_te.fa
    # index
ml samtools/1.16
samtools faidx athaTEref_ru1.classified.fasta


# curated TE library from Dfam
wget https://dfam.org/releases/Dfam_3.8/families/Dfam_curatedonly.embl.gz
gzip -d Dfam_curatedonly.embl.gz

SCRIPT=/workspace/cflthc/bin/embl_to_fasta
input_embl=$OUT/Dfam_curatedonly.embl
output_fa=Dfam_curatedonly.fa
output_dir=$OUT

# use python to convert .embl to .fa
ml conda
conda activate TEgenomeSimulator
python

#/usr/bin/env python
```python
from Bio import SeqIO
import sys
infile = "/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_ref/Dfam_curatedonly.embl"
intype = "embl"
outfile = "/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_ref/Dfam_curatedonly.fasta"
outtype = "fasta"
SeqIO.convert(infile, intype, outfile, outtype)
exit()
```

    # index
ml samtools/1.16
samtools faidx Dfam_curatedonly.fasta

    # create a file for collecting the full headers in Dfam
grep ">" Dfam_curatedonly.fasta | sed 's/>//g' > Dfam_curatedonly.fasta.header

    # combine .fai and the full header
paste Dfam_curatedonly.fasta.fai Dfam_curatedonly.fasta.header > Dfam_curatedonly.fasta.fai.fullheader