#!/bin/bash


# genome
# Extract chr1 of TAIR10
WRK=$(pwd)
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna |
grep -A1 -E 'chromosome 1' > $WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna

# extract chr1 TE annotation for visulisation
cd $WRK/pub/input
fullRMout=../output/TEgenomeSimulator_tair10_ch1_m2_result/repeatmasker/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.out
cp $fullRMout .
    # convert .out to .gff
in_rmout=GCF_000001735.4_TAIR10.1_genomic_chr1.fna.out
in_fai=athrep.updated.nonredun.noCenSatelli.fasta.stitched.fai
out_gff=GCF_000001735.4_TAIR10.1_genomic_chr1.fna.out.gff
col2label="RepeatMasker"
$WRK/pub/script/RMout2gff3v2.sh $in_rmout $in_fai $out_gff $col2label

# repeat
# # cp RM.out (no need)
# cd $WRK/pub/input
# cp ../../../../scratch/TEGE/TEgenomeSimulator/tair10/output/TEgenomeSimulator_tair10_run6_result/repeatmasker/GCF_000001735.4_TAIR10.1_genomic.fna.out ./
# cp RM.align
cp ../../../../scratch/TEGE/TEgenomeSimulator/tair10/output/TEgenomeSimulator_tair10_run6_result/repeatmasker/GCF_000001735.4_TAIR10.1_genomic.fna.align ./

# TRF.out
ml trf/4.07
cd $WRK/pub/input
trf GCF_000001735.4_TAIR10.1_genomic_chr1.fna 2 5 7 80 10 50 2000 -h -d | trf -ngs -

JOB=TRF_tair10
TIME="3-00:00:00"
THREADS=1
MEM=15G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "trf $WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna 2 5 7 80 10 50 2000 -h -d | trf -ngs -"
# Submitted batch job 8797341

# Convert TRF output .dat file to the format required for Garlic
../script/trfConverter4Garlic.sh GCF_000001735.4_TAIR10.1_genomic_chr1.fna.2.5.7.80.10.50.2000.dat > GCF_000001735.4_TAIR10.1_genomic_chr1.fna.2.5.7.80.10.50.2000.dat.reformated
../script/trfConverter4Garlic.sh GCF_000001735.4_TAIR10.1_genomic.fna.2.5.7.80.10.50.2000.dat > GCF_000001735.4_TAIR10.1_genomic.fna.2.5.7.80.10.50.2000.dat.reformated

# JOB=TRF_tair10_chr1
# TIME="06:00:00"
# THREADS=1
# MEM=15G
# sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "trf GCF_000001735.4_TAIR10.1_genomic_chr1.fna 2 5 7 80 10 50 2000 -h -d"
# Submitted batch job 8793320 # failed without informative error message # use `-ngs`
# Submitted batch job 8795571 # completed but no output file; remove `-d` and retry
# Job ID: 8795571
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:01:32
# CPU Efficiency: 97.87% of 00:01:34 core-walltime
# Job Wall-clock time: 00:01:34
# Memory Utilized: 372.39 MB
# Memory Efficiency: 3.64% of 10.00 GB

# Submitted batch job 8795573 # completed but no output file; back to use `-d`
# Job ID: 8795573
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:01:32
# CPU Efficiency: 98.92% of 00:01:33 core-walltime
# Job Wall-clock time: 00:01:33
# Memory Utilized: 343.17 MB
# Memory Efficiency: 2.23% of 15.00 GB


# Gene table in UCSC ensGenes.txt format
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.61.gff3.gz
gzip -d -k Arabidopsis_thaliana.TAIR10.61.gff3.gz
JOB=gff32ensGene_tair10_chr1
TIME="2-00:00:00"
THREADS=1
MEM=10G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python ../../../../bin/gff32ensGene.py -i Arabidopsis_thaliana.TAIR10.61.gff3"
# Submitted batch job 8793784
# Job ID: 8793784
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 01:02:58
# CPU Efficiency: 22.88% of 04:35:09 core-walltime
# Job Wall-clock time: 04:35:09
# Memory Utilized: 645.75 MB
# Memory Efficiency: 6.31% of 10.00 GB

# perl createModel.pl -m myOrg -f myOrg.fa -r RM.out -t TRF.out -g Genes.table
OUT=$WRK/pub/output/Garlic
mkdir -p $OUT
cd $OUT
GarlicDir=../../../../../bin/Garlic/bin
genome=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
RMout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna.out
TRFout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.2.5.7.80.10.50.2000.dat
ensGene=$WRK/pub/input/ensGene_table.tsv
JOB=Garlic_tair10
TIME="6-00:00:00"
THREADS=1
MEM=10G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "perl $GarlicDir/createModel.pl -m tair10 -d $OUT -f $genome -r $RMout -t $TRFout -g $ensGene"
# Submitted batch job 8795557 # failed
# Error message: Simple repeat file /powerplant/workspace/cflthc/script/TEgenomeSimulator/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.2.5.7.80.10.50.2000.dat is not in the expected format!
# Go back to check the code of Garlic to find which line would print this error message, and then figure out what the input format should be based on the code.
# Developed a script 'trfConverter4Garlic.sh' to convert the .dat file from TRF. Check ~line 20.

OUT=$WRK/pub/output/Garlic
mkdir -p $OUT
cd $OUT
GarlicDir=../../../../../bin/Garlic/bin
genome=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
RMout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna.out
TRFout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.2.5.7.80.10.50.2000.dat.reformated
ensGene=$WRK/pub/input/ensGene_table.tsv
JOB=Garlic_tair10
TIME="6-00:00:00"
THREADS=1
MEM=10G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "perl $GarlicDir/createModel.pl -m tair10 -d $OUT -f $genome -r $RMout -t $TRFout -g $ensGene"
# Submitted batch job 
# Job ID: 8796372
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:00:23
# CPU Efficiency: 92.00% of 00:00:25 core-walltime
# Job Wall-clock time: 00:00:25
# Memory Utilized: 55.54 MB
# Memory Efficiency: 0.54% of 10.00 GB

# Check the length of tair10 chr1
ml samtools
cd $WRK/pub/input/
samtools faidx $WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna 

# Create simulated sequence
OUT=$WRK/pub/output/Garlic
cd $OUT
GarlicDir=../../../../../bin/Garlic/bin
repbase2021=../../../../../../ComparativeDataSources/RepBase/RepBase26.11/RepBase26.11.embl 
size=$(awk '{print $2}' ../../input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.fai)
JOB=GarlicSim_tair10
TIME="6-00:00:00"
THREADS=1
MEM=20G
    # copy repbase embl format
mkdir -p tair10/RepBase
#cp $repbase2021 tair10/repbase/ # should not copying the .ref file, and the folder should called "RepBase" instead of "repbase"
cp -r $repbase2021 tair10/RepBase/
# have to change the name of copied "RepBase26.11.embl" to "RepeatMaskerLib.embl"
mv tair10/RepBase/RepBase26.11.embl tair10/RepBase/RepeatMaskerLib.embl
# The repbase folder should be one level up
mkdir RepBase
cp -vr tair10/RepBase/RepeatMaskerLib.embl ./RepBase/
rm -rf tair10/RepBase
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "perl $GarlicDir/createFakeSequence.pl -m tair10 -s $size -d . -o tair10_garlic_sim_chr1.fa"
# Submitted batch job 8796788 #Fatal error: cannot open tair10/tair10/tair10.GCt.W1000.data
# Submitted batch job 8797170
# Job ID: 8797170
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:00:32
# CPU Efficiency: 96.97% of 00:00:33 core-walltime
# Job Wall-clock time: 00:00:33
# Memory Utilized: 64.58 MB
# Memory Efficiency: 0.32% of 20.00 GB
# But one of the output file "*.insert" has no info

# Maybe I should convert the TE lib I used to embl format and replace the RepBase db
#tefasta=$WRK/pub/input/athrep.updated.nonredun.fasta
infa=../../../../scratch/TEGE/TEgenomeSimulator/tair10/output/TEgenomeSimulator_tair10_run6_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
cp $infa ./
outembl=$WRK/pub/input/athrep.updated.nonredun.noCenSatelli.fasta.stitched.embl
python3 - <<EOF
from Bio import SeqIO
records = []
for r in SeqIO.parse("$infa", "fasta"):
    r.annotations["molecule_type"] = "DNA"
    records.append(r)
SeqIO.write(records, "$outembl", "embl")
EOF

OUT=$WRK/pub/output/Garlic
cd $OUT
rm -rf RepBase/*
cp $outembl RepBase/
mv RepBase/athrep.updated.nonredun.noCenSatelli.fasta.stitched.embl RepBase/RepeatMaskerLib.embl
# Now rerun $GarlicDir/createFakeSequence.
GarlicDir=../../../../../bin/Garlic/bin
size=$(awk '{print $2}' ../../input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.fai)
#JOB=GarlicSim_tair10
#TIME="6-00:00:00"
#THREADS=1
#MEM=20G
perl $GarlicDir/createFakeSequence.pl -m tair10 -s $size -d . -o tair10_garlic_sim_chr1.fa
# still nothing in .inserts file

###############################################################
# Rerun createModel.pl using reformated RepeatMasker out file #
###############################################################
OUT=$WRK/pub/output/Garlic
mkdir -p $OUT
cd $OUT
GarlicDir=../../../../../bin/Garlic/bin
genome=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
RMout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna.out.reformated
TRFout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna.2.5.7.80.10.50.2000.dat.reformated
ensGene=$WRK/pub/input/ensGene_table.tsv
perl $GarlicDir/createModel.pl -m tair10 -d $OUT -f $genome -r $RMout -t $TRFout -g $ensGene

# Try to use alignment file instead # Turns out it should be .align file to be used as input instead of .out file!!!!
# Confused by this line "perl createModel.pl -m myOrg -f myOrg.fa -r RM.out -t TRF.out -g Genes.table" in Garlic's README!!!!!!!
OUT=$WRK/pub/output/Garlic
mkdir -p $OUT
cd $OUT
GarlicDir=../../../../../bin/Garlic/bin
genome=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
RMout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna.align
TRFout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna.2.5.7.80.10.50.2000.dat.reformated
ensGene=$WRK/pub/input/ensGene_table.tsv
perl $GarlicDir/createModel.pl -m tair10 -d $OUT -f $genome -r $RMout -t $TRFout -g $ensGene

# Rerun createFakeSequence.pl
# The embl file converted from A. thaliana TE lib was used as database, instead of the one downloaded from repbase.
# The converted TE lib in embl format has to be stored under the folder "RepBase" and named as "RepeatMaskerLib.embl"
OUT=$WRK/pub/output/Garlic
cd $OUT
GarlicDir=../../../../../bin/Garlic/bin
size=$(awk '{print $2}' ../../input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.fai)
perl $GarlicDir/createFakeSequence.pl -m tair10 -s $size -d . -o tair10_garlic_sim_chr1.fa -v 2> createFakeSequence.log
# Simple repeats were inserted, but TEs were not.
# Seem to be some problem in the embl file.
# example stdout with -v option:
# inserting repeat elements
# Trying to add 23% in repeats
# TE fraction = 0
# selected: ATLANTYS3:LTR/Ty3:-:29.63:2.80:18.55:495:1
# sequence for ATLANTYS3 (LTR/Ty3) not found!
#  - evolveRepeat returned bad status
# TE fraction = 0
# selected: ATLINE1A:LINE/unknown:-:12.61:0.00:0.00:110:1
# sequence for ATLINE1A (LINE/unknown) not found!
#  - evolveRepeat returned bad status
# TE fraction = 0
# selected: ATHILA5:LTR/Ty3:+:21.80:0.71:4.03:545:1
# sequence for ATHILA5 (LTR/Ty3) not found!
#  - evolveRepeat returned bad status
# TE fraction = 0
# selected: ATHILA7:LTR/Ty3:+:25.77:8.65:4.67:299:1
# sequence for ATHILA7 (LTR/Ty3) not found!
#  - evolveRepeat returned bad status

# Try if replacing '#' to ' ' works
sed 's/#/ /g' ./RepBase/RepeatMaskerLib.embl > ./RepBase/temp.embl
mv ./RepBase/temp.embl ./RepBase/RepeatMaskerLib.embl

# Rerun createFakeSequence.pl
# The embl file converted from A. thaliana TE lib was used as database, instead of the one downloaded from repbase.
# The converted TE lib in embl format has to be stored under the folder "RepBase" and named as "RepeatMaskerLib.embl"
JOB=GarlicSim_tair10
TIME="6-00:00:00"
THREADS=1
MEM=20G
OUT=$WRK/pub/output/Garlic
cd $OUT
GarlicDir=../../../../../bin/Garlic/bin
size=$(awk '{print $2}' ../../input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.fai) # echo $size # 30427671
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "perl $GarlicDir/createFakeSequence.pl -m tair10 -s $size -d . -o tair10_garlic_sim_chr1.fa -v 2> createFakeSequence.log"
# Submitted batch job 8800141 # aklppg32
# Job ID: 8800141
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:01:07
# CPU Efficiency: 98.53% of 00:01:08 core-walltime
# Job Wall-clock time: 00:01:08
# Memory Utilized: 172.16 MB
# Memory Efficiency: 0.84% of 20.00 GB



###########################
# Use size of non-TE chr1 #
###########################
# It seems that the specified size in input refers to the lenth of the non-repeat sequences to be synthesized.
# Beter try the size of nonTE chr1
# Ended up the previous one that use size of the whole chr1 is correct
ml samtools
cd $WRK/pub/input
cp ../../../../scratch/TEGE/TEgenomeSimulator/tair10/output/TEgenomeSimulator_tair10_run7_result/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed ./
samtools faidx GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed

# Rerun createModel.pl so I can save output to a different folder
OUT=$WRK/pub/output/Garlic
mkdir -p $OUT
cd $OUT
GarlicDir=../../../../../bin/Garlic/bin
genome=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
RMout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna.align
TRFout=$WRK/pub/input/GCF_000001735.4_TAIR10.1_genomic.fna.2.5.7.80.10.50.2000.dat.reformated
ensGene=$WRK/pub/input/ensGene_table.tsv
perl $GarlicDir/createModel.pl -m tair10_2 -d $OUT -f $genome -r $RMout -t $TRFout -g $ensGene

# Rerun createFakeSequence.pl using a different synthesized size
# The embl file converted from A. thaliana TE lib was used as database, instead of the one downloaded from repbase.
# The converted TE lib in embl format has to be stored under the folder "RepBase" and named as "RepeatMaskerLib.embl"
JOB=GarlicSim_tair10_2
TIME="6-00:00:00"
THREADS=1
MEM=20G
OUT=$WRK/pub/output/Garlic
cd $OUT
GarlicDir=../../../../../bin/Garlic/bin
size=$(grep "NC_003070.9" ../../input/GCF_000001735.4_TAIR10.1_genomic.fna.masked.reformatted.nonTE.emptfixed.fai | awk '{print $2}' ) # echo $size # 27195677
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "perl $GarlicDir/createFakeSequence.pl -m tair10_2 -s $size -d . -o tair10_2_garlic_sim_chr1.fa -v 2> createFakeSequence2.log"
# Submitted batch job 8800142 # aklppb43
# Job ID: 8800142
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:22:31
# CPU Efficiency: 99.48% of 00:22:38 core-walltime
# Job Wall-clock time: 00:22:38
# Memory Utilized: 3.28 GB
# Memory Efficiency: 16.42% of 20.00 GB


ml samtools
samtools faidx tair10_garlic_sim_chr1.fa.fasta
samtools faidx tair10_2_garlic_sim_chr1.fa.fasta


# The records of TE and Simple Repeat insertions in the .inserts file are in different format. Need to saperate them before importing them to R.
# Note: Garlic analyses CDS profile but doesn't simulate CDS in the synthetic genome.
WRK=$(pwd)
OUT=$WRK/pub/output/Garlic
cd $OUT
grep -v "#" tair10_garlic_sim_chr1.fa.inserts | grep -v "SIMPLE" > tair10_garlic_sim_chr1.fa.inserts.te
grep "SIMPLE" tair10_garlic_sim_chr1.fa.inserts > tair10_garlic_sim_chr1.fa.inserts.simple
grep -v "#" tair10_2_garlic_sim_chr1.fa.inserts | grep -v "SIMPLE" > tair10_2_garlic_sim_chr1.fa.inserts.te
grep "SIMPLE" tair10_2_garlic_sim_chr1.fa.inserts > tair10_2_garlic_sim_chr1.fa.inserts.simple

# Convert .insert output to gff format
cd $OUT
../../script/garlic2tegff.sh tair10_garlic_sim_chr1.fa.inserts.te > tair10_garlic_sim_chr1.fa.inserts.te.gff
../../script/garlic2simprepgff.sh tair10_garlic_sim_chr1.fa.inserts.simple > tair10_garlic_sim_chr1.fa.inserts.simple.gff
../../script/garlic2tegff.sh tair10_2_garlic_sim_chr1.fa.inserts.te > tair10_2_garlic_sim_chr1.fa.inserts.te.gff
../../script/garlic2simprepgff.sh tair10_2_garlic_sim_chr1.fa.inserts.simple > tair10_2_garlic_sim_chr1.fa.inserts.simple.gff

# The .inserts output looks like:
### ARTIFICIAL SEQUENCE 1 ###
#(INI   END]-zero_based NUM     REPEAT  REPEAT_EVOL
#20462816        20463231        1       VANDAL3:DNA/MuDR:+:21.88:4.65:3.88:411:1;1:9.00:10.22:4.62:3.65[$inserted_TE_sequence] # Note the coordinates are zero-based, so length = END - INI. 
# The followings are information from the script `createFakeSequence.pl` (from line 1156)
#push @inserts,
      #    "$pos\t$pos_end\t$urep\t$new\[$seq\]\t$info\n$con\n$mat\n$mut\n";
      # RMH: The transitions/transversion/insertions/deletions are now calcualted
      #      by evolveRepeat and reported after the canonical repeat pattern.
      #      The format is:
      #         prototype_details;instance_details[sequence],...
      #      Where prototype details are:
      #         family:class:orientation:%div:%ins:%del:fraglen:break
      #
      #         family: The family identifier for the TE
      #         class: The TE classification
      #         orientation: '+','-'
      #         %div: The percent divergence of the prototype being simulated
      #         %ins: The percent insertion of the prorotype being simulated
      #         %del: The percent deletion of the prototype being simulated
      #         fraglen: The length of the TE fragment being simulated           # Note: this is the length before adding small insertion and deletion!
      #         break: The number of fragments for this prototype
      #
      #      and instance details are:
      #         level:%transitions:%transversions:%insertions:%deletions
      #
      #         level: For fragmented insertions this indicates the level of 
      #                the fragment (starting from 1).  e.g An insertion of a
      #                MLT1 in an AluSx would have three fragments listed as:
      #                AluSx:::::::;1::::
      #                MLT11:::::::;2::::
      #                AluSx:::::::;1::::
      #
      #         %transitions: The percent transitions (relative to the consensus
      #                       fragment length).
      #         %transversions: The percent transversions (relative to the consensus
      #                       fragment length).
      #         %deletions: The percent deletions (relative to the consensus
      #                       fragment length).
      #         %insertions: The percent insertions (relative to the consensus
      #                       fragment length).
      # e.g:
      # AluSz:SINE/Alu:-:16.4:0.0:0.6:308:1;1:7.14:8.12:0.00:0.32[GGCCGGGGGCGG...]
      #    
# 