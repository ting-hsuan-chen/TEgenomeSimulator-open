#!/bin/bash
############
## mode 2 ##
############
# Navigate to the cloned repository of TEgenomeSimulator
ml apptainer
WRK=$(pwd)
TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/pub/input
OUT=$WRK/pub/output
mkdir -p $IN $OUT

prefix=tair10_ch1_m2
genome=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta

JOB=TEGS_tair10_ch1_m2
TIME="06:00:00"
THREADS=10
MEM=5G
cd $OUT
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "apptainer run $TEgenomeSimulator -M 2 -p $prefix -g $genome -r $repeat -r2 $repeat -o $OUT -t $THREADS"
# Submitted batch job 8806406
# Job ID: 8806406
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 10
# CPU Utilized: 00:04:31
# CPU Efficiency: 18.44% of 00:24:30 core-walltime
# Job Wall-clock time: 00:02:27
# Memory Utilized: 572.30 MB
# Memory Efficiency: 11.18% of 5.00 G

############
## mode 0 ##
############
# Navigate to the cloned repository of TEgenomeSimulator
WRK=$(pwd)
# Use user-specified TE table
ml conda
#TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
IN=$WRK/pub/input
OUT=$WRK/pub/output
prefix=tair10_ch1_m0
repeat=$IN/athrep.updated.nonredun.noCenSatelli.fasta
tetable=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/TElib_sim_list_mode2.table
# Create chr index file
#cd $IN
ml samtools
conda activate hraaxt_emboss # use emboss' tool infoseq to extract gc%
cd $OUT/TEgenomeSimulator_tair10_ch1_m2_result 
samtools faidx GCF_000001735.4_TAIR10.1_genomic_chr1.fna.masked.reformatted.nonTE.emptfixed
nontesize=$(awk '{print $2}' GCF_000001735.4_TAIR10.1_genomic_chr1.fna.masked.reformatted.nonTE.emptfixed.fai)
gc=$(cat GCF_000001735.4_TAIR10.1_genomic_chr1.fna.masked.reformatted.nonTE.emptfixed | infoseq -auto -only -name -length -pgc stdin | awk 'NR>1{print $3}')
echo "chr1,$nontesize,$gc" > $IN/tair10_chr1_index_nonte.csv
conda deactivate

cd $OUT
conda activate TEgenomeSimulator
TEgenomeSimulator=$WRK/TEgenomeSimulator/TEgenomeSimulator.py
indcsv=$IN/tair10_chr1_index_nonte.csv
repeat=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
tetable=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/TElib_sim_list_mode2.table
JOB=TEGS_tair10_ch1_m0
TIME="06:00:00"
THREADS=1
MEM=15G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 0 -p $prefix -c $indcsv -r $repeat -o $OUT --te_table $tetable --frag_mode 2"     
# Submitted batch job 8818984 
# Job ID: 8818984
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:00:26
# CPU Efficiency: 55.32% of 00:00:47 core-walltime
# Job Wall-clock time: 00:00:47
# Memory Utilized: 177.46 MB
# Memory Efficiency: 1.16% of 15.00 GB

##########################
## mode 0 (full mode 0) ##
##########################
# Navigate to the cloned repository of TEgenomeSimulator
WRK=$(pwd)
# Use user-specified TE table
ml conda
#TEgenomeSimulator=$WRK/TEgenomeSimulator_v1.0.0.sif
conda activate TEgenomeSimulator
IN=$WRK/pub/input
OUT=$WRK/pub/output
cd $OUT
prefix=tair10_ch1_m0_full
TEgenomeSimulator=$WRK/TEgenomeSimulator/TEgenomeSimulator.py
indcsv=$IN/tair10_chr1_index_nonte.csv # Use chr index file created in previous code block
repeat=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
tetable=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/TElib_sim_list_mode2.table
maxcp=$(awk -F'\t' 'NR > 1 {if ($4 > max) max=$4} END {print max}' $tetable) #429
mincp=$(awk -F'\t' 'NR>1 {min=$4} NR>2 {if ($4 < min) min=$4} END {print min}' $tetable) #5
maxidn=$(awk -F'\t' 'NR > 1 {if ($5 > max) max=$5} END {print max}' $tetable) #100
minidn=$(awk -F'\t' 'NR>1 {min=$5} NR>2 {if ($5 < min) min=$5} END {print min}' $tetable) #79
maxsd=$(awk -F'\t' 'NR > 1 {if ($6 > max) max=$6} END {print max}' $tetable) #30
minsd=$(awk -F'\t' 'NR>1 {min=$6} NR>2 {if ($6 < min) min=$6} END {print min}' $tetable) #8
intfreq=$(awk -F'\t' 'NR>1 {
    gsub(/[\[\]]/, "", $NF)
    n = split($NF, arr, ",")
    for (i=1; i<=n; i++) {
        gsub(/^ +| +$/, "", arr[i])
        total++
        if (arr[i] == 1 || arr[i] == 1.0) count++
    }
}
END {
    if (total > 0) {
        print count / total
    } else {
        print 0
    }
}' $tetable)

JOB=TEGS_tair10_ch1_m0_full
TIME="06:00:00"
THREADS=1
MEM=15G
sbatch --cpus-per-task=$THREADS --mem $MEM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $TEgenomeSimulator -M 0 -p $prefix -c $indcsv -r $repeat -m $maxcp -n $mincp --maxidn $maxidn --minidn $minidn --maxsd $maxsd --minsd $minsd -a 0.5 -b 0.6 -i $intfreq -o $OUT --frag_mode 0"     
# Submitted batch job 8818985 #aklppb40
# Job ID: 8818985
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:25:41
# CPU Efficiency: 99.04% of 00:25:56 core-walltime
# Job Wall-clock time: 00:25:56
# Memory Utilized: 1.10 GB
# Memory Efficiency: 7.35% of 15.00 GB