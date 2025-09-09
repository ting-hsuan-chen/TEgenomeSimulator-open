#!/bin/bash
# Navigate to the cloned repository of TEgenomeSimulator
ml conda
conda create --name cflthc_tesimeval python=3.9 
conda activate cflthc_tesimeval
pip install biopython numpy scipy scikit-learn pandas matplotlib

WRK=$(pwd)
IN=$WRK/pub/input
OUT=$WRK/pub/output

############################################
# Run repeatMasker for Recovery assessment #
############################################
ml singularity/3
TOOL=/workspace/cflthc/bin/TETools
mkdir -p $OUT/TE_sim_eval
cd $OUT/TE_sim_eval
THREADS=10
RAM=15G
TIME="06:00:00"
GRP=1

# TEGSm2
prefix=tair10_ch1_m2
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
telib=$IN/athrep.updated.nonredun.noCenSatelli.fasta
JOB=RMask_$prefix
mkdir -p $OUT/TE_sim_eval/${prefix}/RepeatMasker
cd $OUT/TE_sim_eval/${prefix}/RepeatMasker
ln -sf $genome sim_genome.fa
ln -sf $telib sim_used_telib.fa 
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -x -q -no_is -norna -nolow -div 40 -lib sim_used_telib.fa -cutoff 225 sim_genome.fa"
# Submitted batch job 8818970
# Job ID: 8818970
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 10
# CPU Utilized: 00:03:12
# CPU Efficiency: 19.20% of 00:16:40 core-walltime
# Job Wall-clock time: 00:01:40
# Memory Utilized: 497.64 MB
# Memory Efficiency: 3.24% of 15.00 GB

# TEGSm0
prefix=tair10_ch1_m0
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
telib=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
JOB=RMask_$prefix
mkdir -p $OUT/TE_sim_eval/${prefix}/RepeatMasker
cd $OUT/TE_sim_eval/${prefix}/RepeatMasker
ln -sf $genome sim_genome.fa
ln -sf $telib sim_used_telib.fa 
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -x -q -no_is -norna -nolow -div 40 -lib sim_used_telib.fa -cutoff 225 sim_genome.fa"
# Submitted batch job 8819110
# Job ID: 8819110
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:02:23
# CPU Efficiency: 26.68% of 00:08:56 core-walltime
# Job Wall-clock time: 00:08:56
# Memory Utilized: 453.97 MB
# Memory Efficiency: 2.96% of 15.00 GB

# TEGSm0_full
prefix=tair10_ch1_m0_full
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
telib=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
JOB=RMask_$prefix
mkdir -p $OUT/TE_sim_eval/${prefix}/RepeatMasker
cd $OUT/TE_sim_eval/${prefix}/RepeatMasker
ln -sf $genome sim_genome.fa
ln -sf $telib sim_used_telib.fa 
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -x -q -no_is -norna -nolow -div 40 -lib sim_used_telib.fa -cutoff 225 sim_genome.fa"
# Submitted batch job 8819112
# Job ID: 8819112
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:34:44
# CPU Efficiency: 45.14% of 01:16:57 core-walltime
# Job Wall-clock time: 01:16:57
# Memory Utilized: 2.61 GB
# Memory Efficiency: 17.40% of 15.00 GB

# Garlic
prefix=Garlic
genome=$WRK/pub/output/${prefix}/tair10_garlic_sim_chr1.fa.fasta
telib=$OUT/TEgenomeSimulator_tair10_ch1_m2_result/athrep.updated.nonredun.noCenSatelli.fasta.stitched
JOB=RMask_$prefix
mkdir -p $OUT/TE_sim_eval/${prefix}/RepeatMasker
cd $OUT/TE_sim_eval/${prefix}/RepeatMasker
ln -sf $genome sim_genome.fa
ln -sf $telib sim_used_telib.fa 
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$TOOL/dfam-tetools.sif RepeatMasker -pa $THREADS -x -q -no_is -norna -nolow -div 40 -lib sim_used_telib.fa -cutoff 225 sim_genome.fa"
# Submitted batch job 8819231
# Job ID: 8819231
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 10
# CPU Utilized: 00:03:52
# CPU Efficiency: 18.56% of 00:20:50 core-walltime
# Job Wall-clock time: 00:02:05
# Memory Utilized: 418.07 MB
# Memory Efficiency: 2.72% of 15.00 GB

################################################
# Run te_sim_eval.py & te_sim_eval_bedtools.sh #
################################################

#### TEGSm2 ####
prefix=tair10_ch1_m2
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
tegff=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_repeat_annotation_out_final.gff
rmout=$OUT/TE_sim_eval/${prefix}/RepeatMasker/sim_genome.fa.out

THREADS=1
RAM=30G
TIME="6-00:00:00"

# bedtools intersect
ml bedtools/2.30.0
out=$OUT/TE_sim_eval/${prefix}/bedtools_cut08
mkdir -p $out
cd $out
# THREADS=1
# RAM=10G
# TIME="06:00:00"
# JOB=tesimeval_${prefix}_bedtools
cut=0.8
$WRK/pub/script/te_sim_eval_bedtools.sh $rmout $tegff $cut
# sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$WRK/pub/script/te_sim_eval_bedtools.sh $rmout $tegff $cut"
# Submitted batch job 8823057


# kmer
out=$OUT/TE_sim_eval/${prefix}/kmer
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_kmer
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --rm_out $rmout --outdir $out --k 6 --win 1000 --eval kmer"
# Submitted batch job 8819177
# Job ID: 8819177
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:00:54
# CPU Efficiency: 94.74% of 00:00:57 core-walltime
# Job Wall-clock time: 00:00:57
# Memory Utilized: 284.32 MB
# Memory Efficiency: 0.93% of 30.00 GB

# entropy
out=$OUT/TE_sim_eval/${prefix}/entropy
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_entropy
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --rm_out $rmout --outdir $out --k 6 --win 1000 --eval entropy"
# Submitted batch job 8819182
# Job ID: 8819182
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:00:14
# CPU Efficiency: 87.50% of 00:00:16 core-walltime
# Job Wall-clock time: 00:00:16
# Memory Utilized: 34.21 MB
# Memory Efficiency: 0.11% of 30.00 GB
# compress
out=$OUT/TE_sim_eval/${prefix}/compress
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_compress
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --rm_out $rmout --outdir $out --k 6 --win 1000 --eval compress"
# Submitted batch job 8819183
# Job ID: 8819183
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:01:06
# CPU Efficiency: 94.29% of 00:01:10 core-walltime
# Job Wall-clock time: 00:01:10
# Memory Utilized: 425.36 MB
# Memory Efficiency: 1.38% of 30.00 GB
# repeatmasker
out=$OUT/TE_sim_eval/${prefix}/rm
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_rm
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --rm_out $rmout --outdir $out --k 6 --win 1000 --eval rm"
# Submitted batch job 8819184
# Job ID: 8819184
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:40:35
# CPU Efficiency: 99.47% of 00:40:48 core-walltime
# Job Wall-clock time: 00:40:48
# Memory Utilized: 262.79 MB
# Memory Efficiency: 0.86% of 30.00 GB

#### TEGSm0 ####
prefix=tair10_ch1_m0
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
tegff=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_repeat_annotation_out_final.gff
rmout=$OUT/TE_sim_eval/${prefix}/RepeatMasker/sim_genome.fa.out

THREADS=1
RAM=30G
TIME="6-00:00:00"

for eval in kmer entropy compress rm
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --rm_out $rmout --outdir $out --k 6 --win 1000 --eval ${eval}"
done
# Submitted batch job 8819214 # completed
# Submitted batch job 8819215 # completed
# Submitted batch job 8819216 # completed
# Submitted batch job 8819217
# Job ID: 8819217
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:44:53
# CPU Efficiency: 98.57% of 00:45:32 core-walltime
# Job Wall-clock time: 00:45:32
# Memory Utilized: 509.36 MB
# Memory Efficiency: 1.66% of 30.00 GB

# bedtools intersect
ml bedtools/2.30.0
out=$OUT/TE_sim_eval/${prefix}/bedtools_cut08
mkdir -p $out
cd $out
#THREADS=1
#RAM=10G
#TIME="06:00:00"
#JOB=tesimeval_${prefix}_bedtools
cut=0.8
$WRK/pub/script/te_sim_eval_bedtools.sh $rmout $tegff $cut
# sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$WRK/pub/script/te_sim_eval_bedtools.sh $rmout $tegff $cut"
# Submitted batch job 8822988

#### TEGSm0_full ####
prefix=tair10_ch1_m0_full
genome=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_genome_sequence_out_final.fasta
tegff=$OUT/TEgenomeSimulator_${prefix}_result/${prefix}_repeat_annotation_out_final.gff
rmout=$OUT/TE_sim_eval/${prefix}/RepeatMasker/sim_genome.fa.out

THREADS=1
RAM=30G
TIME="6-00:00:00"

for eval in kmer entropy compress rm
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --rm_out $rmout --outdir $out --k 6 --win 1000 --eval ${eval}"
done
# Submitted batch job 8819237 # completed
# Submitted batch job 8819238 # completed
# Submitted batch job 8819239 # completed
# Submitted batch job 8819240 # cancelled

# bedtools intersect
# ml bedtools/2.30.0
# out=$OUT/TE_sim_eval/${prefix}/bedtools
# mkdir -p $out
# cd $out
# THREADS=1
# RAM=10G
# TIME="06:00:00"
# JOB=tesimeval_${prefix}_bedtools
# sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "$WRK/pub/script/te_sim_eval_bedtools.sh $rmout $tegff"
# Submitted batch job 8822991


#### Garlic ####
prefix=Garlic
genome=$WRK/pub/output/${prefix}/tair10_garlic_sim_chr1.fa.fasta
tegff=$WRK/pub/output/${prefix}/tair10_garlic_sim_chr1.fa.inserts.te.gff
rmout=$OUT/TE_sim_eval/${prefix}/RepeatMasker/sim_genome.fa.out

THREADS=1
RAM=30G
TIME="6-00:00:00"

for eval in kmer entropy compress rm
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --rm_out $rmout --outdir $out --k 6 --win 1000 --eval ${eval}"
done
# Submitted batch job 8819492 #completed
# Submitted batch job 8819493 #completed
# Submitted batch job 8819494 #completed
# Submitted batch job 8819495
# Job ID: 8819495
# Cluster: powerplant
# User/Group: cflthc/cflthc
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 09:32:44
# CPU Efficiency: 99.38% of 09:36:18 core-walltime
# Job Wall-clock time: 09:36:18
# Memory Utilized: 434.60 MB
# Memory Efficiency: 1.41% of 30.00 GB

# bedtools intersect
ml bedtools/2.30.0
out=$OUT/TE_sim_eval/${prefix}/bedtools_cut08
mkdir -p $out
cd $out
sed 's/^chr1/artificial_sequence_1/g' $tegff > truth_tegff_chrmodified.gff
cut=0.8
$WRK/pub/script/te_sim_eval_bedtools.sh $rmout truth_tegff_chrmodified.gff $cut

#### Original ####
prefix=original
genome=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
tegff=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.out.gff

THREADS=1
RAM=30G
TIME="6-00:00:00"

for eval in kmer entropy compress 
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --rm_out $rmout --outdir $out --k 6 --win 1000 --eval ${eval}"
done
# Submitted batch job 8819623 # completed
# Submitted batch job 8819624 # completed
# Submitted batch job 8819625 # completed

###### no run ######
prefix=original
genome=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna
tegff=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.out.gff
rmout=$IN/GCF_000001735.4_TAIR10.1_genomic_chr1.fna.out

THREADS=1
RAM=30G
TIME="6-00:00:00"

for eval in rm 
do
out=$OUT/TE_sim_eval/${prefix}/${eval}
mkdir -p $out
cd $out
JOB=tesimeval_${prefix}_${eval}
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "python $WRK/pub/script/te_sim_eval.py --genome $genome --te_gff $tegff --rm_out $rmout --outdir $out --k 6 --win 1000 --eval ${eval}"
done
