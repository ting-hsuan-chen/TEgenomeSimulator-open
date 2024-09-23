#!/bin/bash
#Ting-Hsuan Chen
#2024-04-02
#Show kmer profile

## 21 mer
module load jellyfish/2.2.10
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/kmer_analysis
genomename="tair10_sim"
genome=$WRK/10.01_TEgenomeSimulator_output/sim_tair10_genome/tair10_genome_sequence_out_nest.fasta

JOB="jellyfish_bc"
RAM=16G
THREADS=20
TIME="06:00:00"
mkdir -p $OUT/tair10_sim
cd $OUT/tair10_sim
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish bc -m 21 -s 100G -t 20 -o $genomename.bc $genome"
#Submitted batch job 5185070
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:08:52
CPU Efficiency: 5.33% of 02:46:20 core-walltime
Job Wall-clock time: 00:08:19
Memory Utilized: 260.77 GB
Memory Efficiency: 1629.82% of 16.00 GB

JOB="jellyfish_count"
RAM=280G
THREADS=16
TIME="06:00:00"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish count -m 21 -s 100G -t 16 --bc $genomename.bc $genome"
Submitted batch job 5185183
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 16
CPU Utilized: 00:33:51
CPU Efficiency: 34.68% of 01:37:36 core-walltime
Job Wall-clock time: 00:06:06
Memory Utilized: 588.46 GB
Memory Efficiency: 210.16% of 280.00 GB

JOB="jellyfish_hist"
RAM=1G
THREADS=4
TIME="06:00:00"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish histo -t 16 mer_counts.jf > kmer_histo.out"
Submitted batch job 5185225
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 4
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 00:00:16 core-walltime
Job Wall-clock time: 00:00:04
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 1.00 GB (1.00 GB/node)

## 16 mer
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/kmer_analysis
genomename="tair10_sim"
genome=$WRK/10.01_TEgenomeSimulator_output/sim_tair10_genome/tair10_genome_sequence_out_nest.fasta
kmer=16
JOB="jellyfish_"$kmer"_mer"
RAM=250G
THREADS=20
TIME="06:00:00"
mkdir -p $OUT/$genomename
cd $OUT/$genomename

sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ${JOB}.out
#SBATCH -e ${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

jellyfish bc -m $kmer -s 100G -t 20 -o $genomename'_'$kmer'm'.bc $genome
jellyfish count -m $kmer -s 100G -t 16 --bc $genomename'_'$kmer'm'.bc -o $genomename'_'$kmer'm'.jf $genome
jellyfish histo -t 16 $genomename'_'$kmer'm'.jf > $genomename'_histo_'$kmer'm'.out

EOF
Submitted batch job 5197628
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:20:22
CPU Efficiency: 4.18% of 08:06:40 core-walltime
Job Wall-clock time: 00:24:20
Memory Utilized: 268.39 GB
Memory Efficiency: 107.36% of 250.00 GB

### The original tair10 genome
## 21 mer
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/kmer_analysis
genomename="tair10_ori"
genome=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_ref/AthaGenome.fa

JOB="jellyfish_bc"
RAM=16G
THREADS=20
TIME="06:00:00"
mkdir -p $OUT/$genomename
cd $OUT/$genomename
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish bc -m 21 -s 100G -t 20 -o $genomename.bc $genome"
#Submitted batch job 5186330
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:08:20
CPU Efficiency: 2.92% of 04:45:20 core-walltime
Job Wall-clock time: 00:14:16
Memory Utilized: 260.77 GB
Memory Efficiency: 1629.82% of 16.00 GB


JOB="jellyfish_count"
RAM=280G
THREADS=16
TIME="06:00:00"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish count -m 21 -s 100G -t 16 --bc $genomename.bc $genome"
Submitted batch job 5190166
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 16
CPU Utilized: 00:30:07
CPU Efficiency: 6.48% of 07:44:48 core-walltime
Job Wall-clock time: 00:29:03
Memory Utilized: 588.46 GB
Memory Efficiency: 210.16% of 280.00 GB


JOB="jellyfish_hist"
RAM=1G
THREADS=4
TIME="06:00:00"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish histo -t 16 mer_counts.jf > kmer_histo.out"
Submitted batch job 5190170
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 4
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 00:00:12 core-walltime
Job Wall-clock time: 00:00:03
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 1.00 GB (1.00 GB/node)



## 16 mer
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/kmer_analysis
genomename="tair10_ori"
genome=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_ref/AthaGenome.fa

JOB="jellyfish_bc"
RAM=16G
THREADS=20
TIME="06:00:00"
mkdir -p $OUT/$genomename
cd $OUT/$genomename
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish bc -m 16 -s 100G -t 20 -o $genomename'_16m'.bc $genome"
#Submitted batch job 5197625
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:09:44
CPU Efficiency: 7.28% of 02:13:40 core-walltime
Job Wall-clock time: 00:06:41
Memory Utilized: 260.77 GB
Memory Efficiency: 1629.82% of 16.00 GB


JOB="jellyfish_count"
RAM=280G
THREADS=16
TIME="06:00:00"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish count -m 16 -s 100G -t 16 --bc $genomename'_16m'.bc $genome"
Submitted batch job 5197626
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 16
CPU Utilized: 00:06:54
CPU Efficiency: 13.99% of 00:49:20 core-walltime
Job Wall-clock time: 00:03:05
Memory Utilized: 268.39 GB
Memory Efficiency: 95.85% of 280.00 GB

JOB="jellyfish_hist"
RAM=1G
THREADS=4
TIME="06:00:00"
sbatch --cpus-per-task=$THREADS --mem $RAM -t $TIME -J $JOB -o $JOB.out -e $JOB.err --wrap "jellyfish histo -t 16 mer_counts.jf > $genomename'_histo_16m'.out"

Submitted batch job 5197627
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 4
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 00:00:04 core-walltime
Job Wall-clock time: 00:00:01
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 1.00 GB (1.00 GB/node)



## 17 mer
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/kmer_analysis
genomename="tair10_ori"
genome=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_ref/AthaGenome.fa
kmer=17
JOB="jellyfish_"$kmer"_mer"
RAM=250G
THREADS=20
TIME="06:00:00"
mkdir -p $OUT/$genomename
cd $OUT/$genomename

sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ${JOB}.out
#SBATCH -e ${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

jellyfish bc -m $kmer -s 100G -t 20 -o $genomename'_'$kmer'm'.bc $genome
jellyfish count -m $kmer -s 100G -t 16 --bc $genomename'_'$kmer'm'.bc -o $genomename'_'$kmer'm'.jf $genome
jellyfish histo -t 16 $genomename'_'$kmer'm'.jf > $genomename'_histo_'$kmer'm'.out

EOF
Submitted batch job 5197631
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:19:40
CPU Efficiency: 6.18% of 05:18:00 core-walltime
Job Wall-clock time: 00:15:54
Memory Utilized: 291.25 GB
Memory Efficiency: 116.50% of 250.00 GB

## 22 mer
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/kmer_analysis
genomename="tair10_ori"
genome=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_ref/AthaGenome.fa
kmer=22
JOB="jellyfish_"$kmer"_mer"
RAM=250G
THREADS=20
TIME="06:00:00"
mkdir -p $OUT/$genomename
cd $OUT/$genomename

sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ${JOB}.out
#SBATCH -e ${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

jellyfish bc -m $kmer -s 100G -t 20 -o $genomename'_'$kmer'm'.bc $genome
jellyfish count -m $kmer -s 100G -t 16 --bc $genomename'_'$kmer'm'.bc -o $genomename'_'$kmer'm'.jf $genome
jellyfish histo -t 16 $genomename'_'$kmer'm'.jf > $genomename'_histo_'$kmer'm'.out

EOF
Submitted batch job 5197632
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:47:20
CPU Efficiency: 13.08% of 06:02:00 core-walltime
Job Wall-clock time: 00:18:06
Memory Utilized: 619.17 GB
Memory Efficiency: 247.67% of 250.00 GB

## 23 to 28 mer
for i in {23..28}
do
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/kmer_analysis
genomename="tair10_ori"
genome=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_ref/AthaGenome.fa
kmer=$i
JOB="jellyfish_"$kmer"_mer"
RAM=250G
THREADS=20
TIME="06:00:00"
mkdir -p $OUT/$genomename
cd $OUT/$genomename

sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ${JOB}.out
#SBATCH -e ${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

jellyfish bc -m $kmer -s 100G -t 20 -o $genomename'_'$kmer'm'.bc $genome
jellyfish count -m $kmer -s 100G -t 16 --bc $genomename'_'$kmer'm'.bc -o $genomename'_'$kmer'm'.jf $genome
jellyfish histo -t 16 $genomename'_'$kmer'm'.jf > $genomename'_histo_'$kmer'm'.out

EOF
done
Submitted batch job 5197642 #failed 23mer
Submitted batch job 5197643 #completed 24mer
Submitted batch job 5197644 #completed 25mer
Submitted batch job 5197645 #failed 26mer
Submitted batch job 5197646 #failed 27mer
Submitted batch job 5197647 #failed 28mer


## 31 mer
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/kmer_analysis
genomename="tair10_ori"
genome=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_ref/AthaGenome.fa
kmer=31
JOB="jellyfish_"$kmer"_mer"
RAM=500G
THREADS=20
TIME="06:00:00"
mkdir -p $OUT/$genomename
cd $OUT/$genomename

sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ${JOB}.out
#SBATCH -e ${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

jellyfish bc -m $kmer -s 100G -t 20 -o $genomename'_'$kmer'm'.bc $genome
jellyfish count -m $kmer -s 100G -t 16 --bc $genomename'_'$kmer'm'.bc -o $genomename'_'$kmer'm'.jf $genome
jellyfish histo -t 16 $genomename'_'$kmer'm'.jf > $genomename'_histo_'$kmer'm'.out

EOF
Submitted batch job 5197649
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:45:44
CPU Efficiency: 16.92% of 04:30:20 core-walltime
Job Wall-clock time: 00:13:31
Memory Utilized: 919.07 GB
Memory Efficiency: 183.81% of 500.00 GB

## 100 mer
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/10.01_TEgenomeSimulator_output/kmer_analysis
genomename="tair10_ori"
genome=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_ref/AthaGenome.fa
kmer=100
JOB="jellyfish_"$kmer"_mer"
RAM=500G
THREADS=20
TIME="06:00:00"
mkdir -p $OUT/$genomename
cd $OUT/$genomename

sbatch << EOF
#!/bin/bash
#SBATCH -J ${JOB}
#SBATCH -o ${JOB}.out
#SBATCH -e ${JOB}.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$RAM
#SBATCH --time=$TIME

jellyfish bc -m $kmer -s 100G -t 20 -o $genomename'_'$kmer'm'.bc $genome
jellyfish count -m $kmer -s 100G -t 16 --bc $genomename'_'$kmer'm'.bc -o $genomename'_'$kmer'm'.jf $genome
jellyfish histo -t 16 $genomename'_'$kmer'm'.jf > $genomename'_histo_'$kmer'm'.out

EOF
Submitted batch job 5197660