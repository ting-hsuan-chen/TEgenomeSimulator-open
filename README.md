# TEgenomeSimulator
A tool to simulate TE mutation and insertion into a random-synthesised or user-provided genome.

## Introduction
TEgenomeSimulator was created based on [Matias Rodriguez & Wojciech Makałowski. Software evaluation for de novo detection of transposons. 2022 Mobil DNA](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-022-00266-2). The original scripts were from [denovoTE-eval](https://github.com/IOB-Muenster/denovoTE-eval). Some existing functions had been modified and several new functions had been created for TEgenomeSimulator.

### TEgenomeSimulator flowchart
<img src="TEgenomeSimulator.png" width=720>

### Two modes
TEgenomeSimulator has two modes to be specified by user using the argument `-M` or `--mode`:
1) `-M 0` or`--mode 0`: **Random Synthesized Genome Mode**. This mode would synthesize a random genome with multiple chromosomes that are later inserted with TEs randomly. With this mode, user must use `-c` or `--chridx` to provide a chromosome index text file in the formate like the file [rabdin_genome_chr_index.txt](./test/input/random_genome_chr_index.txt), which contains 3 columns seperated by commas, representing the desired chromosome id, chromosme length, and GC content.
2) `-M 1` or`--mode 1`: **Custom Genome Mode**. It utilises a user-provided genome containing multiple chromosomes where TE bases had been removed, providing a customised genome canvas for random TE insertion. When this mode is enabled, the user must use `-g` or `--genome` to provide a genome fasta file, e.g. [Dondhond.chromosomes.only.fa.nonTE.2chrs](./test/input/Donghong.chromosomes.only.fa.nonTE.2chrs.gz), which contains TE-depleted chromosome 1 and 2 from the published [_Actinidia chinensis_ var. 'Donghong'](https://doi.org/10.1016/j.molp.2022.12.022). In this fasta file, the TE sequences on the original chromosome 1 and 2 have been detected and removed. 

### Other required input files/information
- `-p` or `--prefix`: a prefix for your simulation
- `-r` or `--repeat`: a fasta file containing TE family sequences (i.e. a TE library) to be mutated and then inserted into genome. The TE lib in our test [combined_curated_TE_lib_ATOSZM_selected.fasta](./test/input/combined_curated_TE_lib_ATOSZM_selected.fasta.gz) was created by concatanate TE library from a) _Arabidopsis thaliana_, b) _Oryza sativa_, and c) _Zea maize_ before using bbmap2 to remove exact duplicates and using CDHIT to generate a non-redundat TE lib for simulation test. The A. thaliana and Z. maize TE liraries were obtained from [EDTA repo (Ou et al. 2019)](https://github.com/oushujun/EDTA/tree/master/database) and the O. sativa TE lib was from [RepeatModeller2 (Flynn et al. 2020)](https://github.com/jmf422/TE_annotation/tree/master/benchmark_libraries/RM2).
- `-m` or `--maxcp`: an interger representing the maximum copy to be simulated for a TE family
- `-n` or `--mincp`: an integer representing the minimum copy to be simulated for a TE family
- `-c` or `--chridx`: the chromosome index file if **Random Synthesized Genome Mode** (`--mode 0`) is enabled. See the previous section **Two modes**.
- `-g` or `--genome`: the genome fasta file if **Custom Genome Mode** (`--mode 1`) is enabled. See the previous section **Two modes**.
- `-s` or `--seed`: an interger to specify the seed value. Default is 1.
- `-o` or `--outdir`: the path to the output directory

### Two steps
1) **Step 1**: The first step of TEgenomeSimulator is to perform non-overlap random TE insertion into the random synthesized or custom genome sequence.
2) **Step 2**: The second step is to perform nested TE insertion.

### Output files from TEgenomeSimulator
- **TElib_sim_list.table**: A tab-delemited table storing the parameters for simulating TE mutation. See example at ./test/output/ 
- **TEgenomeSimulator_{$prefix}.yml**: A configuration file for simulation. See example at ./test/output/
- **{$prefix}_genome_sequence_out.fasta**: The simulated genome fasta file with non-overlap random TE insertion (from **Step 1**).
- **{$prefix}_"_repeat_sequence_out.fasta**: The simulated TE sequences that had been inserted into the genome (from **Step 1**).
- **{$prefix}_repeat_annotation_out.gff**: The coordinates of inserted TE after non-overlap random TE insertion (from **Step 1**).
- **{$prefix}_genome_sequence_out_final.fasta**: The final simulated genome fasta file with non-overlap and nested random TE insertion (from **Step 2**).
- **{$prefix}_repeat_sequence_out_final.fasta**: The final simulated TE sequences, including the nested TEs, that had been inserted into the genome (from **Step 2**).
- **{$prefix}_repeat_annotation_out_final.gff**: The final coordinates of all inserted TE adjusted after nested TE insertion (from **Step 2**)
- **TEgenomeSimulator.log**: The log file of the simulation.


## Intsllation steps
### 1. Clone this repo
```bash
cd $MYDIR
git clone git@github.com:PlantandFoodResearch/TEgenomeSimulator.git
```

### 2. Navigate to the cloned folder and install
```bash
cd TEgenomeSimulator
pip install .
```

### 3. Usage
After installation, you can take a look at the arguments of TEgenomeSimulator by typing `tegnomesimulator --help`. 
```text
usage: tegenomesimulator [-h] -M {0,1} -p PREFIX -r REPEAT -m MAXCP -n MINCP [-c CHRIDX] [-g GENOME] [-s SEED] -o OUTDIR

main arguments of TEgenomeSimulator to simulate TE mutation and insertion into genome.

optional arguments:
  -h, --help            show this help message and exit
  -M {0,1}, --mode {0,1}
                        Mode for genome simulation (either 0 or 1).
  -p PREFIX, --prefix PREFIX
                        Prefix for output files.
  -r REPEAT, --repeat REPEAT
                        TE family fasta file.
  -m MAXCP, --maxcp MAXCP
                        Maximum copies of TE family.
  -n MINCP, --mincp MINCP
                        Minimum copies of TE family.
  -c CHRIDX, --chridx CHRIDX
                        Chromosome index file if mode 0 is selected.
  -g GENOME, --genome GENOME
                        Genome fasta file if mode 1 is selected.
  -s SEED, --seed SEED  Random seed (default is 1).
  -o OUTDIR, --outdir OUTDIR
                        Output directory.
```

## Run TEgenomeSimulator using test data
### Random Synthesized Genome Mode 
```bash
cd TEgenomeSimulator
prefix=test_random
chridx="../test/input/random_genome_chr_index.txt" 
repeat="../test/input/combined_curated_TE_lib_ATOSZM_selected.fasta"
max=5
min=1
outdir="../test/output"
mkdir -p $outdir

python3 tegenomesimulator.py -M 0 -p $prefix -c $chridx -r $repeat -m $max -n $min -o $outdir
```

### Custom Genome Mode
```bash
cd TEgenomeSimulator
prefix=test_custom
genome="../test/input/Donghong.chromosomes.only.fa.nonTE.2chrs"
repeat="../test/input/combined_curated_TE_lib_ATOSZM_selected.fasta"
max=5
min=1
outdir="../test/output"
mkdir -p $outdir

python3 tegenomesimulator.py -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $outdir
```

## References:
- [Matias Rodriguez & Wojciech Makałowski. Software evaluation for de novo detection of transposons. 2022 Mobil DNA](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-022-00266-2)
- [denovoTE-eval: https://github.com/IOB-Muenster/denovoTE-eval](https://github.com/IOB-Muenster/denovoTE-eval)
- [Ou et al. Benchmarking transposable element annotation methods for creation of a streamlined, comprehensive pipeline. 2019 Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1905-y)
- [Flynn et al. RepeatModeler2 for automated genomic discovery of transposable element families. 2020 PNAS](https://doi.org/10.1073/pnas.1921046117)

## Funding and Acknowledgement
- This work is part of the Kiwifruit Royalty Investment Programme (KRIP)-funded _Genome Landscape_ objective, led by Susan Thomson.
- The first release of TEgenomeSimulator v0.1.0 was built with the joint efforts from the Bioinformatics Analysis Team and Bioinformatics Engineering Team of Plant & Food Research. We appreciate all the advices and comments from our team members.

