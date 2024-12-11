# TEgenomeSimulator

TEgenomeSimulator was created based on [Matias Rodriguez & Wojciech Makałowski. Software evaluation for de novo detection of transposons. 2022 Mobil DNA](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-022-00266-2). The original scripts were from [denovoTE-eval](https://github.com/IOB-Muenster/denovoTE-eval). Some existing functions had been modified and several new functions had been created for TEgenomeSimulator.\

TEgenomeSimulator allows:
1) the simulation of a random genome with multiple chromosomes that are later inserted with TEs randomly;
2) utilise a user-provided genome containing multiple chromosomes where TE bases had been removed, providing a customised genome canvas for random TE insertion.


## Intsllation
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

## Functions of TEgenomeSimulator
After installation, you can take a look at the arguments of TEgenomeSimulator by typing `tegnomesimulator --help`. 
```
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

## Run TEgenomeSimulator
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

python3 TEgenomeSimulator.py -M 0 -p $prefix -c $chridx -r $repeat -m $max -n $min -o $outdir
```

```bash
# Custom Genome Mode
ml conda
conda activate TEgenomeSimulator

cd $MYDIR/TEgenomeSimulator/TEgenomeSimulator
prefix=test_custom
genome="../test/input/Donghong.chromosomes.only.fa.nonTE.2chrs"
repeat="../test/input/combined_curated_TE_lib_ATOSZM_selected.fasta"
max=5
min=1
outdir="../test/output"
mkdir -p $outdir

python3 TEgenomeSimulator.py -M 1 -p $prefix -g $genome -r $repeat -m $max -n $min -o $outdir
```

Note: the conda environment of TEgenomeSimulator was build as follow:

```bash
module load conda
conda create --name TEgenomeSimulator python=3.9
conda activate TEgenomeSimulator
pip3 install biopython pyyaml numpy
pip install pandas
```

