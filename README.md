# TEgenomeSimulator

TEgenomeSimulator was created based on [Matias Rodriguez & Wojciech Makałowski. Software evaluation for de novo detection of transposons. 2022 Mobil DNA](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-022-00266-2). The original scripts were from [denovoTE-eval](https://github.com/IOB-Muenster/denovoTE-eval). Some existing functions had been modified and several new functions had been created for TEgenomeSimulator.\

TEgenomeSimulator allows:
1) the simulation of a random genome with multiple chromosomes that are later inserted with TEs randomly;
2) utilise a user-provided genome containing multiple chromosomes where TE bases had been removed, providing a customised genome canvas for random TE insertion.


## Test
```bash
# Clone this repo
cd $MYDIR
git clone git@github.com:PlantandFoodResearch/TEgenomeSimulator.git
```

```bash
# Random Genome Mode
ml conda
conda activate TEgenomeSimulator

cd $MYDIR/TEgenomeSimulator/TEgenomeSimulator
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

