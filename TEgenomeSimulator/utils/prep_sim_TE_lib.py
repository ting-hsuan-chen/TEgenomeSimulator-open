import os
import sys
import argparse
import random
import re
from Bio import SeqIO
from pathlib import Path


# parse argument
parser = argparse.ArgumentParser()
parser.add_argument("-p","--prefix", type=str,
                    help="project prefix of the simulation")
parser.add_argument("-r","--repeat", type=str,
                    help="path to the repeat fasta file")
parser.add_argument("-m", "--maxcp", type=int,
                    help="maximum copy number to be simulated for a TE family")
parser.add_argument("-n", "--mincp", type=int,
                    help="minimum copy number to be simulated for a TE family")
parser.add_argument('-s', '--seed', type=int, 
                    default=1, help="Random seed (default is 1).")
parser.add_argument("-o", "--outdir", type=str,
                    help="output directory")

args = parser.parse_args()
prefix = args.prefix
te_fa = args.repeat
max_cp = args.maxcp
min_cp = args.mincp
seed = args.seed
out_dir = args.outdir

print("\n")
print("#########################################################")
print("### Prepare TE library table with simulation settings ###")
print("#########################################################")
print(f"Using repeat fasta file {args.repeat}")
print(f"Output directory set as {args.outdir}")

#Set seed
if seed:
    random.seed(seed)

# Load files
te_lib = SeqIO.to_dict(SeqIO.parse(te_fa ,"fasta"))

# TE family
te_family = []
for te_id in te_lib:
    te_family.append(te_id)

# TE superfamily
te_superfamily = []
for te_id in te_family:
    superfamily = re.sub(".*#", "", te_id)
    te_superfamily.append(superfamily)
    
# TE subclass
ltr_retro = ['LTR/Copia', 'LTR/Gypsy', 'LTR/Ty3', 'LTR/Solo', 'LTR/unknown']
line = ['LINE/unknown', 'LINE/L1']
sine = ['SINE/unknown', 'SINE/tRNA']
tir = ['DNA/hAT', 'DNAnona/hAT', 'DNAauto/hAT', 'DNA/CACTA', 'DNAnona/CACTA', 'DNAauto/CACTA', 'DNA/Harbinger', 'DNA/MuDR', 'DNAnona/MULE', 'DNAauto/MULE', 'DNA/Mariner']
helitron = ['DNA/Helitron', 'DNAnona/Helitron', 'DNAauto/Helitron']
mite = ['MITE/Stow', 'MITE/Tourist']

te_subclass = []
for sup_fam in te_superfamily:
    if sup_fam in ltr_retro:
        subclass = 'LTR_retrotransposon'
    if sup_fam in line:
        subclass = 'LINE_retrotransposon'
    if sup_fam in sine:
        subclass = 'SINE_retrotransposon'
    if sup_fam in tir:
        subclass = 'TIR_transposon'
    if sup_fam in helitron:
        subclass = 'Helitron'
    if sup_fam in mite:
        subclass = 'MITE'
    te_subclass.append(subclass)
     
# count
print("\n")
print("## Random chosing copy numbe for each TE family ##")
copy_number = []
minimum = min_cp
maximum = max_cp
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    copy_number.append(random_num)
print(f"Maximum copy number set by user: {args.maxcp}")
print(f"Minimum copy number set by user: {args.maxcp}")

# identity
print("\n")
print("## Random chosing the averaged sequence identity for each TE family ##")
identity = []
minimum = 80
maximum = 95
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    identity.append(random_num)
print(f"Maximum averaged sequence identity set by default: {maximum}")
print(f"Minimum averaged sequence identity set by default: {minimum}")

# standard deviation of identity
print("\n")
print("## Random chosing the standard deviation of averaged sequence identity for each TE family ##")
sd = []
minimum = 1
maximum = 20
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    sd.append(random_num)
print(f"Maximum standard deviation of averaged sequence identity set by default: {maximum}")
print(f"Minimum standard deviation of averaged sequence identity set by default: {minimum}")
    
# indel (as a proportion to total SNP = substitution + indel)
print("\n")
print("## Random chosing the proportion of INDEL to total SNP (dependant on sequence identity) for each TE family ##")
indel = []
minimum = 5
maximum = 20
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    indel.append(random_num)
print(f"Maximum INDEL proportion set by default: {maximum}")
print(f"Minimum INDEL proportion set by default: {minimum}")
    
# tsd (use prior knowledge)
print("\n")
print("## Setting the length of std based on prior knowledge ##")
ltr_retro = ['LTR/Copia', 'LTR/Gypsy', 'LTR/Ty3', 'LTR/Solo', 'LTR/unknown'] 
line = ['LINE/unknown', 'LINE/L1']
sine = ['SINE/unknown', 'SINE/tRNA']
DTA = ['DNA/hAT', 'DNAnona/hAT', 'DNAauto/hAT']
DTC = ['DNA/CACTA', 'DNAnona/CACTA', 'DNAauto/CACTA']
DTH = ['DNA/Harbinger']
DTM = ['DNA/MuDR', 'DNAnona/MULE', 'DNAauto/MULE']
DTT = ['DNA/Mariner']
helitron = ['DNA/Helitron', 'DNAnona/Helitron', 'DNAauto/Helitron']
mite = ['MITE/Stow', 'MITE/Tourist']

tsd = []
for sup_fam in te_superfamily:
    if sup_fam in ltr_retro:
        tsd_range = '5,5'
    if sup_fam in line:
        tsd_range = '5,20'
    if sup_fam in sine:
        tsd_range = '5,20'
    if sup_fam in DTA:
        tsd_range = '5,8'
    if sup_fam in DTC:
        tsd_range = '2,4'
    if sup_fam in DTH:
        tsd_range = '3,3'
    if sup_fam in DTM:
        tsd_range = '8,9'
    if sup_fam in DTT:
        tsd_range = '2,2'
    if sup_fam in helitron:
        tsd_range = '0,0'
    if sup_fam in mite:
        tsd_range = '2,10'
    tsd.append(str(tsd_range))

print("Length range of TSD for LTR retrotransposon set to: 5 - 5")
print("Length range of TSD for LINE set to: 5 - 20")
print("Length range of TSD for SINE set to: 5 - 20")
print("Length range of TSD for DTA set to: 5 - 8")
print("Length range of TSD for DTC set to: 2 - 4")
print("Length range of TSD for DTH set to: 3 - 3")
print("Length range of TSD for DTM set to: 8 - 9")
print("Length range of TSD for DTT set to: 2 - 2")
print("Length range of TSD for Helitron set to: 0")
print("Length range of TSD for MITE set to: 2 - 10")

# length
print("\n")
print("## Extracting the length of each TE family ##")
length = []
for te_id in te_lib:
    te_len = len(te_lib[te_id].seq)
    length.append(te_len)

# fragmented TE loci (as a proportion to total TE loci of the family)
print("\n")
print("## Setting the proportion of fragmented TE loci of each TE family ##")
fragment = []
minimum = 50
maximum = 98
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    fragment.append(random_num)
print(f"Maximum proportion of fragmented TE loci of each TE family set by default: {maximum}")
print(f"Minimum proportion of fragmented TE loci of each TE family set by default: {minimum}")

# nest TE loci (as a proportion to total TE loci of the family; only apply for Copia and Gypsy)
print("\n")
print("## Setting the proportion of nested TE insertion of each Copia or Gypsy family ##")
nested = []
minimum = 0
maximum = 30
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    nested.append(random_num)
print(f"Maximum proportion of nested TE insertion of each Copia or Gypsy family set by default: {maximum}")
print(f"Minimum proportion of nested TE insertion of each Copia or Gypsy family set by default: {minimum}")

# create the table
print("\n")
print("## Printing the TE library table ##")

final_out = str(out_dir) + "/TEgenomeSimulator_" + str(prefix) + "_result"
#Path(out_dir, "TEgenomeSimulator_" + prefix + "_result").mkdir(parents=True, exist_ok=True)
#os.chdir(Path(out_dir, "TEgenomeSimulator_" + file_prefix + "_result"))

table_out = open(Path(final_out, "TElib_sim_list.table"), "w")
table_out.write("\t".join(["#TE_family", "superfamily", "subclass", "count", "idn", "sd", "indels", "tsd", "length", "frag", "nest" + "\n"]))
n = len(te_family)
for i in range(n - 1):
    table_out.write("\t".join([str(te_family[i]), str(te_superfamily[i]), str(te_subclass[i]), str(copy_number[i]), str(identity[i]), str(sd[i]), str(indel[i]), str(tsd[i]), str(length[i]), str(fragment[i]), str(nested[i]) + "\n"]))
table_out.close()    

print(f"Generated the TE library table for simulation. File saved as {out_dir}/TElib_sim_list.table")
print("\n")