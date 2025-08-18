import os
import sys
import argparse
import random
import re
from Bio import SeqIO
from pathlib import Path
import numpy as np


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
parser.add_argument('--maxidn', type=int, 
                    help="The upper bound of mean sequence identity to be sampled for each TE family.")
parser.add_argument('--minidn', type=int, 
                    help="The lower bound of mean sequence identity to be sampled for each TE family.")
parser.add_argument('--maxsd', type=int, 
                    help="The upper bound of standard deviation of mean identity to be sampled for each TE family.")
parser.add_argument('--minsd', type=int, 
                    help="The lower bound of standard deviation of mean identity to be sampled for each TE family.")
parser.add_argument("-i", "--intact", type=float, help="maximum proportion of intact TEs per TE family (default is 0.001, i.e. 0.1%)")
parser.add_argument('-s', '--seed', type=int, 
                    default=1, help="Random seed (default is 1).")
parser.add_argument("-o", "--outdir", type=str,
                    help="output directory")

args = parser.parse_args()
prefix = args.prefix
te_fa = args.repeat
max_cp = args.maxcp
min_cp = args.mincp
max_idn = args.maxidn
min_idn = args.minidn
max_sd = args.maxsd
min_sd = args.minsd
intact = args.intact
seed = args.seed
#out_dir = args.outdir
final_out = args.outdir

print("\n", flush=True)
print("#########################################################", flush=True)
print("### Prepare TE library table with simulation settings ###", flush=True)
print("#########################################################", flush=True)
print(f"Using repeat fasta file {args.repeat}", flush=True)
print(f"Output directory set as {args.outdir}", flush=True)

#Set seed
if seed:
    random.seed(seed)
    np.random.seed(seed)

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
te_subclass = []
for sup_fam in te_superfamily:
    subclass = re.sub("/.*", "", sup_fam)
    #sup_fam.replace("/.*", "", regex=True)
    te_subclass.append(subclass)

# te_subclass = []
# ltr_retro = ['LTR/Copia', 'LTR/Gypsy', 'LTR/Ty3', 'LTR/Solo', 'LTR/unknown']
# line = ['LINE/unknown', 'LINE/L1']
# sine = ['SINE/unknown', 'SINE/tRNA']
# tir = ['DNA/hAT', 'DNAnona/hAT', 'DNAauto/hAT', 'DNA/CACTA', 'DNAnona/CACTA', 'DNAauto/CACTA', 'DNA/Harbinger', 'DNA/MuDR', 'DNAnona/MULE', 'DNAauto/MULE', 'DNA/Mariner']
# helitron = ['DNA/Helitron', 'DNAnona/Helitron', 'DNAauto/Helitron']
# mite = ['MITE/Stow', 'MITE/Tourist']
# 
# te_subclass = []
# for sup_fam in te_superfamily:
#     if sup_fam in ltr_retro:
#         subclass = 'LTR_retrotransposon'
#     if sup_fam in line:
#         subclass = 'LINE_retrotransposon'
#     if sup_fam in sine:
#         subclass = 'SINE_retrotransposon'
#     if sup_fam in tir:
#         subclass = 'TIR_transposon'
#     if sup_fam in helitron:
#         subclass = 'Helitron'
#     if sup_fam in mite:
#         subclass = 'MITE'
#     te_subclass.append(subclass)
     
# count
print("\n", flush=True)
print("## Random chosing copy numbe for each TE family ##", flush=True)
copy_number = []
minimum = min_cp
maximum = max_cp
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    copy_number.append(random_num)
print(f"Maximum copy number set by user: {args.maxcp}", flush=True)
print(f"Minimum copy number set by user: {args.maxcp}", flush=True)

# identity
print("\n", flush=True)
print("## Random chosing the averaged sequence identity for each TE family ##", flush=True)
identity = []
minimum = min_idn #default 80
maximum = max_idn #default 95
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    identity.append(random_num)
print(f"Maximum averaged sequence identity: {maximum}", flush=True)
print(f"Minimum averaged sequence identity: {minimum}", flush=True)

# standard deviation of identity
print("\n", flush=True)
print("## Random chosing the standard deviation of averaged sequence identity for each TE family ##", flush=True)
sd = []
minimum = min_sd #default 1
maximum = max_sd #default 20
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    sd.append(random_num)
print(f"Maximum standard deviation of averaged sequence identity: {maximum}", flush=True)
print(f"Minimum standard deviation of averaged sequence identity: {minimum}", flush=True)
    
# indel (as a proportion to total SNP = substitution + indel)
print("\n", flush=True)
print("## Random chosing the proportion of INDEL to total SNP (dependant on sequence identity) for each TE family ##", flush=True)
indel = []
minimum = 5
maximum = 20
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    indel.append(random_num)
print(f"Maximum INDEL proportion set by default: {maximum}", flush=True)
print(f"Minimum INDEL proportion set by default: {minimum}", flush=True)
    
# tsd (use prior knowledge)
print("\n", flush=True)
print("## Setting the length of std based on prior knowledge ##", flush=True)
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
    else:
        tsd_range = '0,0'
    tsd.append(str(tsd_range))

print("Length range of TSD for LTR retrotransposon set to: 5 - 5", flush=True)
print("Length range of TSD for LINE set to: 5 - 20", flush=True)
print("Length range of TSD for SINE set to: 5 - 20", flush=True)
print("Length range of TSD for DTA set to: 5 - 8", flush=True)
print("Length range of TSD for DTC set to: 2 - 4", flush=True)
print("Length range of TSD for DTH set to: 3 - 3", flush=True)
print("Length range of TSD for DTM set to: 8 - 9", flush=True)
print("Length range of TSD for DTT set to: 2 - 2", flush=True)
print("Length range of TSD for Helitron set to: 0", flush=True)
print("Length range of TSD for MITE set to: 2 - 10", flush=True)
print("Length range of TSD for else set to: 0", flush=True)

# length
print("\n", flush=True)
print("## Extracting the length of each TE family ##", flush=True)
length = []
for te_id in te_lib:
    te_len = len(te_lib[te_id].seq)
    length.append(te_len)

# fragmented TE loci (as a proportion to total TE loci of the family)
print("\n", flush=True)
print("## Setting the proportion of fragmented TE loci of each TE family ##", flush=True)
#min_intact_ratio = 0
#max_intact_ratio = float(intact)
# the value of "intact" represent the maximum chance of 
fragment = []
minimum = 100 - intact
maximum = 100
n = len(te_family)
for i in range(n):
    random_num = np.random.uniform(minimum, maximum)
    #random_num = random.choice(range(minimum, maximum + 1))
    fragment.append(random_num)
print(f"Maximum chance of keeping a TE insertion intact as 100% integrity in each TE family: {intact}", flush=True)
print(f"Minimum chance of keeping a TE insertion intact in each TE family: 0", flush=True)
print(f"Maximum proportion of fragmented TE loci of each TE family: {maximum}", flush=True)
print(f"Minimum proportion of fragmented TE loci of each TE family: {minimum}", flush=True)

# nest TE loci (as a proportion to total TE loci of the family; only apply for Copia and Gypsy)
print("\n", flush=True)
print("## Setting the proportion of nested TE insertion of each Copia or Gypsy family ##", flush=True)
nested = []
minimum = 0
maximum = 30
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    nested.append(random_num)
print(f"Maximum proportion of nested TE insertion of each Copia or Gypsy family set by default: {maximum}", flush=True)
print(f"Minimum proportion of nested TE insertion of each Copia or Gypsy family set by default: {minimum}", flush=True)

# create the table
print("\n", flush=True)
print("## Printing the TE library table ##", flush=True)

table_out = open(Path(final_out, "TElib_sim_list.table"), "w")
table_out.write("\t".join(["#TE_family", "subclass", "superfamily", "count", "idn", "sd", "indels", "tsd", "length", "frag", "nest" + "\n"]))
n = len(te_family)
for i in range(n - 1):
    table_out.write("\t".join([str(te_family[i]), str(te_subclass[i]), str(te_superfamily[i]),  str(copy_number[i]), str(identity[i]), str(sd[i]), str(indel[i]), str(tsd[i]), str(length[i]), str(fragment[i]), str(nested[i]) + "\n"]))
table_out.close()    

print(f"Generated the TE library table for simulation. File saved as {final_out}/TElib_sim_list.table", flush=True)
print("\n", flush=True)