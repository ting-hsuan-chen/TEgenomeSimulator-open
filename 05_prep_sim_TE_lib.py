import os
import sys
import random
import re
from Bio import SeqIO

# Load files
#te_fa = str(sys.argv[1])
te_fa = "/workspace/cflthc/scratch/2022_Actinidia_TE/09.10_BM_modules/curated_TElib/combined_curated_TE_lib_ATOSZM_selected.fasta"
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
copy_number = []
minimum = 10
maximum = 1000
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    copy_number.append(random_num)

# identity
identity = []
minimum = 80
maximum = 95
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    identity.append(random_num)

# standard deviation of identity
sd = []
minimum = 1
maximum = 20
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    sd.append(random_num)
    
# indel (as a proportion to total SNP = substitution + indel)
indel = []
minimum = 5
maximum = 20
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    indel.append(random_num)
    
# tsd (use prior knowledge)
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

# length
length = []
for te_id in te_lib:
    te_len = len(te_lib[te_id].seq)
    length.append(te_len)

# fragmented TE loci (as a proportion to total TE loci of the family)
fragment = []
minimum = 50
maximum = 90
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    fragment.append(random_num)

# nest TE loci (as a proportion to total TE loci of the family; only apply for Copia and Gypsy)
nested = []
minimum = 0
maximum = 30
n = len(te_family)
for i in range(n):
    random_num = random.choice(range(minimum, maximum + 1))
    nested.append(random_num)

# create the table
table_out = open("TElib_sim_list.table", "w")
table_out.write("\t".join(["#TE_family", "subclass", "superfamily", "count", "idn", "sd", "indels", "tsd", "length", "frag", "nest" + "\n"]))
n = len(te_family)
for i in range(n - 1):
    table_out.write("\t".join([str(te_family[i]), str(te_superfamily[i]), str(te_subclass[i]), str(copy_number[i]), str(identity[i]), str(sd[i]), str(indel[i]), str(tsd[i]), str(length[i]), str(fragment[i]), str(nested[i]) + "\n"]))
table_out.close()    

