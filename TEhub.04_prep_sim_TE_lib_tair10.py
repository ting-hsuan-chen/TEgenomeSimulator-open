import os
import sys
import random
import re
import pandas as pd
from Bio import SeqIO

# Load files
#te_fa = str(sys.argv[1])
rm_outsum = "/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator/model_species_RM/tair10/AthaGenome.fa.out.summary"
te_sum = pd.read_table(rm_outsum, sep='\t')

# TE family is the column 'original_family'

# TE superfamily is the colmn 'superfamily'

# TE subclass
te_sum['subclass'] = te_sum['superfamily'].replace("/.*", "", regex=True)  
# count is the column 'copynumber'

# identity
te_sum['idn'] = round(te_sum['mean_identity'] * 100, 0)
te_sum['idn'] = te_sum['idn'].astype(int)

# standard deviation of identity
te_sum['sd'] = round(te_sum['sd_identity'] * 100, 0)
te_sum['sd'] = te_sum['sd'].astype(int)    

# indel (as a proportion to total SNP = substitution + indel)
te_sum['indels'] = round(te_sum['mean_indel'] * 100, 0)
te_sum['indels'] = te_sum['indels'].astype(int)

# tsd (use prior knowledge)
te_superfamily = list(te_sum['superfamily'])
ltr_retro = ['LTR/Copia', 'LTR/Gypsy', 'LTR/Ty3', 'LTR/Solo', 'LTR/unknown'] 
line = ['LINE/unknown', 'LINE/L1', 'LINE/LINE_Unknown']
sine = ['SINE/unknown', 'SINE/tRNA', 'SINE/SINE_Unknown']
DTA = ['DNA/hAT', 'DNAnona/hAT', 'DNAauto/hAT', 'DNA/DTA']
DTC = ['DNA/CACTA', 'DNAnona/CACTA', 'DNAauto/CACTA', 'DNA/DTC']
DTH = ['DNA/Harbinger', 'DNA/DTH']
DTM = ['DNA/MuDR', 'DNAnona/MULE', 'DNAauto/MULE', 'DNA/DTM']
DTT = ['DNA/Mariner', 'DNA/DTT']
helitron = ['DNA/Helitron', 'DNAnona/Helitron', 'DNAauto/Helitron', 'RC/Helitron']
mite = ['MITE/Stow', 'MITE/Tourist']

tsd = []
for sup_fam in te_superfamily:
    if sup_fam in ltr_retro:
        tsd_range = '5,5'
    elif sup_fam in line:
        tsd_range = '5,20'
    elif sup_fam in sine:
        tsd_range = '5,20'
    elif sup_fam in DTA:
        tsd_range = '5,8'
    elif sup_fam in DTC:
        tsd_range = '2,4'
    elif sup_fam in DTH:
        tsd_range = '3,3'
    elif sup_fam in DTM:
        tsd_range = '8,9'
    elif sup_fam in DTT:
        tsd_range = '2,2'
    elif sup_fam in helitron:
        tsd_range = '0,0'
    elif sup_fam in mite:
        tsd_range = '2,10'
    else:
        tsd_range = '0,0' 
    tsd.append(str(tsd_range))

te_sum['tsd'] = pd.DataFrame({'tsd': tsd})

# length is the column 'full_length_bp'

# fragmented TE loci (as a proportion to total TE loci of the family)
# need adjust the propostion of fragmented loci
te_sum['new_frag'] = (te_sum['fragment_ratio'] + 1) / 2
te_sum['frag'] = round(te_sum['new_frag'] * 100, 0)
te_sum['frag'] = te_sum['frag'].astype(int)

# nest TE loci (as a proportion to total TE loci of the family; only apply for Copia and Gypsy)
copia_gypsy = ['LTR/Copia', 'LTR/Gypsy', 'LTR/Ty3'] 

nested = []
minimum = 0
maximum = 30
for sup_fam in te_superfamily:
    if sup_fam in copia_gypsy:
        nest_value = random.choice(range(minimum, maximum + 1))
    else:
        nest_value = int(0)
    nested.append(nest_value)

te_sum['nest'] = pd.DataFrame({'nest': nested})
 
# create the table
df = te_sum[['original_family', 'subclass', 'superfamily', 'copynumber', 'idn', 'sd', 'indels', 'tsd', 'full_length_bp', 'frag', 'nest']]
df = df.rename(columns={'original_family': '#TE_family', 'copynumber': 'count', 'full_length_bp': 'length'})
df.to_csv("TElib_sim_list_tair10.table", sep='\t', index=False, header=True)
