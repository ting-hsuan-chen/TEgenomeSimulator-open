import os
import sys
import random
import re
import pandas as pd
from Bio import SeqIO

# Read args
prefix = sys.argv[1]
rm_outsum = sys.argv[2]
#intact = sys.argv[3]
seed = sys.argv[3]
out_dir = sys.argv[4]

print("\n")
print("#########################################################")
print("### Prepare TE library table with simulation settings ###")
print("#########################################################")
print(f"Using repeatmasker output summary file {rm_outsum}")
print(f"Output directory set as {out_dir}")

#Set seed
if seed:
    random.seed(seed)

# Load file
te_sum = pd.read_table(rm_outsum, sep='\t')

# TE family is the column 'original_family'

# TE superfamily is the colmn 'superfamily'

# TE subclass
te_sum['subclass'] = te_sum['superfamily'].replace("/.*", "", regex=True)  

# allow a factor to adjust copynumber 
# count is the column 'copynumber'
#factor = 0.7
#te_sum['copynumber'] = (te_sum['copynumber']) * factor
#te_sum['copynumber'] = te_sum['copynumber'].astype(int)   

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
#te_sum['new_frag'] = (te_sum['fragment_ratio'] + 1) / 2
#te_sum['frag'] = round(te_sum['new_frag'] * 100, 0)
#te_sum['frag'] = te_sum['frag'].astype(int)
te_sum['frag'] = te_sum['fragment_ratio'].astype(float)
te_sum['frag'] = te_sum['frag'] * 100

# integrity list
# integrities = te_sum['integrity_lst']

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
print("\n")
print("## Printing the TE library table for mode 2 ##") 

df = te_sum[['original_family', 'subclass', 'superfamily', 'copynumber', 'idn', 'sd', 'indels', 'tsd', 'full_length_bp', 'frag', 'nest', 'integrity_lst']]
df = df.rename(columns={'original_family': '#TE_family', 'copynumber': 'count', 'full_length_bp': 'length'})

# write to output
final_out = str(out_dir) + "/TEgenomeSimulator_" + str(prefix) + "_result"
output_file = os.path.join(final_out, "TElib_sim_list_mode2.table")
df.to_csv(output_file, sep='\t', index=False, header=True)

print(f"Generated the TE library table for simulation. File saved as {out_dir}/TElib_sim_list_mode2.table")
print("\n")