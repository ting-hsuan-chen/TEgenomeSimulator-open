import os
import sys
import random
import yaml
import re
from Bio import SeqIO, Seq
from pathlib import Path
from custom_genome_seq_TEs import parse_random_genome_yaml, parse_custom_genome_yaml, load_repeats_chr, generate_mismatches, add_indels, add_base_changes
from custom_genome_seq_TEs import get_identity, create_TSD, fragment
#from random_sequence_TEs import load_repeats, generate_mismatches, add_indels, add_base_changes
#from random_sequence_TEs import get_identity, create_TSD, fragment

#Load params from YAML file config.yml in same directory
def parse_yaml():
    params = yaml.load(open('config.yml', 'r'), Loader=yaml.FullLoader)
    return params

#Turn gff table into a big list where each item represent a row in gff table
def load_gff(gff_file):
    gff = []
    with open(gff_file) as gff_fh:
        for i in gff_fh:
            e = i.strip().split("\t")
            gff.append(e)
    return gff

#Modify the coordinates of TE loci locating downstream of the nested insertion 
def modify_coords(offset, index, new_gff):
    new_gff_aux = []
    for i in range(index, len(new_gff)):
            start = int(new_gff[i][3]) + offset
            new_gff[i][3] = str(start)
            end = int(new_gff[i][4]) + offset
            new_gff[i][4] = str(end)
    return new_gff

#Modify the coordinates of TE loci locating downstream of the nested insertion for genome with multiple chr (function created by THC)
def modify_genome_coords(offset, index, new_gff, chr_id):
    new_gff_aux = []
    for i in range(index, len(new_gff)):
        if new_gff[i][0] == chr_id:
            start = int(new_gff[i][3]) + offset
            new_gff[i][3] = str(start)
            end = int(new_gff[i][4]) + offset
            new_gff[i][4] = str(end)
    return new_gff

def filter_nonest(gff):
    c = 0
    vec_cand = []
    for i in gff:
        if "Alu" in i[8] or "SINE" in i[8]:
            pass
        else:
            vec_cand.append(c)
        c+=1
    return vec_cand

#Turn the inserted repeat fasta file into a dictionary with te_id as key (function created by THC)
def load_isrt_te_fa(inserted_te_fasta_file):
    isrt_te_dict=SeqIO.to_dict(SeqIO.parse(inserted_te_fasta_file,"fasta"))    
    keys_copy = list(isrt_te_dict.keys())
    for key in keys_copy:
        new_key = re.sub(".*_TE", "TE", key)
        new_key = re.sub("#.*","", new_key)
        isrt_te_dict[new_key] = isrt_te_dict.pop(key)    
    return isrt_te_dict

#Generate vector of coords for base_changes and indels
def generate_nests(repeats, gff, seq, prefix):
    rep_count = []
    new_seq = ""
    #Create a list containing the TE family and the required copies for nested insertion
    for i in repeats:
        rep_count += [i]*int(repeats[i].num_rep*(repeats[i].nest/100.0))
    #Shuffle the list
    random.shuffle(rep_count)
    #Acquire the index of TEs from gff file (excluding Alu and SINE)
    vec_cand = filter_nonest(gff)
    #Randomly pick the index of TE loci to be inserted with nested TEs
    #insert_index = random.sample(range(0, len(gff)-1), len(rep_count))
    insert_index = random.sample(vec_cand, len(rep_count))
    new_gff = gff
    sorted_table = sorted(zip(insert_index, rep_count))
    counter = 0
    n = 0
    for j,k in zip(sorted_table, rep_count):
        gff_sel = new_gff[j[0] + counter]
        start = int(gff_sel[3])
        end = int(gff_sel[4])
        length = end - start +1
        
        #Decide on the position of nested insertion in association to the length of host TE locus 
        pct_pos = random.randint(40,60)
        ins_pos = int(round((pct_pos/100.0) * length))
        
        #Create SNPs, indels and TSD for the nested TE
        nest_seq = repeats[k].sequence
        nest_identity = get_identity(repeats[k].identity, repeats[k].sd)
        nest_identity_fix = nest_identity + (100 - nest_identity) * 0.5
        nest_indels = repeats[k].indels
        base_changes_vec, indels_changes_vec = generate_mismatches(nest_seq, nest_identity_fix, nest_indels)
        nest_seq_mismatches = add_base_changes(nest_seq, base_changes_vec)
        new_nest_seq = add_indels(nest_seq_mismatches, indels_changes_vec)

        new_nest_seq_tsd = new_nest_seq
        tsd_5_len = tsd_3_len = 0
        if repeats[k].tsd:
            tsd_seq_5, tsd_seq_3 = create_TSD(nest_identity_fix, nest_indels)
            new_nest_seq_tsd = tsd_seq_5 + new_nest_seq + tsd_seq_3
            tsd_5_len = len(tsd_seq_5)
            tsd_3_len = len(tsd_seq_3)

        #Fragment weighted (applying two-third chance of fragmentation)
        isFrag = random.choice([1,1,0])
        new_nest_seq_tsd_frag = new_nest_seq_tsd
        if isFrag:
            new_nest_seq_tsd_frag, frag, cut = fragment(new_nest_seq_tsd)

        nest_len = len(new_nest_seq_tsd_frag)
        nest_name = repeats[k].name

        #Calculate the coordinates of the host TEs and nested TEs after insertion
        new_end_1 = start + ins_pos
        new_start_2 = new_end_1 + nest_len 
        new_end_2 = new_start_2 + (length - ins_pos)

        #Apply strand sense
        strands = ["+", "-"]
        strand = random.choice(strands)
        new_nest_seq_str = new_nest_seq_tsd_frag
        
        if strand == "-":
            new_nest_seq_str = str(Seq.Seq(new_nest_seq_tsd_frag).reverse_complement())

        #Prepare updated content to be put into gff list    
        frag_note = ""
        if isFrag:
            frag_note = ";fragment=" + str(frag)
        ori_name_1 = [gff_sel[8].replace(";ide","_1;ide") + ";note=cut_" + str(pct_pos)]
        ori_line_1 = gff_sel[:3] + [start] + [new_end_1] + gff_sel[5:8] + ori_name_1

        nest_name_in = ["ID=" + nest_name + "_n" + str(n) + ";identity=" + str(nest_identity) + ";note=nested" + frag_note ]
        nested_line = gff_sel[:3] + [new_end_1+1+tsd_5_len] + [new_start_2-tsd_3_len] + [".\t" + strand + "\t."] + nest_name_in
        
        ori_name_2 = [gff_sel[8].replace(";ide","_2;ide") + ";note=cut_" + str(100-pct_pos)]
        ori_line_2 = gff_sel[:3] + [new_start_2] + [new_end_2 -1] + gff_sel[5:8] + ori_name_2
        
        n += 1
        
        #Update the entire gff list
        index = j[0] + counter
        new_gff.pop(index)
        new_gff_aux = new_gff[:index]
        new_gff_aux.append(ori_line_1)
        new_gff_aux.append(nested_line)
        new_gff_aux.append(ori_line_2)
        new_gff_aux += new_gff[index:]
        new_gff = new_gff_aux
        new_gff = modify_coords(nest_len, index+3, new_gff)
        counter += 2
        
        #Update the genome sequence
        new_seq = seq[:new_end_1] + new_nest_seq_str + seq[new_end_1:] 
        seq = new_seq
    #Output fasta files of the genome, repeats and the gff file
    print_data(prefix, seq, new_gff)

#Generate genome with nested insertions (function created by THC)
def generate_genome_nests(repeats, isrt_te_dict, gff, genome):
    rep_count = []
    new_seq = ""
    nest_te_dict = {}
    
    #Create a list containing the TE family and the required copies for nested insertion
    for i in repeats:
        rep_count += [i]*int(repeats[i].num_rep*(repeats[i].nest/100.0))
    
    #Acquire the index of TEs from gff file (excluding Alu and SINE)
    vec_cand = filter_nonest(gff)
    
    #Randomly pick the index of TE loci to be inserted with nested TEs
    insert_index = random.sample(vec_cand, len(rep_count))
    
    new_gff = gff
    sorted_table = sorted(zip(insert_index, rep_count))
    counter = 0
    n = 1
    
    for j,k in zip(sorted_table, rep_count):
        gff_sel = new_gff[j[0] + counter]
        chr_id = gff_sel[0]
        start = int(gff_sel[3])
        end = int(gff_sel[4])
        length = end - start +1
        
        #Decide on the position of nested insertion in association to the length of host TE locus 
        pct_pos = random.randint(40,60)
        ins_pos = int(round((pct_pos/100.0) * length))
                        
        #Get family name, subclass and superfamily of nested TE
        nest_name = repeats[k].name
        nest_subclas = repeats[k].subclass
        nest_superfam = repeats[k].superfamily 
        
        #Create SNPs, indels and TSD for the nested TE
        nest_seq = repeats[k].sequence
        nest_identity = get_identity(repeats[k].identity, repeats[k].sd)
        nest_identity_fix = nest_identity + (100 - nest_identity) * 0.5
        nest_indels = repeats[k].indels
        base_changes_vec, indels_changes_vec = generate_mismatches(nest_seq, nest_identity_fix, nest_indels)
        nest_seq_mismatches = add_base_changes(nest_seq, base_changes_vec)
        new_nest_seq = add_indels(nest_seq_mismatches, indels_changes_vec)
        
        new_nest_seq_tsd = new_nest_seq
        tsd_5_len = tsd_3_len = 0
        if repeats[k].tsd != [0, 0]:
            TSD_min = repeats[k].tsd[0]
            TSD_max = repeats[k].tsd[1]
            tsd_seq_5, tsd_seq_3 = create_TSD(TSD_min, TSD_max, nest_identity_fix, nest_indels)
            new_nest_seq_tsd = tsd_seq_5 + new_nest_seq + tsd_seq_3
            tsd_5_len = len(tsd_seq_5)
            tsd_3_len = len(tsd_seq_3)
        
        #Fragment weighted (applying two-third chance of fragmentation)
        isFrag = random.choice([1,1,0])
        new_nest_seq_tsd_frag = new_nest_seq_tsd
        
        if isFrag:
            new_nest_seq_tsd_frag, frag, cut = fragment(new_nest_seq_tsd)
        else:
            frag = 100
            
        nest_len = len(new_nest_seq_tsd_frag)
        #nest_name = repeats[k].name
    
        #Calculate the coordinates of the host TEs and nested TEs after insertion
        new_end_1 = start + ins_pos
        new_start_2 = new_end_1 + nest_len 
        new_end_2 = new_start_2 + (length - ins_pos)
    
        #Apply strand sense
        strands = ["+", "-"]
        strand = random.choice(strands)
        new_nest_seq_str = new_nest_seq_tsd_frag
        
        if strand == "-":
            new_nest_seq_str = str(Seq.Seq(new_nest_seq_tsd_frag).reverse_complement())
    
        #Prepare updated content to be put into gff list    
        nested_te_id = str(n).zfill(6) #prints at least 6 characters wide; i.e. at most 999,999 nested TE insertions
        nested_te_id = "TEn" + nested_te_id
        frag_note = ""
        frag_note = ";Integrity=" + str((frag/100))
        
        ori_seq_te_id = re.sub(".*_TE", "TE", gff_sel[8])
        ori_seq_te_id = re.sub(";Name.*", "", ori_seq_te_id)
        
        ori_name_1 = [gff_sel[8].replace(";Name","_1;Name").replace(";Clas","_1;Clas") + ";Cut_at=" + str((pct_pos/100)) + ";Cut_by=" + nest_name + "_" + nested_te_id]
        ori_line_1 = gff_sel[:3] + [start] + [new_end_1] + gff_sel[5:8] + ori_name_1
        
        nest_name_in = ["ID=" + nest_name + "_" + nested_te_id + ";Name=" + nested_te_id + ";Classification=" + nest_superfam + ";Identity=" + str((nest_identity/100)) + frag_note + ";Nest_in=" + ori_seq_te_id]
        
        # Nested insertion that has undergone fragmentation does not have 5' tsd
        if frag == 100:
            nested_line = gff_sel[:2] + [nest_subclas] + [new_end_1+1+tsd_5_len] + [new_start_2-tsd_3_len] + [".\t" + strand + "\t."] + nest_name_in
        else:
            nested_line = gff_sel[:2] + [nest_subclas] + [new_end_1+1] + [new_start_2-tsd_3_len] + [".\t" + strand + "\t."] + nest_name_in
     
        ori_name_2 = [gff_sel[8].replace(";Name","_2;Name").replace(";Clas","_2;Clas") + ";Cut_at=" + str((pct_pos/100)) + ";Cut_by=" + nest_name + "_" + nested_te_id]
        ori_line_2 = gff_sel[:3] + [new_start_2] + [new_end_2 -1] + gff_sel[5:8] + ori_name_2
    
        n += 1
    
        #Update the entire gff list
        index = j[0] + counter
        new_gff.pop(index)
        new_gff_aux = new_gff[:index]
        new_gff_aux.append(ori_line_1)
        new_gff_aux.append(nested_line)
        new_gff_aux.append(ori_line_2)
        new_gff_aux += new_gff[index:]
        new_gff = new_gff_aux
        new_gff = modify_genome_coords(nest_len, index+3, new_gff, chr_id)
        counter += 2
        
        #Update the all inserted repeat seq dictionary      
        ori_seq_name_old = isrt_te_dict[ori_seq_te_id].description
        ori_seq_name_old_head = re.sub("-.*", "-", ori_seq_name_old)
        ori_seq_name_old_tail = re.sub(".*;I", ";I", ori_seq_name_old)
        ori_seq_name_old_tail = re.sub("]", "", ori_seq_name_old_tail)
        nest_note = ";Cut_at=" + str((pct_pos/100)) + ";Cut_by=" + nest_name + "_" + nested_te_id + "]"
        ori_seq_name_new = ori_seq_name_old_head + str(new_end_2 -1) + ori_seq_name_old_tail + nest_note

        nested_seq_name = nest_name + "_" + nested_te_id + "#" + nest_superfam + " [Location=" + chr_id + ":" + str(new_end_1+1+tsd_5_len) + "-" + str(new_start_2-tsd_3_len) + ";Identity=" + str(nest_identity/100) + frag_note + ";Nest_in=" + ori_seq_te_id + "]"
        
        isrt_te_dict[ori_seq_te_id].description = ori_seq_name_new
        
        #Create and update dictionary for nested TE seq
        nest_te_dict[nested_te_id] = {'id': nested_seq_name, 'seq': new_nest_seq_str[tsd_5_len:(nest_len-tsd_3_len)]}
                
        #Update the genome sequence                
        seq = str(genome[chr_id].seq)
        new_seq = seq[:new_end_1] + new_nest_seq_str + seq[new_end_1:] 
        genome[chr_id].seq = new_seq
    
    return genome, isrt_te_dict, nest_te_dict, new_gff
    
    
#Print final sequence to stdout
def print_data(prefix, seq, new_gff):
    fasta_out = open(prefix + "_out_sequence_nest.fasta", "w")
    fasta_out.write( ">sequence_nest\n" )
    for n in range(0,len(seq),100):
        fasta_out.write(str(seq[n:n+100]) + "\n")
    fasta_out.close()
    gff_out = open(prefix + "_out_repeats_nest.gff", "w")
    for i in new_gff:
        gff_out.write("\t".join(map(str,i)) + "\n")
    gff_out.close()

#Print final sequence to files (function created by THC)
def print_genome_nest_data(genome, isrt_te_dict, nest_te_dict, new_gff, params):
    #Setup output directory   
    file_prefix = str(params['prefix'])
    Path("./result/sim_" + file_prefix +"_genome").mkdir(parents=True, exist_ok=True)
    os.chdir("./result/sim_" + file_prefix +"_genome")
    
    #Output file prefix # temp modification with "v2"
    file_prefix = str(params['prefix'])
    genome_fa = file_prefix + "_genome_sequence_out_nest.fasta"
    te_fa = file_prefix + "_repeat_sequence_out_nest.fasta"
    te_gff = file_prefix + "_repeat_annotation_out_nest.gff"
    
    #For genome fasta file
    fasta_out = open(genome_fa, "w")
    for chromosome in genome:
        seq = str(genome[chromosome].seq)
        fasta_out.write(">" + chromosome + "\n" + seq + "\n")
    fasta_out.close()
    
    #For all inserted TE sequences (including nested TEs)
    te_fa_out = open(te_fa, "w")
    for te in isrt_te_dict:
        header = str(isrt_te_dict[te].description)
        seq = str(isrt_te_dict[te].seq)
        te_fa_out.write(">" + header + "\n" + seq + "\n")
    for nested_te in nest_te_dict:
        header = nest_te_dict[nested_te]['id']
        seq = nest_te_dict[nested_te]['seq']
        te_fa_out.write(">" + str(header) + "\n" + str(seq) + "\n")
    te_fa_out.close()
    
    #For the new gff file
    gff_out = open(te_gff, "w")
    for i in new_gff:
        gff_out.write("\t".join(map(str,i)) + "\n")
    gff_out.close()
    
def main():
    ##For custom genome with multiple chromosoms
    ##Load parameters
    params_chr = parse_custom_genome_yaml('config_custom_genome_DH_5_1000.yml')
    
    seed = params_chr['seed']
    if seed:
        random.seed(seed)
        
    #Specify files created from previous step that generats non-overlapping insertions
    file_prefix = str(params_chr['prefix'])
    gff_file = "./result/sim_" + file_prefix +"_genome/" + file_prefix + "_repeat_annotation_out.gff"
    fasta_file = "./result/sim_" + file_prefix +"_genome/" + file_prefix + "_genome_sequence_out.fasta"
    isrt_te_fasta = "./result/sim_" + file_prefix +"_genome/" + file_prefix + "_repeat_sequence_out.fasta"
    
    #Load fasta files into dictionaries
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file,"fasta"))
    isrt_te = load_isrt_te_fa(isrt_te_fasta)
    
    #Load non-redundant TE library and the gff file containing all TE insertions
    repeats_dict = load_repeats_chr(params_chr)
    gff = load_gff(gff_file)
    
    #Generate nested insertion
    genome, isrt_te_dict, nest_te_dict, new_gff = generate_genome_nests(repeats_dict, isrt_te, gff, genome)
    
    #Output new genome fasta, all inserted TE fasta, and GFF after nested insertion.
    print_genome_nest_data(genome, isrt_te_dict, nest_te_dict, new_gff, params_chr)

if __name__ == "__main__":
    main()
