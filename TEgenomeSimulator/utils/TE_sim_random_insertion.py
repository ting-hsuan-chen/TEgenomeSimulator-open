import os
import sys
import argparse
import random
import yaml
import numpy
from Bio import SeqIO
from Bio import Seq
from pathlib import Path
import matplotlib.pyplot as plt
import ast

class Repeat:
    """Store the information of the TE families.
    This class has been modified by THC"""
    def __init__(self, name, subclass, superfamily, sequence, num_rep, identity, sd, indels, tsd, frag, nest, integrity_list):
        self.name = name
        self.subclass = subclass
        self.superfamily = superfamily
        self.sequence = sequence
        self.num_rep = num_rep
        self.identity = identity
        self.sd = sd
        self.indels = indels
        self.tsd = tsd
        self.frag = frag
        self.nest = nest
        self.integrities = integrity_list

##Load params_chr from YAML file config_genome.yml in same directory (modified by THC)
def parse_random_genome_yaml(config_file):
    params_chr = yaml.load(open(config_file, 'r'), Loader=yaml.FullLoader)
    return params_chr

##Load params_chr from YAML file config_custum_genome.yml in same directory (function created by THC)
def parse_custom_genome_yaml(config_file):
    params_chr = yaml.load(open(config_file, 'r'), Loader=yaml.FullLoader)
    return params_chr

##Load existing non-repeat genome sequences (function created by THC)
def load_custom_genome(params_chr):
    chr_id = list(params_chr['chrs'].keys())
    chrs_dict = {}
    fasta = SeqIO.index(params_chr['genome_fasta'], "fasta")
    for chrid in chr_id:
        chrs_dict[chrid] = fasta[chrid].seq
    return chrs_dict

##Load collection of repeats and params
def load_repeats(params):
    repeats_dict = {}
    fasta = SeqIO.index(params['rep_fasta'], "fasta")
    with open(params['rep_list'], 'r') as repeats_file:
        next(repeats_file)
        for line in repeats_file:
            elem = line.rstrip().split()
            name = elem[0]
            sequence = str(fasta[name].seq).upper()
            num_rep = int(elem[1])
            identity = int(elem[2])
            sd = int(elem[3])
            indels = int(elem[4])
            tsd = True if elem[5] == "y" else False
            frag = float(elem[7])
            nest = int(elem[8])
            #integrity_list = elem[9]
            repeat = Repeat(name, sequence, num_rep, identity, sd, indels, tsd, frag, nest)
            repeats_dict[name] = repeat
    return repeats_dict

##Load collection of repeats and params for chrs simulation for mode 0 and 1 (function created by THC)
def load_repeats_chr(params_chr):
    chr_id = list(params_chr['chrs'].keys())
    repeats_dict = {}
    fasta = SeqIO.index(params_chr['rep_fasta'], "fasta")
    with open(params_chr['rep_list'], 'r') as repeats_file:
        next(repeats_file)
        for line in repeats_file:
            elem = line.rstrip().split("\t")
            name = elem[0]
            subclass = elem[1]
            superfamily = elem[2]
            sequence = str(fasta[name].seq).upper()
            num_rep = int(elem[3])
            identity = int(elem[4])
            sd = int(elem[5])
            indels = int(elem[6])
            tsd = [int(elem[7].rstrip().split(",")[0]), int(elem[7].rstrip().split(",")[1])]
            frag = float(elem[9])
            nest = int(elem[10])
            integrity_list = []
            repeat = Repeat(name, subclass, superfamily, sequence, num_rep, identity, sd, indels, tsd, frag, nest, integrity_list)
            repeats_dict[name] = repeat
    return repeats_dict

##Load collection of repeats and params for chrs simulation for mode 2 (function created by THC)
def load_repeats_chr_m2(params_chr):
    chr_id = list(params_chr['chrs'].keys())
    repeats_dict = {}
    fasta = SeqIO.index(params_chr['rep_fasta'], "fasta")
    with open(params_chr['rep_list'], 'r') as repeats_file:
        next(repeats_file)
        for line in repeats_file:
            elem = line.rstrip().split("\t")
            name = elem[0]
            subclass = elem[1]
            superfamily = elem[2]
            sequence = str(fasta[name].seq).upper()
            num_rep = int(elem[3])
            identity = int(elem[4])
            sd = int(elem[5])
            indels = int(elem[6])
            tsd = [int(elem[7].rstrip().split(",")[0]), int(elem[7].rstrip().split(",")[1])]
            frag = float(elem[9])
            nest = int(elem[10])
            integrity_list = ast.literal_eval(elem[11])
            repeat = Repeat(name, subclass, superfamily, sequence, num_rep, identity, sd, indels, tsd, frag, nest, integrity_list)
            repeats_dict[name] = repeat
    return repeats_dict

##Generate single synthetic sequence
def generate_random_sequence(params):
    #Create DNA alphabet for random sequence
    alphabet = ["A", "T", "G", "C"]
    gc = float(params['gc_content'])/100
    weights = [(1-gc)/2, (1-gc)/2, gc/2, gc/2]
    #Generate random sequence that is going to separate the repeats
    base_sequence = "".join(numpy.random.choice(alphabet, params['seq_length'], p=weights, replace=True))
    return base_sequence

##Generate random sequence for each chromosome (function created by THC)
def generate_random_chr_sequence(params_chr):
    #Create DNA alphabet for random chr sequence
    chrs_dict = {} 
    for chr_id in list(params_chr['chrs'].keys()):
        params = params_chr['chrs'][chr_id]
        chr_seq = generate_random_sequence(params)
        chrs_dict[chr_id] = chr_seq
    return chrs_dict

##Randomly select chrs and positions to insert all the repeats (function created by THC)
def assign_chr_coord_repeats(params_chr, repeats_dict):
    #Total number of non-nested repeats 
    total_num_rep = 0
    for rep in repeats_dict:
        num_rep = int(repeats_dict[rep].num_rep - repeats_dict[rep].num_rep*(repeats_dict[rep].nest / 100.0))
        total_num_rep += num_rep
    #Total genome bp
    total_genome_bp = sum([params_chr['chrs'][chrid]['seq_length'] for chrid in params_chr['chrs']])
    #Create a dictonary containing chr as key and the cummulated chr length as value
    chr_len = list(params_chr['chrs'][chrid]['seq_length'] for chrid in params_chr['chrs'])
    chr_len_cumm = list(sum(chr_len[0:x:1]) for x in range(0, len(chr_len)+1))[1:]
    chr_id = list(params_chr['chrs'][chrid]['prefix'] for chrid in params_chr['chrs'])
    chr_cumm_dict = {}
    for i in range(0, len(chr_id)):
        ch = chr_id[i]
        end_cumm = chr_len_cumm[i]
        if i == 0:
            start_cumm = 1
        else:
            start_cumm = chr_len_cumm[i - 1] + 1
        chr_cumm_dict[ch] = (start_cumm, end_cumm)
    #Create random positions on the stitched genome and sort
    random_coords = random.sample(range(total_genome_bp-1), total_num_rep)
    random_coords.sort()
    #Bin the insertion posisions to corresponding chr
    random_repeats_coords = []
    for num in random_coords:
        for key, values in chr_cumm_dict.items():
            start_cumm = values[0]
            end_cumm = values[1]
            if num <= end_cumm:
                random_repeats_coords.append((key, num - start_cumm + 1))
                break
    return random_repeats_coords

#Shuffle the order of TE families to be inserted
def shuffle_repeats(repeats_dict):
    allnames = []
    allpositions = []
    for rep in repeats_dict:
        num_rep = int(repeats_dict[rep].num_rep - repeats_dict[rep].num_rep*(repeats_dict[rep].nest / 100.0))
        names = num_rep * [repeats_dict[rep].name]
        n_frags = int(((num_rep * repeats_dict[rep].frag) / 100 ))
        positions = [0] * num_rep
        sample_changes = random.sample(range(len(positions)), n_frags)
        for f in sample_changes:
            positions[f] += 1
        allnames += names
        allpositions += positions
    name_pos = list(zip(allnames, allpositions))
    random.shuffle(name_pos)
    return name_pos

##Get identity using a normal distribution (modified by THC)
def get_identity(mean, sd):
    #identity = int(numpy.random.normal(mean, sd, 1)) # int an array is deprecated in NumPy 1.25
    #to prevent the warning message, use the following instead
    identity = numpy.random.normal(mean, sd, 1)
    identity = int(identity[0])
    while identity > 100:
        #identity = int(numpy.random.normal(mean, sd, 1))
        identity = numpy.random.normal(mean, sd, 1)
        identity = int(identity[0])
    return identity

##Generate vector of coords for base_changes and indels
def generate_mismatches(sequence, identity, indels):
    alphabet = ["T", "G", "C", "A"]
    seq_len = len(sequence)
    seq = sequence
    #Calculate number of nucleotides that need to be changed (i.e. SNPs)
    num_changes = seq_len - int(round((seq_len * identity/100.0)))
    #Generate a vector containing locations for SNPs
    pos_changes_vec = random.sample(range(seq_len), num_changes)
    #Calculate number of SNPs that need to be changed as indels
    num_indels = int(round(num_changes * (indels/100.0)))
    #Generate a vector containing SNP locations to be changed to indels
    indel_changes_vec = random.sample(pos_changes_vec, num_indels)
    #Generate a vector containing SNP locations (removal of indel locations)
    base_changes_vec = list(set(pos_changes_vec) - set(indel_changes_vec))
    base_changes_vec.sort()
    indel_changes_vec.sort()
    return base_changes_vec, indel_changes_vec

##reate SNPs
def add_base_changes(repeat_seq, base_changes_vec):
    alphabet = ["T", "G", "C", "A"]
    repeat_seq_list = list(repeat_seq)
    for pos in base_changes_vec:
        new_base = random.choice(list(set(alphabet) - set(repeat_seq_list[pos])))
        repeat_seq_list[pos] = new_base
    new_repeat_seq =  "".join(repeat_seq_list)
    return new_repeat_seq

##Add indels
def add_indels(repeat_seq, indels_changes_vec):
    alphabet = ["T", "G", "C", "A"]
    repeat_seq_list = list(repeat_seq)
    for i in range(len(indels_changes_vec)):
        #Randomly choose 1 for insertion and 0 for deletion
        if random.choice([0,1]):
            new_base = random.choice(alphabet)
            pos = indels_changes_vec[i]
            repeat_seq_list.insert(pos, new_base)
            for j in range(len(indels_changes_vec)):
                #Increment the indel position by 1 after adding one base
                indels_changes_vec[j] +=1
        else:
            repeat_seq_list.pop(i)
            for j in range(len(indels_changes_vec)):
                #Substract the indel position by 1 after removing one base 
                indels_changes_vec[j] -=1
    new_repeat_seq =  "".join(repeat_seq_list)
    return new_repeat_seq

##Add TSD if required (modified by THC to allow customised TSD length range)
def create_TSD(tsd_min, tsd_max, identity, indels):
    alphabet = ["T", "G", "C", "A"]
    tsd_seq_5 = "".join([random.choice(alphabet) for i in range(random.choice(range(tsd_min, tsd_max+1)))])
    tsd_len = len(tsd_seq_5)
    tsd_base_changes_vec, tsd_indels_changes_vec  = generate_mismatches(tsd_seq_5, identity, indels)
    tsd_seq_mismatches = add_base_changes(tsd_seq_5, tsd_base_changes_vec)
    tsd_seq_3 = add_indels(tsd_seq_mismatches, tsd_indels_changes_vec)
    return tsd_seq_5, tsd_seq_3

##Fragmentation: beta distribution
def fragment(seq, alpha=0.5, beta=0.7):
    """Fragment TE sequence using a beta distribution."""
    len_seq = len(seq)
    # Sample from beta 
    raw_sample = numpy.random.beta(alpha, beta)
    # Convert to element size
    element_size = int(numpy.clip(raw_sample * 99 + 1, 1, 100))
    # Calculate cut length
    cut_length = int(len_seq * ((100 - element_size) / 100.0))
    return seq[cut_length:], element_size, cut_length

##Fragmentation for mode 2
def fragment_m2(seq, integrity_lst):
    """Fragment TE sequence using prior info from repeatmasker output."""
    len_seq = len(seq)
    length_ratio = random.choice(integrity_lst)
    # Make sure length_ratio <= 1
    if length_ratio > 1:
        length_ratio = 1    
    # Calculate the length of the TE sequence in bp to be removed
    cut_length = int(len_seq * (1 - length_ratio))
    return seq[cut_length:], length_ratio*100 , cut_length

##Generate new sequence including the repeats in the random one (modified by THC)
def generate_sequence(repeats_dict, rand_rep_pos, rand_seq, total_names_rep, alpha, beta, mode):
    seq = ""
    tsd_seq_5= ""
    tsd_seq_3= ""
    pre_n = 0
    n=0
    new_repeats_coord = []
    for n,m in zip(rand_rep_pos, total_names_rep):
        #Get sequence of repeat
        repeat_seq = repeats_dict[m[0]].sequence
        
        #Get family name, subclass and superfamily
        family = repeats_dict[m[0]].name
        subclas = repeats_dict[m[0]].subclass
        superfam = repeats_dict[m[0]].superfamily        
        
        #Get identity from a normal distribution
        identity = get_identity(repeats_dict[m[0]].identity, repeats_dict[m[0]].sd)
        
        #Get base_changes and indels vectors and identity
        identity_fix = identity + (100 - identity) * 0.5
        base_changes_vec, indels_changes_vec = generate_mismatches(repeats_dict[m[0]].sequence, identity_fix, repeats_dict[m[0]].indels)
        
        #Add mismatches to original repeat sequence
        repeat_seq_mismatches = add_base_changes(repeat_seq, base_changes_vec)
        
        #Add indels to original repeat sequence
        new_repeat_seq = add_indels(repeat_seq_mismatches, indels_changes_vec)

        #Check if TE creates TSDs
        if repeats_dict[m[0]].tsd != [0, 0]:
            #Generate TSD
            TSD_min = repeats_dict[m[0]].tsd[0]
            TSD_max = repeats_dict[m[0]].tsd[1]
            tsd_seq_5, tsd_seq_3 = create_TSD(TSD_min, TSD_max, identity_fix, repeats_dict[m[0]].indels)
       
        #Assign sequence to a random strand
        frag = 100 #Initiate frag size at 100% of the TE seq that has underwent identity/indel/TSD check
        #cut_length = 0 #Initiate cut size at 0%
        strands = ["+", "-"]
        strand = random.choice(strands)
       # new_repeat_seq_tsd_frag = new_repeat_seq_tsd
        new_repeat_seq_str=""
        new_repeat_seq_frag=""
        if mode == 0 or mode == 1:
            if m[1] == 1:
                new_repeat_seq_frag, frag, cut_length = fragment(new_repeat_seq, alpha, beta)
            else:
                new_repeat_seq_frag = new_repeat_seq
        if mode == 2:
            integrity_lst = repeats_dict[m[0]].integrities
            new_repeat_seq_frag, frag, cut_length = fragment_m2(new_repeat_seq, integrity_lst)
            #if m[1] == 1: # remove this line 20250415
            #    new_repeat_seq_frag, frag, cut_length = fragment(new_repeat_seq, alpha, beta) # remove this line 20250415
            #else: # remove this line 20250415
            #    new_repeat_seq_frag = new_repeat_seq # remove this line 20250415

        #Apply strand sense
        if strand == "-":
            new_repeat_seq_str = str(Seq.Seq(new_repeat_seq_frag).reverse_complement())
        else:
            new_repeat_seq_str = new_repeat_seq_frag
        #Append new repeat sequence to base sequence
        new_repeat_seq_tsd = tsd_seq_5 + new_repeat_seq_str + tsd_seq_3
        seq += rand_seq[pre_n:n] + new_repeat_seq_tsd

        #Get new repeat sequence end coordinate
        repeat_end = len(seq) - len(tsd_seq_3)
        repeat_start = repeat_end - len(new_repeat_seq_str) + 1

        #Append to vector new data about new repeat useful for a GFF
        new_repeats_coord.append([str(repeat_start), str(repeat_end), new_repeat_seq_str, identity, frag, strand, family, subclas, superfam])
        #Sets new end coordinate as start for next roung
        pre_n = n
    #At the last step add the remaining base sequence
    seq += rand_seq[n:]
    #Return sequences and repeat data
    return seq, new_repeats_coord

##Generate new sequence including the repeats in the random chr sequences (function created by THC)
def generate_genome_sequence(repeats_dict, rand_rep_pos, rand_chr_dict, shuffled_repeats, alpha, beta, mode):
    #Create placeholders for the genome seq and TE coordinates
    genome_dict = {}
    new_repeats_coord_dict = {}
    #Iterate through chromosomes
    shuffled_start_index = 0
    for chromosome in rand_chr_dict.keys():
        #Capture chromosome sequence
        rand_seq = rand_chr_dict[chromosome]
        #Capture the random position corresponding to the chromosome
        rand_rep_pos_filter = [coord for chr_id, coord in rand_rep_pos if chr_id == chromosome]
        #Define where to slice the shuffled_repeats list
        shuffled_end_index = shuffled_start_index + len(rand_rep_pos_filter)
        #Extract the slice of the shuffled_repeats for the corresponding chromosome
        shuffled_repeats_sliced = shuffled_repeats[shuffled_start_index:shuffled_end_index]
        #The end of the slice becomes the start for next chromosome
        shuffled_start_index = shuffled_end_index
        #Call generate_sequence function chromosome by chromosome
        sequence, new_repeats_coord = generate_sequence(repeats_dict, rand_rep_pos_filter, rand_seq, shuffled_repeats_sliced, alpha, beta, mode)
        genome_dict[chromosome] = sequence
        new_repeats_coord_dict[chromosome] = new_repeats_coord 
    #Return sequences and repeat data
    return genome_dict, new_repeats_coord_dict

#Print final genome sequence to file (function created by THC)
def print_genome_data(genome_dict, new_repeats_coord_dict, params, final_out):
    #Setup output directory   
    file_prefix = str(params['prefix'])
    os.chdir(final_out)
    
    #Output files
    genome_fa = file_prefix + "_genome_sequence_out.fasta"
    te_fa = file_prefix + "_repeat_sequence_out.fasta"
    te_gff = file_prefix + "_repeat_annotation_out.gff"
    
    #Create genome fasta file
    fasta_out = open(genome_fa, "w")
    for chromosome in genome_dict.keys():
        seq = str(genome_dict[chromosome])
        fasta_out.write(">" + chromosome + "\n" + seq + "\n")
    fasta_out.close()
    
    #Collapse new_repeats_coord_dict to a list with chr information
    new_repeats_coord_list = []
    for chromosome, replists in new_repeats_coord_dict.items():
        for replist in replists:
            new_repeats_coord_list.append([chromosome] + replist)   
    
    #Create repeat fasta and gff files
    fasta_rep = open(te_fa, "w")
    gff_rep = open(te_gff, "w")
    counts = 1
    for n in new_repeats_coord_list:
        chr_id = str(n[0])
        start = str(n[1])
        end = str(n[2])
        repeat_identity = str(n[4]/100)
        integ = str(n[5]/100)
        strand = str(n[6])
        te_family = str(n[7])
        te_subclass = str(n[8])
        te_superfamily = str(n[9])
        te_id = str(counts).zfill(7) #prints at least 7 characters wide; i.e. at most 9,999,999 TE insertions
        repeat_name = ">" + te_family + "_TE" +  te_id + "#" + te_superfamily + " [Location=" + chr_id + ":" + start + "-" + end + ";Identity=" + repeat_identity + ";Integrity=" + integ + "]\n"
        repeat_sequence = str(n[3])+ "\n"

        fasta_rep.write(repeat_name)
        fasta_rep.write(repeat_sequence)

        gff_rep.write("\t".join([chr_id, "TEgenomeSimulator", te_subclass, start, end, ".", strand, ".", "ID=" + te_family + "_TE" + te_id + ";Name=TE" + te_id + ";Classification=" + te_superfamily + ";Identity=" + repeat_identity + ";Integrity=" + integ +"\n"]))
        
        counts += 1
    fasta_rep.close()
    gff_rep.close()
    
def main():
    # For random or custom genome with multiple chromosomes
    # Set up argument parser
    parser = argparse.ArgumentParser(description="arguments of simulating random TE insertions.")

    # Define arguments
    parser.add_argument('-M', '--mode', type=int, help="Mode for genome simulation (either 0 or 1).", required=True)
    parser.add_argument('-p', '--prefix', type=str, help="Prefix for output files.", required=True)
    parser.add_argument('-a', '--alpha', type=float, default=0.5, help="Alpha value for the beta distribution used for fragmentation simulation.")
    parser.add_argument('-b', '--beta', type=float, default=0.7, help="Beta value for the beta distribution used for fragmentation simulation.")
    parser.add_argument('-o', '--outdir', type=str, help="Output directory.", required=True)
     
    # Parse arguments
    args = parser.parse_args()
    mode = args.mode
    prefix = args.prefix
    alpha = args.alpha 
    beta = args.beta
    final_out = args.outdir
    
    print("\n", flush=True)
    print("##############################################################", flush=True)
    print("### Mutate TE sequence and perform non-overlap TE insertion###", flush=True)
    print("##############################################################", flush=True)
    print(f"Using mode {mode} (0 for random genome; 1 for custome genome)", flush=True)

    # Config file
    yml_file = "TEgenomeSimulator_" + str(prefix) + ".yml"
    print(f"Using config file {yml_file}", flush=True)

    # Mode-dependent config file loading
    if args.mode == 0:
        #Load chr parameters from yml file using parse_random_genome_yaml()        
        params_chr = parse_random_genome_yaml(os.path.join(final_out, yml_file))
    elif args.mode == 1 or args.mode == 2:
        #Load chr parameters from yml file using parse_custom_genome_yaml()
        params_chr = parse_custom_genome_yaml(os.path.join(final_out, yml_file))

    #Set seed
    seed = params_chr['seed']
    if seed:
        random.seed(seed)
        numpy.random.seed(seed)

    # Prepare genome for random or custom mode
    if args.mode == 0:
        chrs_dict = generate_random_chr_sequence(params_chr)
    elif args.mode == 1 or args.mode == 2:
        chrs_dict = load_custom_genome(params_chr)

    #Load repeat sequences
    if args.mode == 0 or args.mode == 1:
        repeats_dict = load_repeats_chr(params_chr)
    elif args.mode == 2:
        repeats_dict = load_repeats_chr_m2(params_chr)  
        #repeats_dict = load_repeats_chr(params_chr)  

    #Assign TE coordinates randomly
    repeats_coord = assign_chr_coord_repeats(params_chr, repeats_dict)
    #Shuffle the order the repeats to be inserted into the genome
    shuffled_repeats = shuffle_repeats(repeats_dict)
    #Generate genome sequence with TE insertion
    genome, new_repeats_coord = generate_genome_sequence(repeats_dict, repeats_coord, chrs_dict, shuffled_repeats, alpha, beta, mode)
    #Output to fasta and gff files
    print_genome_data(genome, new_repeats_coord, params_chr, final_out)
if __name__ == "__main__":
    main()
