import sys
import statistics

def extract_sizes(rmasker_tbl):
    with open(rmasker_tbl, 'r') as file:
        lines = file.readlines()
        
        genome_size = int(lines[3].split()[2]) if len(lines) > 3 else None
        masked_size = int(lines[5].split()[2]) if len(lines) > 5 else None
        
        return genome_size, masked_size

def extract_te_family_length(libindex):
    te_data = []
    with open(libindex, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                oldname = parts[0]
                clean_name = oldname.split('#')[0]  # Remove everything after #
                te_data.append((clean_name, parts[1], oldname))
    
    # Sort by the first column
    te_data.sort(key=lambda x: x[0])
    
    # Write to output file
    liblength = "liblength.temp"
    with open(liblength, 'w') as file:
        for entry in te_data:
            file.write('\t'.join(entry) + '\n')
    
    return liblength

def process_repeatmasker_output(rmasker_out):
    forward_temp = "forward.temp"
    reverse_temp = "reverse.temp"
    combined_temp = "combined.temp"
    forward_data = []
    reverse_data = []
    
    with open(rmasker_out, 'r') as file:
        lines = file.readlines()[3:]  # Skip first 3 lines
        for line in lines:
            line = line.replace('(', '').replace(')', '')  # Remove ( and ) characters
            parts = line.strip().split()
            if len(parts) < 15:
                continue
            
            mismatch = float(parts[1]) + float(parts[2]) + float(parts[3])
            mismatch_fraction = mismatch / 100
            match_fraction = (100 - mismatch) / 100
            indel_fraction = (float(parts[2]) + float(parts[3])) / mismatch if mismatch else 0
            length = int(parts[6]) - int(parts[5]) + 1
            
            if parts[8] == "+":
                forward_data.append((parts[9], parts[10], mismatch_fraction, match_fraction, indel_fraction, length, int(parts[11]) - 1, parts[13]))
            elif parts[8] == "C":
                reverse_data.append((parts[9], parts[10], mismatch_fraction, match_fraction, indel_fraction, length, int(parts[13]) - 1, parts[11]))
    
    # Write forward strand
    with open(forward_temp, 'w') as file:
        for entry in forward_data:
            file.write('\t'.join(map(str, entry)) + '\n')
    
    # Write reverse strand
    with open(reverse_temp, 'w') as file:
        for entry in reverse_data:
            file.write('\t'.join(map(str, entry)) + '\n')
    
    # Combine and sort
    combined_data = sorted(forward_data + reverse_data, key=lambda x: x[0])
    with open(combined_temp, 'w') as file:
        for entry in combined_data:
            file.write('\t'.join(map(str, entry)) + '\n')

    return combined_temp

def process_final_output(liblength, combined_temp):
    processed_output = "processed_output.temp"
    family_data = {}
    with open(liblength, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                family_data[parts[0]] = (parts[1], parts[2])
    
    processed_data = []
    with open(combined_temp, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if parts[0] in family_data:
                family_length = int(family_data[parts[0]][0])
                length_ratio = float(parts[5]) / family_length if family_length else 0
                cut_ratio = (int(parts[6]) + int(parts[7])) / family_length if family_length else 0
                processed_data.append(parts + [str(family_length), str(length_ratio), str(cut_ratio), family_data[parts[0]][0]])
    
    # Write processed output
    with open(processed_output, 'w') as file:
        file.write("family\tsuperfamily\tmismatch\tidentity\tindel\tlength\theadcut\ttailcut\tfamily_length\tlength_ratio\tcut_ratio\toriginal_family\n")
        for entry in processed_data:
            file.write('\t'.join(entry) + '\n')

    return processed_output

def summarize_each_family(processed_output):
    family_list = "family_list.temp"
    summary_fam = "summarise_repeatmasker_out_family.txt"
    
    families = set()
    with open(processed_output, 'r') as file:
        next(file)
        for line in file:
            parts = line.strip().split('\t')
            families.add((parts[0], parts[1]))
    
    with open(family_list, 'w') as file:
        for fam, sfam in sorted(families):
            file.write(f"{fam},{sfam}\n")
    
    with open(summary_fam, 'w') as summary_file:
        summary_file.write("family\tsuperfamily\tfull_length_bp\tcopynumber\tmean_identity\tsd_identity\tmean_indel\ttotal_length\tmean_length\tmean_headcut\tmean_tailcut\tmean_cut\tfragment_copy\tfragment_ratio\toriginal_family\n")
        
        for fam, sfam in sorted(families):
            frag_count = 0
            total_identity = total_indel = total_length = total_headcut = total_tailcut = total_cut = count = 0
            identities = []
            
            with open(processed_output, 'r') as file:
                next(file)
                for line in file:
                    parts = line.strip().split('\t')
                    if parts[0] == fam:
                        count += 1
                        identities.append(float(parts[4]))
                        total_identity += float(parts[4])
                        total_indel += float(parts[5])
                        total_length += int(parts[6])
                        total_headcut += int(parts[7])
                        total_tailcut += int(parts[8])
                        total_cut += (int(parts[7]) + int(parts[8]))
                        if float(parts[9]) < 0.9:
                            frag_count += 1
            
            mean_identity = total_identity / count if count else 0
            sd_identity = statistics.stdev(identities) if len(identities) > 1 else 0
            mean_indel = total_indel / count if count else 0
            mean_length = total_length / count if count else 0
            mean_headcut = total_headcut / count if count else 0
            mean_tailcut = total_tailcut / count if count else 0
            mean_cut = total_cut / count if count else 0
            fragment_ratio = frag_count / count if count else 0
            
            summary_file.write(f"{fam}\t{sfam}\t{total_length}\t{count}\t{mean_identity}\t{sd_identity}\t{mean_indel}\t{total_length}\t{mean_length}\t{mean_headcut}\t{mean_tailcut}\t{mean_cut}\t{frag_count}\t{fragment_ratio}\t{fam}\n")
    
    return summary_fam


def create_high_level_summary(summary_fam, genome_size, masked_size):
    highlevel_summary = "summarise_repeatmasker_out_superfamily.txt"
    superfamilies = {}
    
    with open(summary_fam, 'r') as file:
        next(file)
        for line in file:
            parts = line.strip().split('\t')
            superfamily = parts[1]
            loci = int(parts[3])
            length = int(parts[7])
            identity = float(parts[4])
            
            if superfamily not in superfamilies:
                superfamilies[superfamily] = [0, 0, 0, 0]
            superfamilies[superfamily][0] += 1
            superfamilies[superfamily][1] += loci
            superfamilies[superfamily][2] += length
            superfamilies[superfamily][3] += identity
    
    with open(highlevel_summary, 'w') as file:
        file.write("superfamily\tfamily_count\ttotal_loci\ttotal_length\tmean_identity\tpercentage_in_genome\tpercentage_in_masked\n")
        for sfam, values in superfamilies.items():
            mean_identity = values[3] / values[0] if values[0] else 0
            perc_genome = (values[2] * 100) / genome_size if genome_size else 0
            perc_masked = (values[2] * 100) / masked_size if masked_size else 0
            file.write(f"{sfam}\t{values[0]}\t{values[1]}\t{values[2]}\t{mean_identity}\t{perc_genome}\t{perc_masked}\n")
    
    return highlevel_summary


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python summarise_rm_out.py <RMaskerTbl> <LibIndex> <RMaskerOut>")
        sys.exit(1)
    
    rmasker_tbl = sys.argv[1]
    libindex = sys.argv[2]
    rmasker_out = sys.argv[3]
    
    genome_size, masked_size = extract_sizes(rmasker_tbl)
    liblength = extract_te_family_length(libindex)
    combined_temp = process_repeatmasker_output(rmasker_out)
    processed_output = process_final_output(liblength, combined_temp)
    summary_output = summarize_each_family(processed_output)
    create_high_level_summary(summary_output, genome_size, masked_size)
