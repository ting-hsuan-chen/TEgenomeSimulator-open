#!/bin/bash

#Ting-Hsuan Chen
#2024-03-28

################
##### dm6 ######
################
WRK=/workspace/cflthc/script/KRIP_TE/10_TEgenomeSimulator
REF=$WRK/model_species_ref
OUT=$WRK/model_species_RM/dm6
RMaskerOut=$OUT/"dm6.fa.out"
RMaskerTbl=$OUT/"dm6.fa.tbl"
libindex=$REF/"Dfam_curatedonly.fasta.fai"
liblength=$OUT/"Dfam_curatedonly.length"

# extract genome size and masked bases
#genomesize=143726002
#maskedsize=27189125
genomesize=$(awk -F' ' 'NR==4{print $3}' $RMaskerTbl)
maskedsize=$(awk -F' ' 'NR==6{print $3}' $RMaskerTbl)

cd $OUT

# extract length of each TE family
awk -F'\t' '{OFS="\t"; oldname = $1; gsub(/#.*/,"",$1); print $0, oldname}' $libindex | awk '{FS = OFS="\t"}{print $1, $2, $6}' | sort -k1,1 > $liblength

# forward strand of RepeatMasker output
awk 'NR>3 && $9 == "+" {print}' $RMaskerOut |\
sed -e 's/(//g' -e 's/)//g' |\
awk 'mismatch=($2+$3+$4) {OFS="\t"; print $10,$11,mismatch/100,(100-mismatch)/100,($3+$4)/mismatch,$7-$6+1,$12-1,$14}' \
> forward.temp

# reverse strand of RepeatMasker output
awk 'NR>3 && $9 == "C" {print}' $RMaskerOut |\
sed -e 's/(//g' -e 's/)//g' |\
awk 'mismatch=($2+$3+$4) {OFS="\t"; print $10,$11,mismatch/100,(100-mismatch)/100,($3+$4)/mismatch,$7-$6+1,$14-1,$12}' \
> reverse.temp

# combine forward and reverse
cat forward.temp reverse.temp | sort -k1,1 > combined.temp

# join matched family and attach full length as last column 
# calulate length_ratio and cut_ratio
# cut_ratio = (headcut + tailcut)/full_length
# length_ratio = lengh/full_length
# Unify superfamily nomenclature and add column names
awk -F'\t' 'BEGIN{OFS="\t"} {if(FNR==NR) {a[$1]=$2; next} if($1 in a){print $0, a[$1]}}' $liblength combined.temp |\
awk -F'\t' 'BEGIN{OFS="\t"} {if(FNR==NR) {a[$1]=$3; next} if($1 in a){print $0, a[$1]}}' $liblength - |\
awk -F'\t' '{OFS="\t"; print $0, $6/$9, ($7+$8)/$9}' |\
sed -e 's/\sDNA\/CMC-Chapaev/\tDNA\/DTC/g' -e 's/\sDNA\/CMC-EnSpm/\tDNA\/DTC/g' -e 's/\sDNA\/En-Spm/\tDNA\/DTC/g' -e 's/\sDNA\/hAT-Ac/\tDNA\/DTA/g' -e 's/\sDNA\/hAT-Charlie/\tDNA\/DTA/g' -e 's/\sDNA\/hAT-Tag1/\tDNA\/DTA/g' -e 's/\sDNA\/hAT-Tip100/\tDNA\/DTA/g' -e 's/\sDNA\/hAT/\tDNA\/DTA/g' -e 's/\sDNA\/HAT/\tDNA\/DTA/g' -e 's/\sDNA\/MULE-MuDR/\tDNA\/DTM/g' -e 's/\sDNA\/MuDR/\tDNA\/DTM/g' -e 's/\sDNA\/PIF-Harbinger/\tDNA\/DTH/g' -e 's/\sDNA\/Harbinger/\tDNA\/DTH/g' -e 's/\sDNA\/TcMar-Tc1/\tDNA\/DTT/g' -e 's/\sDNA\/TcMar-Pogo/\tDNA\/DTT/g' -e 's/\sDNA\/TcMar-Stowaway/\tDNA\/DTT/g' -e 's/\sDNA\/Mariner/\tDNA\/DTT/g' -e 's/\sDNA\/Pogo/\tDNA\/DTT/g' -e 's/\sDNA\/Tc1/\tDNA\/DTT/g' -e 's/\sLTR\//\tLTRx/' -e 's/\sunknown/\tUnknown/g' |\
sed -e 's/\sLTR\>/\tLTR\/LTR_Unknown/' -e 's/\sLTRx/\tLTR\//g' -e 's/\sLTR\/Unknown/\tLTR\/LTR_Unknown/g' -e 's/\sLINE\/Unknown/\tLINE\/LINE_Unknown/' -e 's/\sLINE?/\tLINE\/LINE_Unknown/' -e 's/\sSINE\>/\tSINE\/SINE_Unknown/' |\
awk '{FS = OFS = "\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $11, $12, $10}' |\
sed 1i"family\tsuperfamily\tmismatch\tidentity\tindel\tlength\theadcut\ttailcut\tfamily_length\tlength_ratio\tcut_ratio\toriginal_family" > $RMaskerOut.processed

# summarise each family
# TE loci shorter than 90% of the TE family sequence are considered fragmented (fragment_copy)
# fragment_ratio is calculted as fragment_copy/copynumber
# standard deviation of identity was calculated as sqrt(((identity - (identity/copy_number))^2)/copy_number)
awk 'NR>1' $RMaskerOut.processed | sort -u -k2,2 -k1,1 | awk '{OFS=","; print $1, $2}' > $RMaskerOut.familylist
rm -f summary.temp

cat $RMaskerOut.familylist | while read line; do
family=$(echo $line | sed 's/,.*//g')
superfamily=$(echo $line | sed 's/.*,//g')
awk -v fam=${family} '{FS = OFS = "\t"}{if($1==fam){print $0}}' $RMaskerOut.processed > family.temp
#grep -e "^${family}\s" $RMaskerOut.processed > family.temp
frag_c=$(awk '$10 < 0.9 { count++ } END { print count }' 'family.temp')
awk -v fam=$family -v sfam=$superfamily -v frag_copy=$frag_c '{FS=OFS="\t"; fulllen=$9; idt+=$4; idl+=$5; len+=$6; headcut+=$7; tailcut+=$8; cut+=($7+$8)} END {Z+=($4-(idt/NR))^2; print fam, sfam, fulllen, NR, idt/NR, sqrt(Z/NR), idl/NR, len, len/NR, headcut/NR, tailcut/NR, cut/NR, frag_copy, frag_copy/NR, $12}' family.temp >> summary.temp
done


# Unify superfamily nomenclature and add column names
#sed -e 's/CMC-Chapaev/DTC/g' -e 's/CMC-EnSpm/DTC/g' -e 's/En-Spm/DTC/g' -e 's/hAT-Ac/DTA/g' -e 's/hAT-Charlie/DTA/g' -e 's/hAT-Tag1/DTA/g' -e 's/hAT-Tip100/DTA/g' -e 's/hAT/DTA/g' -e 's/HAT/DTA/g' -e 's/MULE-MuDR/DTM/g' -e 's/MuDR/DTM/g' -e 's/PIF-Harbinger/DTH/g' -e 's/Harbinger/DTH/g' -e 's/TcMar-Pogo/DTT/g' -e 's/TcMar-Tc1/DTT/g' -e 's/TcMar-Pogo/DTT/g' -e 's/TcMar-Stowaway/DTT/g' -e 's/Mariner/DTT/g' -e 's/Pogo/DTT/g' -e 's/Tc1/DTT/g' -e 's/LTR\//LTRx/' -e 's/RC\/Helitron/DNA\/Helitron/g' -e 's/unknown/Unknown/g' -e 's/Unspecified/Unknown/g' summary.temp |\
#sed -e 's/LTR\>/LTR\/LTR_Unknown/' -e 's/LTRx/LTR\//g' -e 's/LTR\/Unknown/LTR\/LTR_Unknown/g' -e 's/LINE\/Unknown/LINE\/LINE_Unknown/' -e 's/LINE?/LINE\/LINE_Unknown/' -e 's/SINE\>/SINE\/SINE_Unknown/' |\
#sed 1i"family\tsuperfamily\tfull_length_bp\tcopynumber\tmean_identity\tmean_indel\ttotal_length\tmean_length\tmean_headcut\tmean_tailcut\tmean_cut\tfragment_copy\tfragment_ratio\toriginal_family" > $RMaskerOut.summary

sed 1i"family\tsuperfamily\tfull_length_bp\tcopynumber\tmean_identity\tsd_identity\tmean_indel\ttotal_length\tmean_length\tmean_headcut\tmean_tailcut\tmean_cut\tfragment_copy\tfragment_ratio\toriginal_family" summary.temp > $RMaskerOut.summary

# create high level summary
awk 'NR>1{print $2}' $RMaskerOut.summary | sort | uniq > $RMaskerOut.superfamilylist
rm -f sfam_summary.temp
cat $RMaskerOut.superfamilylist | while read superfamily; do
grep -e "${superfamily}" $RMaskerOut.summary > superfamily.temp
awk -v sfam=$superfamily -v gsize=$genomesize -v msize=$maskedsize '{FS=OFS="\t"; loci_sum+=$4; len_sum+=$8; idt_sum+=$5} END {print sfam, NR, loci_sum, len_sum, idt_sum/NR, len_sum*100/gsize, len_sum*100/msize}' superfamily.temp >> sfam_summary.temp
done

#awk 'BEGIN {OFS="\t"; print "superfamily", "family_count", "total_loci", "total_length", "mean_identity", "percentage_in_genome", "percentage_in_masked"}{ print $0, "" }' sfam_summary.temp > $RMaskerOut.highlevel_summary

sed 1i"superfamily\tfamily_count\ttotal_loci\ttotal_length\tmean_identity\tpercentage_in_genome\tpercentage_in_masked" sfam_summary.temp > $RMaskerOut.highlevel_summary

rm -f *.temp

