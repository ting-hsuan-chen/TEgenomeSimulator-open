#!/bin/bash

#Ting-Hsuan Chen
#2023-12-21
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/09.05_RepeatMasker
RMaskerOut="Donghong.chromosomes.only.fa.out"
RMaskerTbl="Donghong.chromosomes.only.fa.tbl"
libindex="Donghong_combined_lib_dedup_uniqname.fa.fai"
liblength="Donghong_combined_lib_dedup_uniqname.length"

# extract genome size and masked bases
#genomesize=608327852
#maskedsize=324908155
genomesize=$(awk -F' ' 'NR==4{print $3}' $RMaskerTbl)
maskedsize=$(awk -F' ' 'NR==6{print $3}' $RMaskerTbl)

cd $OUT

# extract length of each TE family
awk -F'\t' '{OFS=","; gsub(/#/,",",$1);print}' $libindex | awk -F, '{OFS="\t"; print $1, $3}' | sort -k1,1 > $liblength

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
awk -F'\t' 'BEGIN{OFS="\t"} {if(FNR==NR) {a[$1]=$2; next} if($1 in a){print $0, a[$1]}}' $liblength combined.temp |\
awk -F'\t' '{OFS="\t"; print $0, $6/$9, ($7+$8)/$9}' |\
awk 'BEGIN {OFS="\t"; print "#family", "superfamily", "mismatch", "identity", "indel", "length", "headcut", "tailcut", "family_length", "length_ratio", "cut_ratio"}{ print $0, "" }' > $RMaskerOut.processed

# summarise each family
# TE loci shorter than 90% of the TE family sequence are considered fragmented (fragment_copy)
# fragment_ratio is calculted as fragment_copy/copynumber
awk 'NR>1' $RMaskerOut.processed | sort -u -k2,2 -k1,1 | awk '{OFS=","; print $1, $2}' > $RMaskerOut.familylist
rm -f summary.temp
cat $RMaskerOut.familylist | while read line; do
family=$(echo $line | sed 's/,.*//g')
superfamily=$(echo $line | sed 's/.*,//g')
grep -e "${family}" $RMaskerOut.processed > family.temp
frag_c=$(awk '$10 < 0.9 { count++ } END { print count }' 'family.temp')
awk -v fam=$family -v sfam=$superfamily -v frag_copy=$frag_c '{FS=OFS="\t"; fulllen=$9; idt+=$4; idl+=$5; len+=$6; headcut+=$7; tailcut+=$8; cut+=($7+$8)} END { print fam, sfam, fulllen, NR, idt/NR, idl/NR, len, len/NR, headcut/NR, tailcut/NR, cut/NR, frag_copy, frag_copy/NR}' family.temp >> summary.temp
done

# Unify superfamily nomenclature and add column names
sed -e 's/CMC-Chapaev/DTC/g' -e 's/CMC-EnSpm/DTC/g' -e 's/hAT-Ac/DTA/g' -e 's/hAT-Tag1/DTA/g' -e 's/hAT-Tip100/DTA/g' -e 's/MULE-MuDR/DTM/g' -e 's/PIF-Harbinger/DTH/g' -e 's/TcMar-Pogo/DTT/g' -e 's/TcMar-Tc1/DTT/g' -e 's/LTR\//LTRx/' -e 's/RC\/Helitron/DNA\/Helitron/g' -e 's/unknown/Unknown/g' summary.temp |\
sed -e 's/LTR\>/LTR\/LTR_Unknown/' -e 's/LTRx/LTR\//g' -e 's/LTR\/Unknown/LTR\/LTR_Unknown/g' -e 's/LINE\/Unknown/LINE\/LINE_Unknown/' |\
awk 'BEGIN {OFS="\t"; print "#family", "superfamily", "full_length_bp", "copynumber", "mean_identity", "mean_indel", "total_length", "mean_length", "mean_headcut", "mean_tailcut", "mean_cut", "fragment_copy", "fragment_ratio"}{ print $0, "" }' > $RMaskerOut.summary

# create high level summary
awk 'NR>1{print $2}' $RMaskerOut.summary | sort | uniq > $RMaskerOut.superfamilylist
rm -f sfam_summary.temp
cat $RMaskerOut.superfamilylist | while read superfamily; do
grep -e "${superfamily}" $RMaskerOut.summary > superfamily.temp
awk -v sfam=$superfamily -v gsize=$genomesize -v msize=$maskedsize '{FS=OFS="\t"; loci_sum+=$4; len_sum+=$7; idt_sum+=$5} END {print sfam, NR, loci_sum, len_sum, idt_sum/NR, len_sum*100/gsize, len_sum*100/msize}' superfamily.temp >> sfam_summary.temp
done

awk 'BEGIN {OFS="\t"; print "#superfamily", "family_count", "total_loci", "total_length", "mean_identity", "percentage_in_genome", "percentage_in_masked"}{ print $0, "" }' sfam_summary.temp > $RMaskerOut.highlevel_summary

rm -f *.temp