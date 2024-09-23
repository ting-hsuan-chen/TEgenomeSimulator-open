#!/bin/bash

#Ting-Hsuan Chen
#2023-12-21
WRK=/workspace/cflthc/scratch/2022_Actinidia_TE
OUT=$WRK/09.10_BM_modules/curated_TElib
cd $OUT

comb_TElib=combined_curated_TE_lib_ATOSZM_cdhit.fasta
outprefix="combined_curated_TE_lib_ATOSZM_cdhit"

# Convert multi-line fasta to single line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $comb_TElib > temp.fasta

# LTR Copia
grep -A1 "Copia" temp.fasta > $outprefix'_LTRcopia.fasta'
# LTR Gypsy
grep -A1 -e "Gypsy" -e "Ty3" temp.fasta > $outprefix'_LTRgypsy.fasta'
# LTR unknown
grep -A1 "LTR" temp.fasta | grep -A1 "unknown" > $outprefix'_LTRunknown.fasta'
# LTR Solo
grep -A1 "LTR" temp.fasta | grep -A1 "Solo" > $outprefix'_LTRsolo.fasta'
# LINE unknown + L1
grep -A1 "LINE" temp.fasta > $outprefix'_LINE.fasta'
# SINE unknown
grep -A1 "SINE" temp.fasta > $outprefix'_SINE.fasta'
# DNA CACTA
grep -A1 "\/CACTA$" temp.fasta > $outprefix'_DNAcacta.fasta'
# DNA CACTG
grep -A1 "\/CACTG$" temp.fasta > $outprefix'_DNAcactg.fasta'
# DNA hAT
grep -A1 "hAT" temp.fasta > $outprefix'_DNAhat.fasta'
# DNA MuDR (also need to remove entries that only contain tir sequence)
grep -A1 -e "MULE" -e "MuDR" temp.fasta | tr '\n' '@' | sed -e 's/@>/\n>/g' -e 's/@--//g' | grep -v "MULEtir" | sed 's/@/\n/g' | head -n -1 > $outprefix'_DNAmudr.fasta'
# DNA Harbinger
grep -A1 "Harbinger" temp.fasta > $outprefix'_DNAharbinger.fasta'
# DNA Mariner
grep -A1 "Mariner" temp.fasta > $outprefix'_DNAmariner.fasta'
# DNA POLE
grep -A1 "POLE" temp.fasta > $outprefix'_DNApole.fasta'
# DNA PILE
grep -A1 "PILE" temp.fasta > $outprefix'_DNApile.fasta'
# DNA Tourist
grep -A1 "Tourist" temp.fasta > $outprefix'_DNAtourist.fasta'
# DNA Spacer
grep -A1 "spacer" temp.fasta > $outprefix'_DNAspacer.fasta'
# DNA PILE
grep -A1 "PILE" temp.fasta > $outprefix'_DNApile.fasta'
# RC Helitron
grep -A1 "Helitron" temp.fasta > $outprefix'_RChelitron.fasta'
# MITE Tourist
grep -A1 "MITE" temp.fasta | grep -A1 "Tourist" > $outprefix'_MITEtourist.fasta'
# MITE Stow
grep -A1 "MITE" temp.fasta | grep -A1 "Stow" > $outprefix'_MITEstow.fasta'

rm -f temp.fasta

ml samtools/1.16
for file in *_cdhit_*; do
samtools faidx $file
done
