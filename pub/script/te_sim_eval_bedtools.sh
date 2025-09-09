#!/bin/bash
#!/bin/bash

# Inputs
rm_out=$1
truth_gff=$2
cutoff=$3 # minimum reciprocal overlap

# 1. Convert RepeatMasker .out to BED, ensure start < end
awk 'NR>3 {
    start=$6-1; end=$7;
    if (start > end) {tmp=start; start=end; end=tmp}
    print $5 "\t" start "\t" end "\t" $10
}' "$rm_out" > repeatmasker.bed

# 2. Convert truth GFF to BED, ensure start < end
awk '{
    start=$4-1; end=$5;
    if (start > end) {tmp=start; start=end; end=tmp}
    print $1 "\t" start "\t" end "\t" $9
}' "$truth_gff" > truth.bed

# 3. Intersect: with overlap and without overlap
bedtools intersect -u -f $cutoff -r -a truth.bed -b repeatmasker.bed > truth_with_overlap.bed
bedtools intersect -v -f $cutoff -r -a truth.bed -b repeatmasker.bed > truth_without_overlap.bed


# 4. Extract integrity values
extract_integrity () {
  awk -v cat=$1 '{
    integ=""; ident="";
    if ($4 ~ /Integrity=/) {
      match($4,/Integrity=([0-9.]+)/,a); integ=a[1];
    }
    if ($4 ~ /Identity=/) {
      match($4,/Identity=([0-9.]+)/,b); ident=b[1];
    }
    if (integ != "" && ident != "") {
      print cat, integ, ident;
    }
  }' $2
}

    # Run for both categories
echo -e "category\tintegrity\tidentity" > with_overlap_integrity.tsv
extract_integrity with_overlap truth_with_overlap.bed >> with_overlap_integrity.tsv

echo -e "category\tintegrity\tidentity" > without_overlap_integrity.tsv
extract_integrity without_overlap truth_without_overlap.bed >> without_overlap_integrity.tsv

# extract_integrity () {
#   awk -v cat=$1 '{
#     match($4,/Integrity=([0-9.]+)/,a);
#     if (a[1]!="") {
#       bin = int(a[1]*10); if (bin==10) bin=9;
#       print cat, a[1], bin;
#     }
#   }' $2
# }
# 
# extract_integrity with_overlap truth_with_overlap.bed > with_overlap_integrity.tsv
# extract_integrity without_overlap truth_without_overlap.bed > without_overlap_integrity.tsv

# 5. Summarize per bin with labeled intervals
cat with_overlap_integrity.tsv without_overlap_integrity.tsv | \
  awk '
  {
    integ_bin = int($2*10); 
    if (integ_bin==10) integ_bin=9;
    div_bin   = int($3*10);
    if (div_bin==10) div_bin=9;

    key = integ_bin "_" div_bin;
    total[key]++;
    if ($1=="with_overlap") recovered[key]++;
  }
  END{
    print "integrity_interval","divergence_interval","total","recovered","rate";
    for (i=0;i<10;i++) {
      low_i = i/10; high_i = (i+1)/10;
      if (i<9) integ_label = sprintf("[%.1f-%.1f)", low_i, high_i);
      else     integ_label = sprintf("[%.1f-%.1f]", low_i, high_i);
      for (j=0;j<10;j++) {
        low_j = j/10; high_j = (j+1)/10;
        if (j<9) div_label = sprintf("[%.1f-%.1f)", low_j, high_j);
        else     div_label = sprintf("[%.1f-%.1f]", low_j, high_j);
        key = i "_" j;
        tot = (key in total ? total[key] : 0);
        rec = (key in recovered ? recovered[key] : 0);
        rate = (tot>0 ? rec/tot : 0);
        print integ_label, div_label, tot, rec, rate;
      }
    }
  }' OFS="\t" > recovery_by_bin2D.tsv






#cat with_overlap_integrity.tsv without_overlap_integrity.tsv | \
#  awk '
#  {
#    bin = int($3); total[bin]++; if($1=="with_overlap") recovered[bin]++
#  }
#  END{
#    print "integrity_interval","total","recovered","rate";
#    for (i=0;i<10;i++) {
#      low = i/10; high = (i+1)/10;
#      if (i<9) interval = sprintf("[%.1f-%.1f)", low, high);
#      else interval = sprintf("[%.1f-%.1f]", low, high);  # last bin closed
#      rate = (total[i]>0 ? recovered[i]/total[i] : 0);
#      print interval, total[i]+0, recovered[i]+0, rate;
#    }
#  }' OFS="\t" > recovery_by_bin.tsv


## 4. Extract Identity & Integrity
#extract_fields() {
#    file=$1
#    category=$2
#    awk -v cat="$category" '{
#        id = ""; integ = "";
#        if ($4 ~ /Identity=/) {
#            match($4, /Identity=([^;]+)/, a); id=a[1];
#        }
#        if ($4 ~ /Integrity=/) {
#            match($4, /Integrity=([^;]+)/, b); integ=b[1];
#        }
#        if (id != "" && integ != "") {
#            print cat "\t" id "\t" integ
#        }
#    }' "$file"
#}

# 6. Calculate numbers
total=$(wc -l < truth.bed)
with_overlap=$(wc -l < truth_with_overlap.bed)
no_overlap=$(( total - with_overlap ))
ratio=$(echo "$no_overlap / $total" | bc -l)

## 6. Output identity and integrity
#echo -e "category\tidentity\tintegrity" > overlap_eval_vs_seq_decay.txt
#extract_fields truth_with_overlap.bed with_overlap >> overlap_eval_vs_seq_decay.txt
#extract_fields truth_without_overlap.bed without_overlap >> overlap_eval_vs_seq_decay.txt

# 7. Print tab-delimited results
echo -e "Total_truth\tTruth_with_overlap\tTruth_without_overlap\tRatio_non_overlap" > overlap_eval_stat.txt
echo -e "$total\t$with_overlap\t$no_overlap\t$ratio" >> overlap_eval_stat.txt

