#! user/bin/env bash

histone_bed="../data-sets/bed/encode.h3k4me3.hela.chr22.bed.gz"
CTCF="../data-sets/bed/encode.tfbs.chr22.bed.gz"
gzcat $CTCF | awk '($4=="CTCF")'>ctcf-peaks.bed

answer_1=$(bedtools intersect -a ctcf-peaks.bed -b $histone_bed -wo \
    |awk '{print $NF}' \
    |sort -nr |head -n1)
     echo "answer-1:$answer_1"
     echo "answer-1:$answer_1"> answers.yml
     
 #Q2 Calculate GC content

chr22fa="../data-sets/fasta/hg19.chr22.fa"
echo -e "chr22\t19000000\t19000500">region.bed
answer_2=$(bedtools nuc -fi $chr22fa -bed region.bed|grep -v '#'| cut -f5)
echo "answer-2:$answer_2"
echo "answer-2:$answer_2">>answers.yml

#Q3 Calculate the length of CTCF chip seq (largest mean signal)

CTHELA="../data-sets/bedtools/ctcf.hela.chr22.bg"
answer_3=$(bedtools map -a ctcf-peaks.bed -b $CTHELA -c 4 -o mean \
    |sort -k5n \
    |tail -n1 \
    |awk '{print $3-$2}') 

echo "answer-3:$answer_3"
echo "answer-3:$answer_3">>answers.yml

#Q4
genome="../data-sets/genome/hg19.genome"
tss="../data-sets/bed/tss.hg19.chr22.bed"

(bedtools flank -l 1000 -r 0 -s -i $tss -g $genome)>tssflanks.bed 
(bedtools sort -i tssflanks.bed)>tmp.bed 

answer_4=$(bedtools map -a tmp.bed -b $CTHELA -c 4 -o median \
    |sort -k7n \
    |tail -n1 \
    |awk '{print $4}')
echo "answer-4:$answer_4"
echo "answer-4:$answer_4">>answers.yml

# Q5 longest section of chromosome that's not covered by a gene
genesbed="../data-sets/bed/genes.hg19b.bed"

answer_5=$(bedtools complement -i $genesbed -g $genome \
    |awk '{OFS="\t"}{print $3-$2, $1, $2, $3}' \
    |sort -k1n \
    |tail -n1 \
    |awk '{print $2":"$3"-"$4}')

echo "answer-5:$answer_5"
echo "answer-5:$answer_5">>answers.yml
