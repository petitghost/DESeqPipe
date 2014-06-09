#!/bin/sh
##################################################
#Date:2014-05-21
#Author: ghost
#Usage: HTseq_DEseq_pipeline.sh [gtf file] [sorting bam of sample 1] [sorting bam of sample 2]
#Description: (Step 1) Counting read on annotation feature using HTseq; (Step 2) Find differential expression genes using DESeq
##################################################

#Sample 1
python -m HTSeq.scripts.count -f bam -t CDS $2 $1 > readCount.s1 2> ht_err.s1 &
#Sample 2
python -m HTSeq.scripts.count -f bam -t CDS $3 $1 > readCount.s2 2> ht_err.s2 &


wait
#Merge read count of each gene on 2 samples
awk 'NR==FNR{if($1!~/^__/){a[$1]=$2};next}BEGIN{print "GeneID\tSample1\tSample2"}{for(i in a){if(i==$1){print $0"\t"a[i];break}}}' readCount.s1 readCount.s2 > s1_s2.readCount
R -f ~/bin/DEseq.R 
