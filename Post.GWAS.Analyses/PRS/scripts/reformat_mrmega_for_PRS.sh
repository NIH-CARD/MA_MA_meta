#!/bin/bash
# sh reformat_mrmega_for_PRS.sh $FILENAME $PREFIX
FILENAME=$1
PREFIX=$2
cat header.txt > ${PREFIX}.forMungeSumstats.txt
awk 'OFS="\t" {if (NR!=1) print $2, $3, $1, $4, $5, $11, $7}' $FILENAME >> ${PREFIX}.forMungeSumstats.txt 
