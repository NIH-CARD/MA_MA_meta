#!/bin/bash
# Goal: Run MR-MEGA with max PCs for 5 datasets (PCs = 2) or 4 datasets (PCs = 1).
# PC count must be < cohort count - 2
# sh run_MRMEGA_max_PCs.sh $IN_FILE
IN_FILE=$1
NUM_DATASETS=$(cat $IN_FILE | wc -l)
NUM_PCS=$(expr $NUM_DATASETS - 3)
/path/to/MR-MEGA_v0.2/MR-MEGA -i $IN_FILE --pc $NUM_PCS -o multi-ancestry_${IN_FILE}.MAX_PCs
