#!/bin/bash
# Goal: To increase the variant set output by MR-MEGA, run each combination of four input GWAS with a single PC.
# Then merge the outputs to maximize the effective sample size for each variant. 

# sh merge_MRMEGA_outputs.sh $OUTPUT_ALL_DATASETS $LEAVE_ONE_OUT_FILES
OUTPUT_ALL_DATASETS=$1
LEAVE_ONE_OUT_FILES=$2

# Make some temp files
tmp_initial=$(mktemp)
tmp_LeaveOneOut=$(mktemp)
tmp_noSmallCohort=$(mktemp)

# Get path to the .in file with the paths to all datasets input to MR-MEGA
get_in_file () {
    result_file=$1
    echo /path/to/$(echo $result_file | sed 's@multi-ancestry_@@g' | sed 's@\.in.*@@g').in; 
} 

# How many datasets were input to MR-MEGA, in our case this is 5
NUM_TOTAL_DATASETS=$(cat $(get_in_file $OUTPUT_ALL_DATASETS) | wc -l)

# Add extra beta/se columns to MR-MEGA output before concatenating runs with different PCs
if [[ $NUM_TOTAL_DATASETS == 5 ]]; then
   awk -F '\t' 'OFS="\t" {$15=$15"\t""\t"; print }' $OUTPUT_ALL_DATASETS | sed "1s/se_2\t\t/se_2\tbeta_3\tse_3/"> $tmp_initial
else
    cat $OUTPUT_ALL_DATASETS > $tmp_initial
fi

echo There are a total of $NUM_TOTAL_DATASETS datasets included.
echo Merging datasets…

ls $LEAVE_ONE_OUT_FILES | while read FILE
do
    NUM_CURR_DATASETS=$(cat $(get_in_file $FILE) | wc -l)
    STUDIES_REMOVED=$(echo $FILE | sed 's@.*\.no@@g' | sed 's@\.in.*@@g' | sed 's@UKB_Proxy_@UKBProxy@g' | sed 's@_@,@g')
    # Add empty columns for the missing betas
    COLS_FORMATTED=$(mktemp)
    # Add the missing beta and se columns as empty tabs
    if [[ $NUM_CURR_DATASETS == 5 ]]; then
       awk -F '\t' 'OFS="\t" {$15=$15"\t""\t"; print }' ${FILE} > $COLS_FORMATTED
    elif [[ $NUM_CURR_DATASETS == 4 ]]; then
        awk -F '\t' 'OFS="\t" {$13=$13"\t""\t""\t""\t"; print }' ${FILE} > $COLS_FORMATTED
    fi
    # Filter out variants with SmallCohortCount error (Ncohort < 4)
    sed "s/$/\t${STUDIES_REMOVED}/" $COLS_FORMATTED |  awk '(NR==1) || $(NF-1) != "SmallCohortCount"' | tail -n+2 >> $tmp_LeaveOneOut
    echo Dataset missing ${STUDIES_REMOVED} has been merged…
done

# Filter the merged results and add "None" for missing study column
sed "s/$/\tNone/" $tmp_initial |sed "1s/None/Removed_studies/" | awk '(NR==1) || $(NF-1) != "SmallCohortCount"' >> $tmp_noSmallCohort

echo Filtering for unique variants with largest Ncohort…
# Sort by Ncohort column so that those with 5 cohorts are before those with 4
# Select unique markernames based on first occurance
cat $tmp_noSmallCohort $tmp_LeaveOneOut | sort -r -k8,8 | sort -u -k1,1 | sort -n -k2,2 -k3,3 > ${OUTPUT_ALL_DATASETS/in/META}
echo Done with ${OUTPUT_ALL_DATASETS/in/META}.
