#!/bin/bash
# sh FE_METAL.sh $PREFIX $PATTERN
PREFIX=$1
PATTERN=$2

ls $PATTERN > ${PREFIX}_GWAS_FILES.txt
get_colname () {
colnumber=$1
file=$2
head -1 $file |cut -f$colnumber
}

make_metal_file () {
echo """
SCHEME  STDERR
AVERAGEFREQ ON
MINMAXFREQ ON

"""

cat ${PREFIX}_GWAS_FILES.txt | while read line
do
	echo ""
	echo "# === DESCRIBE AND PROCESS INPUT FILE ==="
	echo MARKER $(get_colname 1 $line)
	echo ALLELE $(get_colname 3 $line) $(get_colname 2 $line)
	echo FREQ   $(get_colname 4 $line)
	echo EFFECT $(get_colname 5 $line)
	echo STDERR $(get_colname 6 $line)
	echo PVALUE $(get_colname 7 $line)
	echo PROCESS $line
done
echo """
OUTFILE ${PREFIX}_metal .tbl
ANALYZE HETEROGENEITY
QUIT
"""
}

module load metal
make_metal_file > ${PREFIX}_metal_template.txt
metal ${PREFIX}_metal_template.txt
