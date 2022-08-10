#!/bin/bash
# sh fixP_MRMEGA.sh $FILENAME
# From MR-MEGA documentaion: the current C++ library enables to calculate p-values>1e-14, you can use fixP.r script to recalculate p-values in R  based on chisq and ndf values down to p-values>1e-325. 
INPUT=$1
R --slave --vanilla --args input=${INPUT} < path/to/fixP.r
