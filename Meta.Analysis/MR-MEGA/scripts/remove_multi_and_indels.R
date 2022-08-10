#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Rscript remove_multi_and_indels.R $FILENAME
FILENAME = args[1]
require(dplyr)
require(data.table)
data <- fread(FILENAME,header=T)
multi_allelics <- duplicated(data$MARKERNAME)
data2 <- data[!multi_allelics] 
data3 <- data2 %>% filter(nchar(NEA) == 1) %>% filter(nchar(EA) == 1)
fwrite(data3,file=gsub("txt", "no_multiAllelics_indels.txt", FILENAME),quote=FALSE,row.names=F,sep="\t")
