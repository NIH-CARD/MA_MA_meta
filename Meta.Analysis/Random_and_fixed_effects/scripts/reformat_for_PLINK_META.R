#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Rscript reformat_for_PLINK_META.R $FILENAME $N_CASES $N_CONTROLS
FILENAME = args[1]
N_CASES = args[2]
N_CONTROLS = args[3]
require(dplyr)
require(data.table)
data <- fread(FILENAME,header=T)
data$BETA <- log(as.numeric(data$OR))
data$SE <- (log(as.numeric(data$OR_95U)) - data$BETA)/1.96
data2 <- data %>% select(MARKERNAME, CHROMOSOME, POSITION, EA, NEA, BETA, SE, P)
colnames(data2) <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P")
fwrite(data2,file=gsub("MRMEGA", "PLINK_META", FILENAME),quote=FALSE,row.names=F,sep="\t")
