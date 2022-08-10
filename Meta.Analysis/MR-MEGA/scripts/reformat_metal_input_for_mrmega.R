#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Rscript --vanilla reformat_metal_input_for_mrmega.R $FILENAME $SAMPLE_SIZE
FILENAME = args[1]
SAMPLE_SIZE = args[2]
require(dplyr)
require(data.table)
data <- fread(FILENAME,header=T)
colnames(data) <- c("MARKERNAME", "NEA", "EA", "EAF", "BETA", "SE", "P")
data$OR <- exp(data$BETA)
data$OR_95U <- exp(data$BETA + 1.96*data$SE)
data$OR_95L <- exp(data$BETA - 1.96*data$SE)
data$CHROMOSOME <- gsub(data$MARKERNAME, pattern=":.*",replacement="")
data$POSITION <- gsub(data$MARKERNAME, pattern=".*:",replacement="")
data$EA <- toupper(data$EA)
data$NEA <- toupper(data$NEA)
data$N <- SAMPLE_SIZE
data2 <- data %>% select(MARKERNAME, CHROMOSOME, POSITION, EA, NEA, EAF, OR, OR_95U, OR_95L, N, P)
fwrite(data2,file=gsub("for_METAL.txt", "for_MRMEGA.txt", FILENAME),quote=FALSE,row.names=F,sep="\t")
