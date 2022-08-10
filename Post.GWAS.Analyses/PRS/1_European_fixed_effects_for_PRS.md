# Run fixed-effects meta-analysis using METAL

Use fixed effects meta-analysis of Bellenguez et al. 2022 and FinnGen R6 for PRS

Run METAL script
```
sh METAL.sh "Bellenguez_FinngenR6" "FinngenR6_for_METAL.txt Bellenguez2022_for_METAL.txt"
```
Filter meta-analysis results for HetISq < 80
```
R
require(dplyr)
require(data.table)
data <- fread("Bellenguez_FinngenR6_metal1.tbl",header=T)
data2 <- data %>% filter(HetISq < 80)
fwrite(data2,file="Bellenguez_FinngenR6_metal_HetISq80.txt",quote=FALSE,row.names=F,sep="\t")
q('no')
```
Make same format as pre-MR-MEGA files (although these will all be formatted for PRSice)
```
R
FILENAME = "Bellenguez_FinngenR6_metal_HetISq80.txt"
SAMPLE_SIZE = 625942
require(dplyr)
require(data.table)
data <- fread(FILENAME,header=T)
data$OR <- exp(data$Effect)
data$OR_95U <- exp(data$Effect + 1.96*data$StdErr)
data$OR_95L <- exp(data$Effect - 1.96*data$StdErr)
data$Chr <- gsub(data$MarkerName, pattern=":.*",replacement="")
data$Pos <- gsub(data$MarkerName, pattern=".*:",replacement="")
data$EA <- toupper(data$Allele1)
data$NEA <- toupper(data$Allele2)
data$N <- SAMPLE_SIZE
data2 <- data %>% select(MarkerName, Chr, Pos, EA, NEA, Freq1, OR, OR_95U, OR_95L, N, "P-value")
colnames(data2) <- c("MARKERNAME", "CHROMOSOME", "POSITION", "EA", "NEA", "EAF", "OR", "OR_95U", "OR_95L", "N", "P")  
fwrite(data2,file=gsub("metal_HetISq80.txt", "for_MRMEGA.txt", FILENAME),quote=FALSE,row.names=F,sep="\t")
q('no')
```
Remove multi-allelics and indels
```
Rscript remove_multi_and_indels.R Bellenguez_FinngenR6_for_MRMEGA.txt
```
