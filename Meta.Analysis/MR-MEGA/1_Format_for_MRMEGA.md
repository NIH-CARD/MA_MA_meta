# Format summary statistics for MR-MEGA

## Bellenguez 2022 (European)
Download summary statistics
```
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz
```
Convert from hg38 to hg19
```
Rscript --vanilla liftover_summary_stats.R \
--summary_stats /path/to/GCST90027158_buildGRCh38.tsv.gz \
--from_build "hg38" \
--to_build "hg19" \
--chr "chromosome" \
--pos "base_pair_location" \
--convert_markername "variant_alternate_id" \
--out Bellenguez2022_liftedHg19.txt
```
Format for METAL (fixed effects for PRS)
```
R
require(dplyr)
require(data.table)
data <- fread("Bellenguez2022_liftedHg19.txt",header=T)
data$MARKERNAME <- paste(data$chromosome, data$base_pair_location, sep=":")
nrow(data) # 21058985
data2 <- data %>% select(MARKERNAME, other_allele, effect_allele, effect_allele_frequency, beta, standard_error, p_value)
colnames(data2) <- c("markerID", "alternateAllele", "effectAllele", "effectAlleleFreq", "b", "se", "p")  
fwrite(data2,file="Bellenguez2022_for_METAL.txt",quote=FALSE,row.names=F,sep="\t")
q('no')
```
Format for MR-MEGA
```
Rscript --vanilla reformat_metal_input_for_mrmega.R Bellenguez2022_for_METAL.txt 457511 
```

## Finngen R6 (European Finnish)
Download summary statistics
```
module load google-cloud-sdk 
gsutil cp gs://finngen-public-data-r6/summary_stats/finngen_R6_G6_ALZHEIMER_EXMORE.gz .
gsutil cp gs://finngen-public-data-r6/summary_stats/finngen_R6_G6_ALZHEIMER_EXMORE.gz.tbi .
```
Convert from hg38 to hg19
```
Rscript --vanilla liftover_summary_stats.R \
--summary_stats /path/to/finngen_R6_G6_ALZHEIMER_EXMORE.gz \
--from_build "hg38" \
--to_build "hg19" \
--chr "#chrom" \
--pos "pos" \
--out FinngenR6_liftedHg19.txt
```
Format for METAL (fixed effects for PRS)
```
R
require(dplyr)
require(data.table)
data <- fread("FinngenR6_liftedHg19.txt",header=T)
data$MARKERNAME <- paste(data$"#chrom", data$pos, sep=":")
nrow(data)# 15749231
data2 <- data %>% select(MARKERNAME, ref, alt, af_alt, beta, sebeta, pval)
colnames(data2) <- c("markerID", "alternateAllele", "effectAllele", "effectAlleleFreq", "b", "se", "p")  
fwrite(data2,file="FinngenR6_for_METAL.txt",quote=FALSE,row.names=F,sep="\t")
q('no')
```
Format for MR-MEGA
```
Rscript --vanilla reformat_metal_input_for_mrmega.R FinngenR6_for_METAL.txt 138431
```

## Kunkle 2021 (African American)
Download summary statistics from NIAGADS accession NG00100 and format for MR-MEGA
```
R
require(dplyr)
require(data.table)
data <- fread("path/to/Kunkle2020_ADGC_AA_META_Model1_SummaryStats.withAlleleFreqs.txt",header=T)
data$MARKERNAME <- paste(data$Chr, data$Pos, sep=":")
data$OR <- exp(data$Beta)
data$OR_95U <- exp(data$Beta + 1.96*data$SE)
data$OR_95L <- exp(data$Beta - 1.96*data$SE)
data$N <- 7970 
data2 <- data %>% select(MARKERNAME, Chr, Pos, Effect_allele, Non_Effect_allele, Effect_allele_Freq, OR, OR_95U, OR_95L, N, Pvalue)
colnames(data2) <- c("MARKERNAME", "CHROMOSOME", "POSITION", "EA", "NEA", "EAF", "OR", "OR_95U", "OR_95L", "N", "P")  
fwrite(data2,file="Kunkle2021_for_MRMEGA.txt",quote=FALSE,row.names=F,sep="\t")
q('no')
```

## Shigemizu 2021 (East Asian)
Download summary statistics
```
wget https://humandbs.biosciencedbc.jp/files/hum0237/hum0237_v1_gwas_v1.zip
```
Format for MR-MEGA
```
R
require(dplyr)
require(data.table)
data <- read.table(unzip("/data/CARD/AD/summary_stats/diverse_ancestry/Shigemizu_2021/hum0237_v1_gwas_v1.zip"), header=T)
data$MARKERNAME <- paste(data$CHR, data$BP, sep=":")
N_cases <- 3962
N_controls <- 4074
data$EAF <- (data$MAF_A*N_cases + data$MAF_U*N_controls)/(N_cases + N_controls)
data$N <- N_cases + N_controls
data2 <- data %>% select(MARKERNAME, CHR, BP, A1, A2, EAF, OR, U95, L95, N, P)
colnames(data2) <- c("MARKERNAME", "CHROMOSOME", "POSITION", "EA", "NEA", "EAF", "OR", "OR_95U", "OR_95L", "N", "P")  
fwrite(data2,file="Shigemizu2021_for_MRMEGA.txt",quote=FALSE,row.names=F,sep="\t")
q('no')
```

## Caribbean Hispanic
Convert from hg38 to hg19
```
Rscript liftover_summary_stats.R \
--summary_stats "path/to/CaribbeanHispanic.hybrid" \
--from_build "hg38" \
--to_build "hg19" \
--chr "#CHROM" --pos "POS" \
--convert_markername "ID" \
--out CarHisp_liftedHg19.txt
```
Format for MR-MEGA
```
R
require(dplyr)
require(data.table)
data <- fread("CarHisp_liftedHg19.txt",header=T)
data <- data %>% filter(ERRCODE == ".")
data$MARKERNAME <- paste(data$"#CHROM", data$POS, sep=":")
data$b <- log(data$OR, base = exp(1))
data$OR_95U <- exp(data$b + 1.96*data$"LOG(OR)_SE")
data$OR_95L <- exp(data$b - 1.96*data$"LOG(OR)_SE")
data$NEA<- ifelse(data$A1 == data$REF, data$ALT, data$REF)
data$N <- 2240
data2 <- data %>% select(MARKERNAME, "#CHROM", POS, A1, NEA, A1_FREQ, OR, OR_95U, OR_95L, N, P)
colnames(data2) <- c("MARKERNAME", "CHROMOSOME", "POSITION", "EA", "NEA", "EAF", "OR", "OR_95U", "OR_95L", "N", "P")
fwrite(data2,file="CarHisp_for_MRMEGA.txt",quote=FALSE,row.names=F,sep="\t")
q('no')
```

## Remove multi-allelics/indels and subset MAF > 0.01 prior to MR-MEGA
Remove multi-allelics and indels from all datasets
```
ls *for_MRMEGA.txt | while read line
do
    Rscript remove_multi_and_indels.R $line
done
```
Subset MAF > 0.01
```
ls *for_MRMEGA.no_multiAllelics_indels.txt | while read line
do
    awk '(NR==1) || ($6 > 0.01 && $6 < 0.99)' $line > ${line%%.*}.no_multiAllelics_indels.MAF_0.01.txt
done
```
