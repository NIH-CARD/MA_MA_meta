#!/usr/bin/env Rscript
## start like this
# Rscript find_novel_loci.R \
# --summary_stats path/to/sumstats \
# --chr "Chromosome" \
# --pos "Position" \
# --markername "MarkerName" \
# --P "P.value_association" \
# --out "path/to/summstats/sumstats.loci.csv" \
# --known_loci "known_loci.txt"

library(optparse)
library(data.table)
library(dplyr)
 
option_list = list(
  make_option(c("-s", "--summary_stats"), type="character", default=NULL, metavar="character"),
  make_option(c("-c", "--chr"), type="character", default="Chromosome", metavar="character"),
  make_option(c("-l", "--pos"), type="character", default="Position", metavar="character"),
  make_option(c("-m", "--markername"), type="character", default="MarkerName", metavar="character"),
  make_option(c("-p", "--P"), type="character", default="P.value_association", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", metavar="character"),
  make_option(c("-k", "--known_loci"), type="character", default=NULL, metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

data <- fread(opt$summary_stats,header=T)

## Read in known loci df from Bellenguez et al 2022 Sup Table 5
known_loci <- read.table(opt$known_loci,header=T)
old_names <- c(opt$chr,opt$pos,opt$P,opt$markername)
new_names <- c("Chromosome", "Position", "Pvalue", "MarkerName")
setnames(data,old=old_names, new=new_names)
data[, c("Chromosome", "Position", "Pvalue")] <- lapply(data[, c("Chromosome", "Position", "Pvalue")], as.numeric)
hits <- filter(data, Pvalue < 5e-8)
hits <- hits[order(hits$Pvalue),]

## If the SNP is present in one of the ranges from the known loci list, it is known
## If the SNP is not present in the known loci ranges, create a 500kb range around the SNP
known_locus_func <- function(row) {
   row_Pos <- as.numeric(row["Position"])
   row_Chr <- as.numeric(row["Chromosome"])
   locus <- known_loci %>% filter(chr == row_Chr & start <= row_Pos  & end >= row_Pos)
   out_df <- data.frame(MarkerName=row["MarkerName"], 
	Chromosome=row["Chromosome"], Position=row["Position"], Pvalue=row["Pvalue"])
    if (nrow(locus) == 0) {
        out_df$Locus_start <- row_Pos - 250000
        out_df$Locus_end <- row_Pos + 250000
        out_df$Is_known_locus="NO"
    } else {
        out_df$Locus_start <- locus$start
        out_df$Locus_end <- locus$end
        out_df$Distance_to_known_locus <- 0
        out_df$Is_known_locus="YES" }
    return(out_df)
}

## Assign locus information (start, end) using information from the known loci list
loci <- apply(hits, 1, known_locus_func) %>% bind_rows() %>% data.frame(row.names=NULL)
numeric_cols <- c("Chromosome", "Position", "Pvalue", "Locus_start", "Locus_end")
loci[, numeric_cols] <- lapply(loci[, numeric_cols], as.numeric)

## Group known loci
loci_summary_part1 <- filter(loci, Is_known_locus == "YES") %>% 
	group_by(Chromosome, Locus_start, Locus_end, Is_known_locus, Distance_to_known_locus) %>% 
	summarize(Num_hits=n(), Min_Pvalue=min(Pvalue), Lead_SNP= MarkerName[Pvalue==Min_Pvalue])

## Get distance to nearest known locus
distance_to_known_locus_func <- function(input_Chr, input_Pos) {
    positions_to_compare <- known_loci %>% filter(chr == input_Chr) %>% pull(start, end)
    distance <- min(abs(input_Pos-positions_to_compare))
    return(distance)
}

## Group potentially novel loci into 500kb ranges
part2 <- filter(loci, Is_known_locus == "NO") 
for (i in 1:nrow(part2)) {
   existing_locus <- part2[1:i-1,] %>% filter(Chromosome == part2[i, "Chromosome"] & Locus_start <= part2[i, "Position"] & Locus_end >= part2[i, "Position"]) %>% head(1)
   part2[i, "Lead_SNP"] <- ifelse(nrow(existing_locus) == 0, part2[i, "MarkerName"], existing_locus$MarkerName)
   part2[i, "Distance_to_known_locus"] <- distance_to_known_locus_func(part2[i, "Chromosome"], part2[i, "Position"])
}

loci_summary_part2 <- part2 %>% 
	group_by(Chromosome, Lead_SNP, Is_known_locus) %>% 
	summarize(Num_hits=n(), Min_Pvalue=min(Pvalue), Locus_start=first(Locus_start), Locus_end=first(Locus_end), Distance_to_known_locus=min(Distance_to_known_locus))

## Combine known and potentially novel loci
loci_summary <- bind_rows(loci_summary_part1, loci_summary_part2)
fwrite(loci_summary, file=opt$out, quote=FALSE,row.names=F,col.names=T, sep=",")
