#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Rscript --vanilla munge_for_PRS.R $FILENAME
FILENAME = args[1]
require(MungeSumstats)
require(data.table)
require(dplyr)
format_sumstats(
FILENAME,
ref_genome = "GRCh37",
snp_ids_are_rs_ids = FALSE,
bi_allelic_filter= TRUE,
force_new=TRUE,
sort_coordinates = TRUE,
nThread = 4,
save_path=gsub(pattern="forMungeSumstats", replacement="forPRS", FILENAME),
log_folder_ind = TRUE,
log_mungesumstats_msgs = TRUE,
log_folder = paste("./log_files/", gsub(pattern=".forMungeSumstats.txt", replacement="", FILENAME), sep="")
)
