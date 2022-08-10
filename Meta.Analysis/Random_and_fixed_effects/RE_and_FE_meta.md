
# PLINK meta-analysis

##  Reformat summary statistics for PLINK
```
Rscript --vanilla reformat_for_PLINK_META.R Bellenguez2022_for_MRMEGA.no_multiAllelics_indels.MAF_0.01.txt 85934 401577
Rscript --vanilla reformat_for_PLINK_META.R FinngenR6_for_MRMEGA.no_multiAllelics_indels.MAF_0.01.txt 7329 131102
Rscript --vanilla reformat_for_PLINK_META.R Shigemizu2021_for_MRMEGA.no_multiAllelics_indels.MAF_0.01.txt 3962 4074
Rscript --vanilla reformat_for_PLINK_META.R Kunkle2021_for_MRMEGA.no_multiAllelics_indels.MAF_0.01.txt 2748 5222
Rscript --vanilla reformat_for_PLINK_META.R CarHisp_for_MRMEGA.no_multiAllelics_indels.MAF_0.01.txt 1088 1152
```

## Run random and fixed effects meta-analyses in PLINK
- Use logscale to indicate that betas are used 
- Use the qt handle to output betas instead of odds ratios
- Include study handle to output betas from input studies along with betas from the meta-analysis
```
module load plink # 1.9

GWAS_FILES=$(ls ./{Kunkle2021,CarHisp,Shigemizu2021,Bellenguez2022,FinngenR6}_for_PLINK_META.no_multiAllelics_indels.MAF_0.01.txt | paste -sd " ")

plink --meta-analysis $GWAS_FILES \
+ logscale study qt \
--out multi-ancestry_PLINK
```
