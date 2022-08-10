# AUC Generation

- Example of how to generate AUC plot without APOE
- Same process including APOE, but using the "withAPOE" PRSice outputs.

Load libraries 
```
R
require(data.table)
require(dplyr)
require(metafor)
require(tidyverse)
require(ggplot2)
require(pROC)
require(ggpubr)
```
Pull in covariates and PRS results
```
FE=fread('output_multi-ancestry_PLINK.meta.FE.forPRSice.noapoe.txt.best',header=T)
RE=fread('output_multi-ancestry_PLINK.meta.RE.forPRSice.noapoe.txt.best',header=T)
EUR=fread('output_Bellenguez_FinngenR6.forPRSice.noapoe.txt.best',header=T)
EAS=fread('output_Shigemizu2021.forPRSice.noapoe.txt.best',header=T)
CH=fread('output_CarHisp.forPRSice.noapoe.txt.best',header=T)
AFR=fread('output_Kunkle2021.forPRSice.noapoe.txt.best',header=T)
```
Adjust column names
```
colnames(FE)[4]="FE_PRS"
colnames(RE)[4]="RE_PRS"
colnames(EUR)[4]="EUR_PRS"
colnames(EAS)[4]="EAS_PRS"
colnames(CH)[4]="CH_PRS"
colnames(AFR)[4]="AFR_PRS"
```
Merge all dataframes
```
df_list <- list(CH, AFR, EAS, EUR, FE, RE)
merged <- df_list %>% reduce(full_join, by=c("FID","IID","In_Regression"))

DATA=merged
COVARIATES= 'TANGL_PRS_metadata.txt'
TITLE="Excluding APOE"

covar <- fread(COVARIATES, header=T)
full_table <- left_join(DATA,covar, by = c("FID","IID"))
full_table <- cbind(full_table, pheno = full_table$PHENO-1) # convert phenotype to 0/1 rather than 1/2 as in PLINK
full_table$SEX=full_table$SEX-1 # 0 is male, 1 is female
full_table$PHENO=full_table$PHENO-1 #0 is control, 1 is case
```
Make composite PRS, using EAS as a proxy for Native American
```
full_table$composite_PRS=full_table$AFR_PRS*full_table$`AFRICAN %`+full_table$EUR_PRS*full_table$`EUROPEAN %`+full_table$EAS_PRS*full_table$`NATIVE AMERICAN %`
```

## Plot AUC

Store roc object
```
roc.list <- roc(PHENO ~ RE_PRS + FE_PRS + EUR_PRS + EAS_PRS + AFR_PRS + composite_PRS + CH_PRS, data = full_table, smooth=FALSE)
```
Extract auc
```
data.auc <- roc.list %>% map(~tibble(AUC = .x$auc)) %>% bind_rows(.id = "name")
```
Generate labels
```
data.labels <- data.auc %>% 
  mutate(label_short=c("MULTI (RE)", "MULTI (FE)", "EUR", "EAS", "AFR", "COMPOSITE", "CAR HISP"),
  label_AUC=paste0("AUC = ",paste(round(AUC,2))), label_long=paste0(label_short,": AUC = ",round(AUC,2)))
```
Plot
```
png('AUC_without_APOE.png', width = 5, height = 3, units = "in",res = 500)
ggroc(roc.list, size=0.75) +
scale_color_discrete(labels=data.labels$label_long, direction=-1) +
geom_abline(slope = 1 ,intercept = 1, color="black", size=0.4) + # add identity line
theme(panel.background = element_blank(), 
axis.title.x = element_text(size=10, face="plain"),
axis.title.y = element_text(size=10, face="plain"),
panel.border = element_rect(size=0.75, fill = NA), 
panel.grid.major = element_line(colour="lightgrey", size=0.2),
axis.text.x = element_text(size = 8),
axis.text.y = element_text(size = 8),
plot.title = element_text(size = 11),
axis.ticks.length = unit(.15, "cm"),
axis.ticks=element_line(size=0.25),
legend.title=element_blank(), 
legend.key=element_blank()) +
xlab('Specificity') +
ylab('Sensitivity') +
ggtitle(expression(bold(paste("Excluding ", bolditalic("APOE")))))
dev.off()
```
