# Polygenic Risk Score Analyses

## Reference alignment and QC of base summary statistics

- Align summary statistics to reference genome
- Perform basic quality control (QC) in preparation for PRSice
- Remove APOE region from summary statistics

## Munge summary statistics
Adjust the file headers for MungeSumstats and filter for variants present in at least 3 datasets

Fixed effects
```
awk 'OFS="\t" {if (gsub(/NA/, "") < 3) print $1, $2, $3, $4, $5, $7, $9}' /path/to/multi-ancestry_PLINK.meta > multi-ancestry_PLINK.meta.FE.forMungeSumstats.txt
```
Random effects
```
awk 'OFS="\t" {if (gsub(/NA/, "") < 3) print $1, $2, $3, $4, $5, $8, $10}' /path/to/multi-ancestry_PLINK.meta | sed '1s/P(R)/P/g' | sed '1s/BETA(R)/BETA/g' > /path/to/multi-ancestry_PLINK.meta.RE.forMungeSumstats.txt
```
Reformat individual ancestry files
```
head -1 /path/to/multi-ancestry_PLINK.meta.FE.forMungeSumstats.txt | sed 's/BETA/OR/g' > header.txt
for PREFIX in {"Bellenguez_FinngenR6","Kunkle2021","Shigemizu2021","CarHisp"}
do
   sh reformat_mrmega_prs.sh /path/to/${PREFIX}_for_MRMEGA.no_multiAllelics_indels.MAF_0.01.txt $PREFIX
done
```
Align alleles to the reference
```
ls *forMungeSumstats.txt | while read line
do
    Rscript --vanilla munge_for_PRS.R $line
done
```
Update some positions that were converted to scientific notation
```
ls *forPRS.txt | while read line
do
    awk 'OFS="\t" { print $1,$2,sprintf("%.0f", $3),$4,$5,$6,$7,$8; }' $line | sed '1s/0/BP/g' > ${line/forPRS/forPRSice}
done
```
Remove APOE region prior to PRSice clumping

Start and stop of APOE region (hg19) in random effects analysis according to FUMA: chr19:44771444-45806162
```
R
require(data.table)
require(dplyr)
FILENAMES <- list.files(pattern = 'forPRS.txt$')
for (FILENAME in FILENAMES) {
    file=read.table(FILENAME,header = TRUE)
    file=file[!(file$CHR==19 & file$BP > 44771444 & file$BP < 45806162),]
    fwrite(file, paste(gsub('.{3}$', '', FILENAME),"noapoe.txt",sep=""),quote=FALSE,row.names=F)
}
q('no')
```

## Reference alignment and QC of target dataset (TANGL)
Download hg19 reference fasta file
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
```
Align target dataset to reference genome using PLINK2
```
plink2 --bfile TANGL_PRS \
--ref-from-fa force /path/to/hg19.fa \
--make-bed \
--out TANGL_refAligned
```
Perform QC and make file ready for PRSice

Filter on minor allele frequency, individual missingness, snp missingness, and Hardy Weinburg Equilibrium
```
plink2 --bfile TANGL_refAligned \
--maf 0.01 --mind 0.1 --geno 0.1 --hwe 1e-6 \
--keep /path/to/sampleskeep.txt \
--pheno-name PHENO \
--pheno /path/to/covariate_file.txt \
--set-all-var-ids @:#:'$r':'$a' \
--max-alleles 2 \
--make-bed \
--out TANGL_forPRSice
```

## Calculate PCs and format covariate file
```
module load plink # 1.9
FILENAME=TANGL_forPRSice
```
Filter data
```
plink --bfile $FILENAME --maf 0.01 --geno 0.01 --hwe 5e-6 --autosome --exclude exclusion_regions_hg19.txt --make-bed --out $FILENAME.2  
```
Prune SNPs 
```
plink --bfile $FILENAME.2 --indep-pairwise 1000 10 0.02 --autosome --out pruned_data
```
Extract pruned SNPs for PC calculation
```
plink --bfile $FILENAME.2 --extract pruned_data.prune.in --make-bed --out $FILENAME.3 
```
Calculate PCs
```
flashpca --bfile $FILENAME.3 --suffix _filter_pruned_forPCA.txt --numthreads 19
```
Append PCs to covariate file
```
R
library(data.table)
library(dplyr)
pheno=fread('/path/to/covariate_file.txt',header=T)
eigen=fread('eigenvectors_filter_pruned_forPCA.txt',header=T)
covar=left_join(eigen,pheno, by =c('IID','FID'))
colnames(covar)[20]="AGE_ANALYSIS"
covar=covar[,c(1,2,3,4,5,6,7,8,9,10,15,20)]
fwrite(covar,'covariate_file_with_PCs.txt',quote=F,row.names=F,sep="\t")
q('no')
```
## Run PRSICE without APOE
```
sh run_PRSice.sh multi-ancestry_PLINK.meta.FE.forPRSice.noapoe.txt BETA
sh run_PRSice.sh multi-ancestry_PLINK.meta.RE.forPRSice.noapoe.txt BETA
sh run_PRSice.sh Bellenguez_FinngenR6.forPRSice.noapoe.txt OR
sh run_PRSice.sh CarHisp.forPRSice.noapoe.txt OR
sh run_PRSice.sh Kunkle2021.forPRSice.noapoe.txt OR
sh run_PRSice.sh Shigemizu2021.forPRSice.noapoe.txt OR
```
## Run PRSice with APOE
First make SNP lists including the two APOE SNPs rs429358 and rs7412

The following loop will:
   1) Pull  all the summary statistics into R
   2) Filter these using the pruned SNP lists (end with ".snp") from the analysis excluding the APOE region 
   3) Add in the two APOE SNPs (rs429358 and rs7412) that define the ApoE-ε2, ApoE-ε3, and ApoE-ε4 alleles
```
FILENAMES <- list.files(pattern = 'forPRS.txt$')
for (FILENAME in FILENAMES) {
    file=read.table(FILENAME,header = TRUE)
    snplist=read.table(paste('/path/to/','output_',gsub('.{3}$', '', FILENAME),'noapoe.txt.snp',sep=""),header=TRUE)
    file$SNP=paste(file$CHR,":",file$BP,":",file$A1,":",file$A2,sep="")
    apoe <- file %>% filter(CHR==19 & BP==45411941 | CHR==19 & BP==45412079)
    file <- file %>% filter(SNP %in% snplist$SNP)
    file <- rbind(file,apoe)
    fwrite(file,paste(gsub('.{3}$', '', i),"withapoe.txt",sep=""),quote=FALSE,row.names=F)
}
```
Manually add APOE effects from Schwartzentruber et al. since Bellenguez et al. summary statistics don't have the two APOE SNPs
```
apoe <-... # dataframe with effects from Schwartzentruber et al.
file=read.table(FILENAMES[4],header = TRUE)
snplist=read.table('Bellenguez_FinngenR6.forPRSice.noapoe.txt',header=TRUE)
file$SNP=paste(file$CHR,":",file$BP,":",file$A1,":",file$A2,sep="")
file <- file %>% filter(SNP %in% snplist$SNP)
file <- rbind(file,apoe)
fwrite(file,'Bellenguez_FinngenR6.forPRSice.withapoe.txt',quote=FALSE,row.names=F)
```
Next run PRSice with APOE
```
sh apoe_run_PRSice.sh multi-ancestry_PLINK.meta.FE.forPRSice.withapoe.txt BETA
sh apoe_run_PRSice.sh multi-ancestry_PLINK.meta.RE.forPRSice.withapoe.txt BETA
sh apoe_run_PRSice.sh Bellenguez_FinngenR6.forPRSice.withapoe.txt OR
sh apoe_run_PRSice.sh CarHisp.forPRSice.withapoe.txt OR
sh apoe_run_PRSice.sh Kunkle2021.forPRSice.withapoe.txt OR
sh apoe_run_PRSice.sh Shigemizu2021.forPRSice.withapoe.txt OR
```
