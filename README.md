# MA_MA_meta
Repository for ongoing multi-ancestry meta-analysis efforts at the NIH Center for Alzheimer's and Related Dementias (CARD).

## Abstract
Genome-wide association studies (GWAS) of Alzheimer’s disease are predominantly carried out in European ancestry individuals despite the known variation in genetic architecture and disease prevalence across global populations. We leveraged published and _de novo_ GWAS from European, East Asian, African American, and Caribbean Hispanic populations to perform the largest multi-ancestry GWAS meta-analysis of Alzheimer’s disease to date. This method allowed us to identify two independent novel disease-associated loci on chromosome 3. We also leveraged diverse haplotype structures to fine-map nine loci and globally assessed the heterogeneity of known risk factors across populations. Additionally we compared the generalizability of multi-ancestry- and single-ancestry-derived polygenic risk scores in a three-way admixed Colombian population. Our findings highlight the importance of multi-ancestry representation in uncovering and understanding putative factors that contribute to Alzheimer’s disease risk.

![FIGURE1](https://github.com/NIH-CARD/MA_MA_meta/blob/main/Figures/Figure1.png)

## Links to manuscript

1. [Preprint](https://www.medrxiv.org/content/10.1101/2022.08.04.22278442v1)
2. Publication TBD

## Overview

1. [Meta-analysis](https://github.com/NIH-CARD/MA_MA_meta/tree/main/Meta.Analysis)
    - [MR-MEGA](https://github.com/NIH-CARD/MA_MA_meta/tree/main/Meta.Analysis/MR-MEGA): Performing a multi-ancestry meta-analysis using MR-MEGA v0.2.   
    - [Random and Fixed Effects](https://github.com/NIH-CARD/MA_MA_meta/tree/main/Meta.Analysis/Random_and_fixed_effects): Performing multi-ancestry fixed and random effects meta-analyses using PLINK v1.9.   

2. [Post-GWAS Analyses](https://github.com/NIH-CARD/MA_MA_meta/tree/main/Post.GWAS.Analyses)
    - [Identifying Novel Loci](https://github.com/NIH-CARD/MA_MA_meta/tree/main/Post.GWAS.Analyses/Novel_loci): Identifying potentially novel loci by comparing to a list of known AD loci.  
    - [Fine-mapping](https://github.com/NIH-CARD/MA_MA_meta/tree/main/Post.GWAS.Analyses/Fine-mapping): Fine-mapping using MR-MEGA to create 99% credible sets. 
    - [SMR](https://github.com/NIH-CARD/MA_MA_meta/tree/main/Post.GWAS.Analyses/SMR): Summary-based Mendelian Randomization to nominate susceptibility genes at novel loci.
    - [PRS](https://github.com/NIH-CARD/MA_MA_meta/tree/main/Post.GWAS.Analyses/PRS): Performing polygenic risk score analyses and assessing performance with a receiver operator curve. 
