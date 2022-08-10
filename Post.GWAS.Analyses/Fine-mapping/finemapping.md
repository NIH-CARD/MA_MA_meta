# Fine-mapping

Load packages
```
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import chi2
```
Define fine-mapping function to create 99% credible sets using Bayes' factors from MR-MEGA
```
def finemapLocus(leadSNP, MRMEGA):
    MRMEGA.loc[:,'BF'] = np.exp(MRMEGA.loc[:,'lnBF'])
    metaChrList = [g for _,g in MRMEGA.groupby('Chromosome')]
    credibleSet = list()
    Finemapped = pd.DataFrame()
    for i, row in leadSNP.iterrows():
        row = row.copy()
        chrom = row.chr
        start = row.start
        stop = row.end
        listNum = chrom-1
        chrMeta = metaChrList[listNum]
        locus_selection = (chrMeta.Position >= start) & (chrMeta.Position <= stop)
        locus = chrMeta.loc[locus_selection]
        locus = locus.sort_values(by=['BF'], ascending=False).reset_index()
        denom = locus['BF'].sum()
        #print(locus)
        numer = locus['BF'].iloc[0]
        PP = numer/denom
        count = 1
        while PP < 0.99:
            numer = numer + locus['BF'][count]
            PP = numer/denom
            count = count + 1
        credibleSetDF = locus[0:count]
        credibleSetDF['Locus'] = i+1
        Finemapped_temp = pd.DataFrame({'Locus':i+1, 'NumSNPs':locus.shape[0], 'NumSNPs in 99% Credible Set':count}, index=[i])
        if i==0:
            Finemapped = Finemapped_temp
            credibleSet = credibleSetDF
        else:
            Finemapped = pd.concat([Finemapped, Finemapped_temp])
            credibleSet = pd.concat([credibleSet, credibleSetDF])
    Final = [Finemapped, credibleSet]
    return(Final)

results_dir = 'path/to/summarystats/'
ref_dir = 'path/to/FUMAresults/'
```
Define your risk loci file: we used GenomicRiskLoci.txt from FUMA with R2 > 0.3 and all populations in 1000 Genomes as a reference
```
loci = pd.read_csv(f'{ref_dir}/GenomicRiskLoci.txt', delim_whitespace=True)
sumstat = pd.read_csv('multi-ancestry_MR-MEGA.META.MAX_PCs.result.fixedP', skip_blank_lines=False,sep='\t')
```
Run the finemap function and print both the report and the 99% credible sets
```
finemap_res = finemapLocus(loci, sumstat)
finemap_res[0].to_csv(f"{results_dir}/FineMap.report.txt", sep = "\t", index=False)
finemap_res[1].to_csv(f"{results_dir}/FineMap.99CredibleSets.txt", sep = "\t", index=False)
PP=finemap_res[1].groupby("Locus").apply(lambda x: x['BF']/x['BF'].sum())
finemap_res[1].groupby("Locus").apply(lambda x: x['BF']/x['BF'].sum()).to_csv(f"{results_dir}/PP.txt", sep = "\t", index=False)
```
If you want to make a file with posterior probabilities or counts per set you can do the following:
```finemap_res[1].groupby("Locus").apply(lambda x: x['BF']/x['BF'].sum()).to_csv(f"{results_dir}/PP.txt", sep = "\t", index=False)
finemap_res[1].groupby("Locus").count()["index"].to_csv(f"{results_dir}/count.txt", sep = "\t", index=True)
```
