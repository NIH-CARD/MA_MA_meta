# Run summary-based Mendelian Randomization

- In this repository, we show an example SMR analysis looking at only Meta-Brain eQTLs and Bellenguez et al. 2022 GWAS summary statistics. 
- We also performed SMR using eQTL summary statistics from GTEx, eMeta, and eQTLgen
- The code can be modified to accomodate any of these datasets

Import libraries
```
import zipfile
import pandas as pd
import os
import glob
import numpy as np
import re
```
Save Bellenguez statistics in a compatible format for Meta Brain 
```
gwas_path = 'path/to/gwas/summarystatisticsinbuildhg19'
df = pd.read_csv(gwas_path, sep='\t')
df['n'] = df['n_cases'] + df['n_controls']
df['marker_2'] = df['chromosome'].astype('string') + ':' + df['base_pair_location'].astype('string') + ':' + df['variant_id'] + ':' + df['other_allele'].astype('string') + '_' + df['effect_allele'].astype('string') 
```
Rename columns for SMR
```
df2 = df[['marker_2', 'effect_allele', 'other_allele', 'effect_allele_frequency', 'beta', 'standard_error', 'p_value', 'n']]
df2.columns = ['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'n']
df2 = df2.dropna()
df2.to_csv('/path/to/AD_metabrain_version_smr.txt', sep='\t', header=True, index=False)
```
Convert eQTL stats to hg19
```
Rscript --vanilla liftover_summary_stats.R \
--summary_stats /path/to/eqtl_buildGRCh38.tsv \
--from_build "hg38" \
--to_build "hg19" \
--chr "chromosome" \
--pos "base_pair_location" \
--convert_markername "variant_id" \
--out /path/to/eqtl_buildhg19.txt
```
Run SMR with Meta Brain
- This needs to be run for each chromosome separately
- This example is for the basal ganglia but we ran on all brain regions. 
```
smr --gwas-summary /path/to/AD_metabrain_version_smr.txt --beqtl-summary /path/to/2020-05-26-Basalganglia-EUR-[[INSERT CHROMOSOME NUMBER]]-SMR-besd --out AD_expression_Basalganglia_metaBrain_SMR_chr[[INSERT CHROMOSOME NUMBER]] --bfile /path/to/1kg_eur_1pct_ref_panel_gr38_metaBrain_variant_names --smr-multi
```
