# Identify novel loci
Search through MR-MEGA and random/fixed effects meta-analyses for potentialy novel loci

MR-MEGA
```
Rscript find_novel_loci.R -s multi-ancestry_MR-MEGA.META.MAX_PCs.result.fixedP --p "P.value_association" -o MRMEGA.result.fixedP.novel.csv -k known_loci.txt 
```
Fixed effects
```
Rscript find_novel_loci.R -s multi-ancestry_PLINK.meta -c 'CHR' -l 'BP' -m 'SNP' -p 'P' -o multi-ancestry_PLINK.meta.FE.csv -k known_loci.txt
```
Random effects
```
Rscript find_novel_loci.R -s multi-ancestry_PLINK.meta -c 'CHR' -l 'BP' -m 'SNP' -p '"P\(R\)"' -o multi-ancestry_PLINK.meta.RE.csv -k known_loci.txt
```
