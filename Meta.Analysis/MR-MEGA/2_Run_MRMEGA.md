# Run MR-MEGA

## Download and install MR-MEGA
Download the most recent version from https://genomics.ut.ee/en/tools 
```
wget https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip
unzip MR-MEGA_v0.2.zip 
```
Compile
```
make 
```
Download additional scripts provided by MR-MEGA authors
```
wget https://tools.gi.ut.ee/tools/{fixP.r,manh.r,qq.r}
```
Make fixP.r more efficient by changing read.table/write.table to fread/fwrite
```
sed -i 's/read.table/fread/g' fixP.r sed -i 's/write.table/fwrite/g' fixP.r sed -i 's/data<-/require(data.table)\ndata<-/g' fixP.r
```

## Create input files for MR-MEGA
Input file with paths to all datasets
```
echo "/path/to/Bellenguez2022_for_MRMEGA.no_multiAllelics_indels.MAF_0.01.txt
/path/to/FinngenR6_for_MRMEGA.no_multiAllelics_indels.MAF_0.01.txt
/path/to/Kunkle2021_for_MRMEGA.no_multiAllelics_indels.MAF_0.01.txt
/path/to/Shigemizu2021_for_MRMEGA.no_multiAllelics_indels.MAF_0.01.txt
/path/to/CarHisp_for_MRMEGA.no_multiAllelics_indels.MAF_0.01.txt" > MR-MEGA.in
```
Make input files for "leave-one-out" method with each combination of 4 datasets
```
grep -v Bellenguez MR-MEGA.in > MR-MEGA.noBellenguez.in
grep -v Finngen MR-MEGA.in > MR-MEGA.noFinngen.in
grep -v Kunkle MR-MEGA.in > MR-MEGA.noKunkle.in
grep -v Shigemizu MR-MEGA.in > MR-MEGA.noShigemizu.in
grep -v CarHisp MR-MEGA.in > MR-MEGA.noCarHisp.in
```

## Run MR-MEGA
In reality will want to run in parallel
```
ls *in | while read line
do
    sh run_MRMEGA_max_PCs.sh $line
done
```

## Clean MR-MEGA output
```
ls *result | while read line
do
    sh fixP_MRMEGA.sh $line
done
```

## Concatenate MR-MEGA outputs 
Combine the MR-MEGA output including all studies with the leave-one-out files with 4 input datasets
```
sh merge_leave_one_out.sh "multi-ancestry_MR-MEGA.in.MAX_PCs.result.fixedP" "multi-ancestry_MR-MEGA.no*.in.MAX_PCs.result.fixedP"
```
