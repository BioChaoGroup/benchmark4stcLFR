### check barcodes position
```sh
for i in `ls ../../Source/*/*/*2.fq.gz`; do echo $i; perl $mSeq/src/barcodePos.pl ../../Source/JUNE/barcode.list $i;done > checkBarcode.log
```
### pick 1M reads for test:
```sh
mkdir -p tmp
zcat Aprl_84L02_Ori.1.fq.gz|head -4000000|pigz > tmp/Aprl_84L02_1M.1.fq.gz
zcat Aprl_84L02_Ori.2.fq.gz|head -4000000|pigz > tmp/Aprl_84L02_1M.2.fq.gz

pigz -dc Aprl_84L02_Ori.1.fq.gz|tail -4000000|pigz > tmp/Aprl_84L02t1M.1.fq.gz &
pigz -dc Aprl_84L02_Ori.2.fq.gz|tail -4000000|pigz > tmp/Aprl_84L02t1M.2.fq.gz &
```

### pick 5M reads for test:
```sh
zcat Aprl_84L02_Ori.1.fq.gz|head -20000000|pigz > tmp/Aprl_84L02_5M.1.fq.gz &
zcat Aprl_84L02_Ori.2.fq.gz|head -20000000|pigz > tmp/Aprl_84L02_5M.2.fq.gz &
```

### Load test env:
```
export mSeqOri="/hwfssz1/ST_META/EE/bin/metaSeq"
export mSeqDev="/hwfssz1/ST_META/P18Z10200N0059_stLFR_SZW/USER/fangchao/metaSeq"
export mfastp="/ldfssz1/ST_META/share/User/fangchao/fastp"
export mcOMG="/ldfssz1/ST_META/share/User/fangchao/Omics_pipeline/MetaGenomics"
```

### detect barcodes distribution
```
iPfx="Aprl_84L02_1M"
perl $mcOMG/bin/OAs1 -i tmp/$iPfx.1.fq.gz,tmp/$iPfx.2.fq.gz \
-pfx tmp/$iPfx -stLFR 2 -b ../../Source/JUNE/barcode.list -debug 1 > tmp/$iPfx.log
```
Reports generated by [checkSplit.visulization.Rmd](./checkSplit.visulization.Rmd).

> 1M: Running time: 40 sec.
> 5M: Running time: 193 sec.

### Efficient test

```
iPfx="Aprl_84L02_1M"

#Python
python $mSeqDev/stlfr_split.py -r1 tmp/$iPfx.1.fq.gz -r2 tmp/$iPfx.2.fq.gz -b ../../Source/JUNE/barcode.list -bl 54 -o tmp/$iPfx.py -fastq &> tmp/$iPfx.py.log
#perl
perl $mcOMG/bin/OAs1 -i tmp/$iPfx.1.fq.gz,tmp/$iPfx.2.fq.gz \
-pfx tmp/$iPfx.pl -stLFR 2 -b ../../Source/JUNE/barcode.list -debug 1 > tmp/$iPfx.pl.log
#C++ # 1 threads
$mfastp/fastp --stLFR_barcode_file ../../Source/JUNE/barcode.list \
--in1 tmp/$iPfx.1.fq.gz --in2 tmp/$iPfx.2.fq.gz --disable_adapter_trimming \
--out1 tmp/$iPfx.cp.1.fq.gz --out2 tmp/$iPfx.cp.2.fq.gz --disable_quality_filtering -w 1\
&> tmp/$iPfx.cp1.log
#C++ # 6 threads
$mfastp/fastp --stLFR_barcode_file ../../Source/JUNE/barcode.list \
--in1 tmp/$iPfx.1.fq.gz --in2 tmp/$iPfx.2.fq.gz --disable_adapter_trimming \
--out1 tmp/$iPfx.cp.1.fq.gz --out2 tmp/$iPfx.cp.2.fq.gz --disable_quality_filtering -w 6
&> tmp/$iPfx.cp6.log
```
Results: 1M reads cost:
> python: 52 seconds
> perl:   40 seconds
> c++:    28 seconds (-w 1, actually it comsumed 200% cpu)
> c++:    12 seconds (-w 6, actually it comsumed no more than 600% cpu)

So we can estimated the time for 600 million reads sequences:
> python: 8.6 hr
> perl:   6.7 hr
> c++:    4.7 hr
> c++:    2.0 hr

# RCA test

sample info:
> RCA14  3h  100pg

```bash
SRC14=/path/to/RCA14/rawdata
mkdir -p ../../Results/RCA14/clean
ln -s ../../Results/RCA14

#prep fq
cd RCA14/clean
mkdir tmp
pigz -p 8 -dc $SRC14/split_read.1.fq.gz | perl -e 'while(<>){
  if($.%4==1){@a=split /#|\//;$a[1]=~/(\d+)_(\d+)_(\d+)/;
  $b1 = sprintf("%04d",$1);$b2 = sprintf("%04d",$2);$b3 = sprintf("%04d",$3);
  print "$a[0]/$b1\_$b2\_$b3/1\n"}else{print $_}}' | paste - - - - | \
  sort -T ./tmp -k2,2 -t "/" | tr "\t" "\n" > fastp.sort.1.fq &

pigz -p 8 -dc $SRC14/split_read.2.fq.gz | perl -e 'while(<>){
  if($.%4==1){@a=split /#|\//;$a[1]=~/(\d+)_(\d+)_(\d+)/;
  $b1 = sprintf("%04d",$1);$b2 = sprintf("%04d",$2);$b3 = sprintf("%04d",$3);
  print "$a[0]/$b1\_$b2\_$b3/2\n"}else{print $_}}' | paste - - - - | \
  sort -T ./tmp -k2,2 -t "/" | tr "\t" "\n" > fastp.sort.2.fq &
wait
metabbq randomlyPick.pl -i fastp.sort.1.fq -o fastp.sort.1 -n 10 -m 3 -s 14 -t 16 &
metabbq randomlyPick.pl -i fastp.sort.2.fq -o fastp.sort.2 -n 10 -m 3 -s 14 -t 16 -v
cd ../../
mkdir -p RCA_00/clean
ln -s ../../RCA14/clean/fastp.sort.1_00.fq RCA_00/clean/fastp.sort.1.fq
ln -s ../../RCA14/clean/fastp.sort.2_00.fq RCA_00/clean/fastp.sort.2.fq
```
