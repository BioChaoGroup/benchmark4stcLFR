### check barcodes position
```sh
for i in `ls ../../Source/*/*/*2.fq.gz`; do echo $i; perl $mSeq/src/barcodePos.pl ../../Source/JUNE/barcode.list $i;done > checkBarcode.log
```

### pick 5M reads for test:
```sh
zcat Aprl_84L02_Ori.1.fq.gz|head -20000000|pigz > Aprl_84L02_5M.1.fq.gz
zcat Aprl_84L02_Ori.2.fq.gz|head -20000000|pigz > Aprl_84L02_5M.2.fq.gz
```

### Load test env:
```
export mSeqOri="/hwfssz1/ST_META/EE/bin/metaSeq"
export mSeqDev="/hwfssz1/ST_META/P18Z10200N0059_stLFR_SZW/USER/fangchao/metaSeq"
export mfastp="/hwfssz1/ST_META/P18Z10200N0059_stLFR_SZW/USER/fangchao/fastp"
```

### detect barcodes distribution
```
perl $mSeqDev/src/OAs1 -i Aprl_84L02_5M.1.fq.gz,Aprl_84L02_5M.2.fq.gz \
-pfx Aprl_84L02_5M -stLFR 2 -b ../../Source/JUNE/barcode.list -debug 1
```
> Running time: 267 sec.
