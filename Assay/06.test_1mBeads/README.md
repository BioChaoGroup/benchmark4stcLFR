# CODE NOTE

### prepare
```
ln -s $mSeq2/Snakefile
ln -s $mSeq2/config.yaml
ln -s $mSeq2/instance/REF

ln -s ../../Results/APR842_00_00
ln -s ../../Results/APR841_00_00
ln -s ../../Results/DEC560_d10_00

```

### snakemake
```
snakemake -np APR842_00_00/BC1M/batch.assemble.BC.sh
snakemake -np APR841_00_00/BC1M/batch.assemble.BC.sh
snakemake -np DEC560_d10_00/BC1M/batch.assemble.BC.sh
#Run the executing shells manully
```


### test
```
#cd into BC dir
metabbq reAlign.sh 16 sort.1.fq sort.2.fq scaffolds.fasta /ldfssz1/ST_META/F16HQSB1SY2636_TumorMeta_ZH/fangchao/metaSeq/instance/REF/ncbi.5Ref.fasta reAlign
metabbq anno.pl /ldfssz1/ST_META/F16HQSB1SY2636_TumorMeta_ZH/fangchao/metaSeq/instance/REF/ncbi.5Ref.faID reAlign/scaf2ref.blast6 > reAlign/scaf2ref.blast6.anno
metabbq binWrite best -i reAlign/scaf2ref.blast6.anno -u -o reAlign/scaf2ref.blast6.anno
column -t reAlign/scaf2ref.blast6.anno.best |les

```
