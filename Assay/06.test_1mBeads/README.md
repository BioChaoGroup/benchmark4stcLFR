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
