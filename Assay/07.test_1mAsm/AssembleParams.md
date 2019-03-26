
```bash
SAM0=APR842_00
rTag=SUB_2R100_0D10
mkdir -p $SAM0/AS1M
rsync -atvP --exclude=scaffolds* --exclude=K* ../TT1M/$rTag.Assemble_mashBC $SAM0/AS1M
snakemake -s test.smk
```

# Global Parameters
```bash
SAM0=APR842_00
```

# Run 1
```bash
tag=k77
#use parameters: -k 21,33,55,77
snakemake -s test.smk -j -np $SAM0/AS1M/batch.assemble.BC.sh
cp -r APR842_00/AS1M/{Assemble_mashBC,SUB_2R100_0D10.Assemble_mashBC} # backup
sh $SAM0/AS1M/batch.assemble.BC.sh
```

#Test IDBA
```
#install IDBA from https://github.com/loneknightpy/idba
export PATH="$MOPT/idba-1.1.3/bin/":$PATH
cd APR842_00/AS1M/Assemble_mashBC/BC000006/
#convert
fq2fa --merge --filter sort.1.fq sort.2.fq sort.pair.fa
#run
idba_ud -o idbaud -r sort.pair.fa --mink 21 --step 22 --maxk 121 --min_contig 999
blastn -num_threads 8 -db $refBT -query idbaud/scaffold.fa -out idbaud/scaffolds.F.BLAST.tax.blast6 -outfmt 6
mv APR842_00/AS1M/{Assemble_mashBC,run1.Assemble_mashBC}
```

# Run2
use `SUB_2R100_0D05` to test
```bash
for i in `ls APR842_00/TT1M/SUB_2R100_0D05.Assemble_mashBC/BC000*/sort.{1,2}.fq`;do t=${i/TT1M/AS1M};b=`dirname $t`;mkdir -p $b;cp -v $i $t;done
ln -sf ../TT1M/SUB_2R100_0D05.mash APR842_00/AS1M/mash
ln -sf SUB_2R100_0D05.Assemble_mashBC APR842_00/AS1M/Assemble_mashBC
snakemake -s test.smk -j -np $SAM0/AS1M/batch.assemble.BC.sh
#post
mv APR842_00/AS1M/{Assemble_mashBC,run2.Assemble_mashBC}
#summary
mode="idba";for i in `ls APR842_00/AS1M/run2.Assemble_mashBC/`;do echo $i;awk '$4>999{print}' APR842_00/AS1M/run2.Assemble_mashBC/$i/$mode/scaffolds.$mode.BLAST.tax.blast6.anno.best;done |column -t|les
```

######tmp
```bash
spades-gbuilder input_dataset.yaml testGB/graph -k 77 -t 24  -gfa
```

# Run3

use `ASM_2R50_0D05` to test

```bash
tag="ASM_2R50_0D10"
mkdir $SAM0/$tag
ln -s ../clean $SAM0/$tag/clean
snakemake -s test.smk -j -np $SAM0/$tag/mash/bMin2.raw.dist.anno.stat #once
snakemake -s test.smk --forceall -j -np $SAM0/$tag/batch.assemble.BC.sh
#post
mv APR842_00/AS1M/{Assemble_mashBC,run3.Assemble_mashBC}
#summary
mode="megahit";for i in `ls APR842_00/AS1M/Assemble_mashBC/`;do echo $i;awk '$4>999{print}' APR842_00/AS1M/Assemble_mashBC/$i/$mode/scaffolds.$mode.BLAST.tax.blast6.anno.best;done |column -t|les
mode="megahit";for i in `ls APR842_00/AS1M/run3.Assemble_mashBC/`;do quast.py -i APR842_00/AS1M/run3.Assemble_mashBC/$i/$mode/scaffolds.fa -o APR842_00/AS1M/run3.Assemble_mashBC/$i/$mode/quast;done
```

###### 