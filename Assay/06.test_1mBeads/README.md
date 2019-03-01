# Test per 1 million beads
Key questions:
1. Charactrastic of bead-reads distribution.
2. How many beads with hybrided DNA content.
3. Method to refined the monocolonal beads for assembly

# prepare
```bash
cp $mSeq2/config.yaml ./ # and modified for test
cp $mSeq2/rules/beadCluster_1m.smk test.smk # then modified for test
ln -s $mSeq2/instance/REF

# data source
ln -s ../../Results/APR842_00_00
ln -s ../../Results/APR841_00_00
ln -s ../../Results/DEC560_d10_00
```

**Visualization** code are prepared in `bead1mTest.Rmd`.

## show the beads distribution
```bash
SAM=APR842_00_00

```

# snakemake pipeline execution

### 1st run: beadreads=10 (dist: no selection)
```bash
SAM=APR842_00_00 # or other 2 sample
snakemake -s test.smk --config p_cluster_maxR=10 -j -np \
$SAM/TT1M/batch.assemble.BC.sh # remove -n to execute
#post move
mv $SAM/TT1M/{mash,BeadReads_10_mash}
mv $SAM/TT1M/{Assemble_mashBC,BeadReads_10_Assemble_mashBC}
```

### 2nd run: beadreads=10 && dist <= 0.06
```bash
mkdir -p $SAM/TT1M/mash
ln -s ../BeadReads_10_mash/bMin2.raw.dist  $SAM/TT1M/mash/
awk '$3>0&&$3<=0.06{print}' $SAM/TT1M/mash/bMin2.raw.dist \
> $SAM/TT1M/mash/bMin2.clean.dist
snakemake -s test.smk -j -np $SAM/TT1M/batch.assemble.BC.sh
mv $SAM/TT1M/{mash,BeadReads_10_0D06_mash}
mv $SAM/TT1M/{Assemble_mashBC,BeadReads_10_0D06_Assemble_mashBC}
```

### 3rd run: beadreads=10 && dist <= 0.04
```bash
# third run
mkdir -p $SAM/TT1M/mash
ln -s ../BeadReads_10_mash/bMin2.raw.dist  $SAM/TT1M/mash/
awk '$3>0&&$3<=0.04{print}' $SAM/TT1M/mash/bMin2.raw.dist \
> $SAM/TT1M/mash/bMin2.clean.dist
snakemake -s test.smk -j -np $SAM/TT1M/batch.assemble.BC.sh
mv $SAM/TT1M/{mash,BeadReads_10_0D04_mash}
mv $SAM/TT1M/{Assemble_mashBC,BeadReads_10_0D04_Assemble_mashBC}

snakemake -s test.smk --config p_cluster_maxR=5 p_cluster_minR=5 \
p_cluster_topB=10000 -j -np $SAM/TT1M/batch.assemble.BC.sh
```

### 4th run: beadreads=10-20 && dist <= 0.02 && clusterReads >= 100
```bash
# fourth run
snakemake -s test.smk -j -np $SAM/TT1M/mash/bMin2.raw.dist
awk '$3>0&&$3<=0.02{print}' $SAM/TT1M/mash/bMin2.raw.dist \
> $SAM/TT1M/mash/bMin2.clean.dist
snakemake -s test.smk -j -np $SAM/TT1M/mash/bMin2.bc.tree.target.cluster.count
#keep top rich communities
pfx=$SAM/TT1M/mash/bMin2
awk '($3>=100){print}' $pfx.bc.tree.target.cluster.count |sort -k1,1n \
> $pfx.bc.tree.target.cluster.count.main
perl -e 'open IN,"sort -k1,1n '$pfx'.bc.tree.target.cluster.count.main|";
while(<IN>){@a=split(/\t/,$_);push @ids, $a[0]}; $tag= shift @ids;
while(<>){@a=split /\t/, $_; if($a[0] < $tag){next}elsif($a[0] == $tag){
  print $_}else{$tag= shift @ids;print $_ if $a[0] == $tag } }' $pfx.bc.tree.target.cluster |sort -k1,1n > $pfx.bc.tree.target.cluster.main
ln -sf bMin2.bc.tree.target.cluster.main $pfx.bc.tree.target.cluster

snakemake -s test.smk -j -np $SAM/TT1M/batch.assemble.BC.sh
mv $SAM/TT1M/{mash,BeadReads_10R20_0D02_mash}
mv $SAM/TT1M/{Assemble_mashBC,BeadReads_10R20_0D02_Assemble_mashBC}
```

### 5th run: beadreads=10-20 && dist <= 0.04 && clusterReads >= 100
```bash
mkdir -p $SAM/TT1M/mash
ln -sf ../BeadReads_10R20_0D02_mash/bMin2.raw.dist  $SAM/TT1M/mash/bMin2.raw.dist
awk '$3>0&&$3<=0.04{print}' $SAM/TT1M/mash/bMin2.raw.dist \
> $SAM/TT1M/mash/bMin2.clean.dist
snakemake -s test.smk -j -np $SAM/TT1M/mash/bMin2.bc.tree.target.cluster.count
#keep top rich communities
pfx=$SAM/TT1M/mash/bMin2
awk '($3>=100){print}' $pfx.bc.tree.target.cluster.count |sort -k1,1n \
> $pfx.bc.tree.target.cluster.count.main
perl -e 'open IN,"sort -k1,1n '$pfx'.bc.tree.target.cluster.count.main|";
while(<IN>){@a=split(/\t/,$_);push @ids, $a[0]}; $tag= shift @ids;
while(<>){@a=split /\t/, $_; if($a[0] < $tag){next}elsif($a[0] == $tag){
  print $_}else{$tag= shift @ids;print $_ if $a[0] == $tag } }' $pfx.bc.tree.target.cluster |sort -k1,1n > $pfx.bc.tree.target.cluster.main
ln -sf bMin2.bc.tree.target.cluster.main $pfx.bc.tree.target.cluster

snakemake -s test.smk -j -np $SAM/TT1M/batch.assemble.BC.sh
mv $SAM/TT1M/{mash,BeadReads_10R20_0D02_mash}
mv $SAM/TT1M/{Assemble_mashBC,BeadReads_10R20_0D02_Assemble_mashBC}
```

### 6th run: beadreads=2-20 && dist <= 0.04 && randomly pick 10000 beads
```bash
SAM=APR841_00_00
SAM=APR842_00_00
snakemake -s test.smk --config p_cluster_maxR=20 p_cluster_minR=2 \
p_cluster_topB=10000 -j -np $SAM/TT1M/batch.assemble.BC.sh
mv $SAM/TT1M/{mash,SUB_2R20_0D04_S10k.mash}
mv $SAM/TT1M/{Assemble_mashBC,SUB_2R20_0D04_S10k.Assemble_mashBC}
```

### 7th run: beadreads=2-20 && dist <= 0.04 && randomly pick 10000 beads
```bash
SAM=APR841_00_00
SAM=APR842_00_00
mkdir -p $SAM/TT1M/mash
ln -s ../SUB_2R20_0D04_S10k.mash/bMin2.sort.1.fq $SAM/TT1M/mash/bMin2.sort.1.fq
ln -s ../SUB_2R20_0D04_S10k.mash/bMin2.sort.2.fq $SAM/TT1M/mash/bMin2.sort.2.fq
ln -s ../SUB_2R20_0D04_S10k.mash/bMin2.raw.dist $SAM/TT1M/mash/bMin2.raw.dist
snakemake -s test.smk --config p_cluster_maxR=20 p_cluster_minR=2 \
p_cluster_topB=10000 -j -np $SAM/TT1M/batch.assemble.BC.sh
mv $SAM/TT1M/{mash,SUB_2R20_0D04_S10k.mash}
mv $SAM/TT1M/{Assemble_mashBC,SUB_2R20_0D04_S10k.Assemble_mashBC}
```

### 9th test another sampling APR842
```bash

SAM0=APR842_00
mkdir -p ../../Results/$SAM0/clean
ln -s ../../Results/$SAM0
pigz -p 8 -dc $mSeq2/instance/clean/APR842_00.fps.1.fq.gz > $SAM0/clean/fastp.sort.1.fq &
pigz -p 8 -dc $mSeq2/instance/clean/APR842_00.fps.2.fq.gz > $SAM0/clean/fastp.sort.2.fq &

```
### 9th prepare APR831
> run: beadreads=2-20 && dist <= 0.04 without randomly pick

```bash
SAM=APR831
mkdir -p ../../Results/$SAM/clean
# randomly pick sorted fastq
metabbq randomlyPick.pl -t 8 -n 100 -m 3 -s 831 -i \
$mSeq2/instance/clean/$SAM.fps.1.fq.gz -o ../../Results/$SAM/clean/$SAM.fps.1 &
metabbq randomlyPick.pl -t 8 -n 100 -m 3 -s 831 -i \
$mSeq2/instance/clean/$SAM.fps.2.fq.gz -o ../../Results/$SAM/clean/$SAM.fps.2 &
wait
SAM0=$SAM\_00
mkdir -p ../../Results/$SAM0/clean
ln -s ../../Results/$SAM0
pigz -p 8 -dc ../../Results/$SAM/clean/APR831.fps.1_00.fq.gz > $SAM0/clean/fastp.sort.1.fq &
pigz -p 8 -dc ../../Results/$SAM/clean/APR831.fps.2_00.fq.gz > $SAM0/clean/fastp.sort.2.fq &

# seq index
metabbq smk -j -np APR831_00/clean/BB.stat

snakemake -s test.smk --config p_cluster_maxR=20 p_cluster_minR=2 \
 -j -np $SAM0/TT1M/batch.assemble.BC.sh

pfx=$SAM0/TT1M/mash/bMin2
awk '($3>=100){print}' $pfx.bc.tree.target.cluster.count |sort -k1,1n \
> $pfx.bc.tree.target.cluster.count.main
perl -e 'open IN,"sort -k1,1n '$pfx'.bc.tree.target.cluster.count.main|";
while(<IN>){@a=split(/\t/,$_);push @ids, $a[0]}; $tag= shift @ids;
while(<>){@a=split /\t/, $_; if($a[0] < $tag){next}elsif($a[0] == $tag){
 print $_}else{$tag= shift @ids;print $_ if $a[0] == $tag } }' $pfx.bc.tree.target.cluster |sort -k1,1n > $pfx.bc.tree.target.cluster.main
ln -sf bMin2.bc.tree.target.cluster.main $pfx.bc.tree.target.cluster

snakemake -s test.smk --config p_cluster_maxR=20 p_cluster_minR=2 \
 -j -np $SAM0/TT1M/batch.assemble.BC.sh


mv $SAM0/TT1M/{mash,SUB_2R20_0D04.mash}
mv $SAM0/TT1M/{Assemble_mashBC,SUB_2R20_0D04.Assemble_mashBC}
```

### 10th APR831:
> run: beadreads=2-20 && dist <= 0.1 && reads per cluster >=100

```bash
SAM0=APR831_00
mkdir -p $SAM0/TT1M/mash
ln -s ../SUB_2R20_0D04.mash/bMin2.sort.1.fq $SAM0/TT1M/mash/bMin2.sort.1.fq
ln -s ../SUB_2R20_0D04.mash/bMin2.sort.2.fq $SAM0/TT1M/mash/bMin2.sort.2.fq
ln -s ../SUB_2R20_0D04.mash/bMin2.raw.dist $SAM0/TT1M/mash/bMin2.raw.dist


snakemake -s test.smk --config p_dist_max=0.1 p_rpc_min=100\
 -j -np $SAM0/TT1M/batch.assemble.BC.sh


mv $SAM0/TT1M/{mash,SUB_2R20_0D10.mash}
mv $SAM0/TT1M/{Assemble_mashBC,SUB_2R20_0D10.Assemble_mashBC}
```


# The mash distance distribution.
```bash
awk '{print int($3*1000+0.5)/1000}' $SAM/TT1M/BeadReads_10/mash/bMin2.clean.dist |
sort|uniq -c > $SAM/TT1M/BeadReads_10.dist.stat
```

### test
```
#cd into BC dir
mkdir -p align
blastn -db $REF/fungal5REF.fa -query scaffolds.fasta -out align/scaf2REF.blast6 -outfmt 6
metabbq binWrite best -u -m -i align/scaf2REF.blast6
column -t align/scaf2REF.blast6.best
cut -f1-13 scaffolds.F.BLAST.tax.blast6.anno.best


```
