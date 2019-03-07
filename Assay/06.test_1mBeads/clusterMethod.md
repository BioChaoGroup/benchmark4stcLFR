# Test per 1 million beads
Key questions:
1. Charactrastic of bead-reads distribution.
2. How many beads with hybrided DNA content.
3. Method to refined the monocolonal beads for assembly

# Prepare
First link database and copy config as template:
```bash
cp $mSeq2/config.yaml ./ # and modified for test
cp $mSeq2/rules/beadCluster_1m.smk test.smk # then modified for test
ln -s $mSeq2/instance/REF
```
Specify data source to test:
```bash
#SAM=APR831 # processing now
#SAM=APR832 # Haven't been tested
#SAM=APR841 # Needs update
#SAM=APR842 # Needs update
```
For the first time, randomly pick about 10% beads for test:
```bash
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
pigz -p 8 -dc ../../Results/$SAM/clean/$SAM.fps.1_00.fq.gz > $SAM0/clean/fastp.sort.1.fq &
pigz -p 8 -dc ../../Results/$SAM/clean/$SAM.fps.2_00.fq.gz > $SAM0/clean/fastp.sort.2.fq &
```
Define the default testing parameters:
```yaml
### params
p_cluster_minR: 2
p_cluster_maxR: 20
p_cluster_topB: 0
p_cluster_ranP: 0

p_dist_min: 0
p_dist_max: 0.04

#reads per cluster threashold
p_rpc_min: 0
```
# Run
Each run may conained 3 steps:
1. pre-command, to generate essential data for each sample.
2. command to test.  Mainly organized by `snakemake`. The command provide here are tagged by `-n` for "dry-run" mode, which only preview the pipeline will be executed. To actually run it, remove the `-n` tag.
3. post-command. Since it's a test for a specific step, results are generated in a commonly named directory. To avoid be overwritten by next test, it needs to be moved to another position.

### Run 01
Reproduce the results from Song.  
Key params: `p_dist_max: 0.04`(same with default).  
Optional params: `p_rpc_min=100`.
```bash
#results tag
rTag=SUB_2R20_0D04
# seq index
metabbq smk -j -np $SAM0/clean/BB.stat
# main
snakemake -s test.smk --config p_rpc_min=100 \
 -j -np $SAM0/TT1M/batch.assemble.BC.sh
# post
mv $SAM0/TT1M/{mash,$rTag.mash}
mv $SAM0/TT1M/{Assemble_mashBC,$rTag.Assemble_mashBC}
```

### Run 02
Test the distance range `(0,0.1)`.  
Key params: `p_dist_max: 0.1`.  
Optional params: `p_rpc_min=100`.
```bash
#results tag
rTag=SUB_2R20_0D10
# link previously results
mkdir -p $SAM0/TT1M/mash
ln -s ../SUB_2R20_0D04.mash/bMin2.raw.dist $SAM0/TT1M/mash/bMin2.raw.dist
# main
snakemake -s test.smk --config p_dist_max=0.1 p_rpc_min=100  -j -np $SAM0/TT1M/batch.assemble.BC.sh
# post
mv $SAM0/TT1M/{mash,$rTag.mash}
mv $SAM0/TT1M/{Assemble_mashBC,$rTag.Assemble_mashBC}
# print annotation
for i in `ls -d $SAM0/TT1M/$rTag.Assemble_mashBC/*`;do echo $i;awk '$4>999{print}'t $i/scaffolds.F.BLAST.tax.blast6.anno.best|column -t;done
```
Though long scaffolds increased, the hybridization also increased.

### Run 03
Test the distance range `(0,0.1)` and rpb up to 100.
Key params: `p_dist_max: 0.1` & `p_cluster_maxR: 100`.  
Optional params: none.
```bash
#results tag
rTag=SUB_2R100_0D10
# link previously results
# main
snakemake -s test.smk --config p_dist_max=0.1 p_cluster_maxR=100  -j -np $SAM0/TT1M/mash/bMin2.bc.tree.target.cluster.main
# post
mv $SAM0/TT1M/{mash,$rTag.mash}
```

### Run 03.1
Re-filter the distance range `(0,0.2)` from distance in #Run03.
```bash
#results tag
rTag=SUB_2R100_0D20
# link previously results
mkdir $SAM0/TT1M/mash
ln -s ../SUB_2R100_0D10.mash/bMin2.raw.dist $SAM0/TT1M/mash/
# main
snakemake -s test.smk --config p_dist_max=0.2 p_cluster_maxR=100  -j -np $SAM0/TT1M/mash/bMin2.bc.tree.target.cluster.main
# post
mv $SAM0/TT1M/{mash,$rTag.mash}
```
