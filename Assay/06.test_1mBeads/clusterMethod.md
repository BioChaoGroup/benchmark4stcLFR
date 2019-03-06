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
pigz -p 8 -dc ../../Results/$SAM/clean/$SAM.fps.1_00.fq.gz > $SAM0/clean/fastp.sort.1.fq &
pigz -p 8 -dc ../../Results/$SAM/clean/$SAM.fps.2_00.fq.gz > $SAM0/clean/fastp.sort.2.fq &
```
The default parameters:
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
### Run 01
params: `p_rpc_min:100`
```bash
SAM0=APR831\_00
# seq index
metabbq smk -j -np $SAM0/clean/BB.stat
# main
snakemake -s test.smk --config p_rpc_min=100 \
 -j -np $SAM0/TT1M/batch.assemble.BC.sh
# post
mv $SAM0/TT1M/{mash,SUB_2R20_0D04.mash}
mv $SAM0/TT1M/{Assemble_mashBC,SUB_2R20_0D04.Assemble_mashBC}

```
