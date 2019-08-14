
```bash
ln -s $mSeq/instance/SSB831 ../../Results/
cd ../../Results/SSB831/clean
metabbq randomlyPick.pl -i fastp.sort.1.fq -o fastp.sort.1 -n 10 -m 3 -s 831 -t 4 &
metabbq randomlyPick.pl -i fastp.sort.2.fq -o fastp.sort.2 -n 10 -m 3 -s 831 -t 4 -v
#link to SSB831_00 dir

```
```bash
#stat and check the pre-mash beads number
metabbq smk -p SSB831_00/clean/BB.stat
metabbq smk -np SSB831_00/mash/lv1/tree.cluster.count.main
metabbq beadsWrite3.pl -x --r1 SSB831_00/mash/bMin2.sort.1.fq
metabbq smk -np SSB831_00/summary.megahit.contig.fasta
```
