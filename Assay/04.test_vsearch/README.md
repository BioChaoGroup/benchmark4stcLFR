#Test pipeline based on VSEARCH

### 1 The information of dereplicate
The dereplicate process provide a uc file (unchimera) which list the alignment
of each 2 sequences. Since the CIGAR string always be '*', I guess they are 100%
matched.
Let me stat their overlap between beads.
(The experimental work dir is `wdir=instance/DEC560_d10_00/VSEARCH`)

```bash
cd $wdir
perl beadStat.pl -m uc -i read.merge.derep.uc -o read.merge.derep.uc.bc -v
#grep "^T" read.merge.derep.uc.bc | sort -nrk4,4 > read.merge.derep.uc.bc.s4

perl beadStat.pl -m uc -ident 100 -match 50 -i read.merge.derep.preclustered.uc -o \
read.merge.derep.preclustered.uc.bc -v
#grep "^T" read.merge.derep.preclustered.uc.bc | sort -nrk4,4 > \ read.merge.derep.preclustered.uc.bc.s4

# cat read.merge.derep.uc.bc.s4 read.merge.derep.preclustered.uc.bc.s4| \
# perl beadStat.pl -m tile -o read.merge.derep.bc.clusters
cat read.merge.derep.uc read.merge.derep.preclustered.uc| perl beadStat.pl -m uc -ident 100 -match 50 -o read.merge.derep.2.bc -v

#If it's too huge to run at once:
perl -e '{open I1,"<read.merge.derep.uc.bc.v2";open I2,"<read.merge.derep.preclustered.uc.bc"}'
#cat read.merge.derep.uc read.merge.derep.preclustered.uc| perl beadStat.pl -m uc -o read.merge.derep.2.bc -v
awk '$4>1{print $0}' read.merge.derep.2.bc | sort -nrk4,4 > read.merge.derep.2T.bc.s4
```

prepare clean reads
```bash
cd ../clean/
pigz -dc fastp.sort.1.fq.gz > fastp.sort.1.fq
pigz -dc fastp.sort.2.fq.gz > fastp.sort.2.fq
perl ../../src/beadsWrite3.pl --r1 fastp.sort.1.fq --r2 fastp.sort.2.fq -x -v

perl ../../src/beadsWrite3.pl --r1 fastp.sort.1.fq --r2 fastp.sort.2.fq \
-b ../VSEARCH/read.merge.derep.2T.bc.cluster -o ../Assemble -f fq -v

cd ..
```

Cluster:
```bash
awk '$4>1{print $2"\t"$3"\t"100/$4/$5}' read.merge.derep.2T.bc.s4 > read.merge.derep.2T.bc.dist
perl ../../src/filterMASHoutput.pl -i read.merge.derep.2T.bc.dist -o read.merge.derep.2T.bc.reNum.dist -M read.merge.derep.2T.bc.map
convert -i read.merge.derep.2T.bc.reNum.dist -w read.merge.derep.2T.bc.w -o read.merge.derep.2T.bc.bin
community read.merge.derep.2T.bc.bin -w read.merge.derep.2T.bc.w -l -1 > read.merge.derep.2T.bc.graph.tree
maxLv=`hierarchy read.merge.derep.2T.bc.graph.tree -n|head -1|sed 's/Number of levels: //'`
((maxLv=$maxLv -1 ))
hierarchy read.merge.derep.2T.bc.graph.tree -l $maxLv > read.merge.derep.2T.bc.node\_Lv$maxLv.lst
perl ../../src/change.id.pl -n read.merge.derep.2T.bc.node\_Lv$maxLv.lst \
-m read.merge.derep.2T.bc.map -v -o read.merge.derep.2T.bc.cluster\_Lv$maxLv

cut -d " " -f2 read.merge.derep.2T.bc.node\_Lv$maxLv.lst|sort|uniq -c|sort -nr > read.merge.derep.2T.bc.node\_Lv$maxLv.rank

perl -e 'open IN,"<../clean/fastp.sort.1.fq.idx";while(<IN>){chomp;@a=split;$num=($a[2]-$a[1]+1)/4;$HASH{$a[0]}=$num};while(<>){chomp;@a=split;$HB{$a[0]}{R}+=$HASH{$a[1]};$HB{$a[0]}{C}++};foreach my $c (sort {$a<=>$b} keys %HB){print "$c\t$HB{$c}{C}\t$HB{$c}{R}\n"}' <  read.merge.derep.2T.bc.cluster\_Lv$maxLv >  \
read.merge.derep.2T.bc.cluster\_Lv$maxLv.count

sort -k3,3nr read.merge.derep.2T.bc.cluster\_Lv$maxLv.count| \
awk 'FNR==1{b=$2;r=$3}($2>b/10 && $3>r/10 && FNR<300){b=$2;r=$3;print $0}'|sort -k1,1n \
> read.merge.derep.2T.bc.cluster\_Lv$maxLv.count.main
perl -e 'open IN,"sort -k1,1n read.merge.derep.2T.bc.cluster\_Lv'$maxLv'.count.main|";while(<IN>){@a=split(/\t/,$_);push @ids, $a[0]}; $tag= shift @ids; while(<>){@a=split /\t/, $_; if($a[0] < $tag){next}elsif($a[0] == $tag){print $_}else{$tag= shift @ids;print $_ if $a[0] == $tag } }' < read.merge.derep.2T.bc.cluster\_Lv$maxLv |sort -k1,1n > read.merge.derep.2T.bc.cluster\_Lv$maxLv.main

cd ../
mkdir -p Assemble\_Lv$maxLv
perl ../src/beadsWrite3.pl --r1 ./clean/fastp.sort.1.fq --r2 ./clean/fastp.sort.2.fq \
-b VSEARCH/read.merge.derep.2T.bc.cluster\_Lv$maxLv.main -o Assemble\_Lv$maxLv -f fq -v
```

kaiju
```
kaijudb=/ldfssz1/ST_META/share/database/kaijudb
kaiju -t $kaijudb/nodes.dmp -f $kaijudb/kaiju_db.fmi -i sort.1.fq -j sort.2.fq > kaiju.out
addTaxonNames -t $kaijudb/nodes.dmp -n $kaijudb/names.dmp -i kaiju.out -o kaiju.names.out

```


#Silva
```
cd REF/silva132
zcat source/SILVA_132_LSURef_tax_silva.fasta.gz|grep "^>"|awk -F ';' '{sub(">","");split($1,a," ");gsub(" ","_");if($7==""){$7="unidentified"};print $7"\t"a[1]"\tLSU\tRef\tK__"a[2]";p__"$2";c__"$3";o__"$4";f__"$5";g__"$6";s__"$7}' > SILVA_132_LSURef_tax_silva.fasta.ids

zcat source/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz|grep "^>"|awk -F ';' '{sub(">","");split($1,a," ");gsub(" ","_");if($7==""){$7="unidentified"};print $7"\t"a[1]"\tSSU\tRef\tK__"a[2]";p__"$2";c__"$3";o__"$4";f__"$5";g__"$6";s__"$7}' > SILVA_132_SSURef_Nr99_tax_silva.fasta.ids

```
