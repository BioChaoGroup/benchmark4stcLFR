# Alignment assessment for mock sample

### Get subunits of rRNA for Zymo reference seq
```bash
#bac
barrnap --threads 4 --kingdom bac < $LFR/Source/REF/zymo/D6305.rRNA.fa > $LFR/Source/REF/zymo/D6305.rRNA.barrnap.gff 2> $LFR/Source/REF/zymo/D6305.rRNA.barrnap.log
awk -F "\t|;" '($0!~/^#/){ gsub("Name=", "", $9); \
if($1==arr[1] && $7==arr[7] && (($9 == "16S_rRNA" && arr[9] == "23S_rRNA")|| ($9 == "23S_rRNA" && arr[9] == "16S_rRNA"))){ print $1"\t"arr[5]+1"\t"$4-1"\tITS_rRNA\t1\t"$7"\t"arr[5]+1"\t"$4-1"\t50,50,05";}	\
print $1"\t"$4"\t"$5"\t"$9"\t1\t"$7"\t"$4"\t"$5"\t50,50,05"; \
arr[1]=$1;arr[4]=$4;arr[5]=$5;arr[9]=$9;arr[7]=$7;}' \
< $LFR/Source/REF/zymo/D6305.rRNA.barrnap.gff > $LFR/Source/REF/zymo/D6305.rRNA.barrnap.bed
barrnap --threads 4 --kingdom euk < $LFR/Source/REF/zymo/D6305.rRNA.fa > $LFR/Source/REF/zymo/D6305.rRNA.barrnap.euk.gff 2> $LFR/Source/REF/zymo/D6305.rRNA.barrnap.euk.log

#fungi
makeblastdb -dbtype nucl -in $LFR/Source/REF/fungi/sanger/sanger.mock7.rRNA.fa
barrnap --threads 4 --kingdom euk < $LFR/Source/REF/fungi/sanger/sanger.mock7.rRNA.fa \
> $LFR/Source/REF/fungi/sanger/sanger.mock7.rRNA.barrnap.gff 2> $LFR/Source/REF/fungi/sanger/sanger.mock7.rRNA.barrnap.log
ITSx -t F --cpu 4 -i $LFR/Source/REF/fungi/sanger/sanger.mock7.rRNA.fa -o $LFR/Source/REF/fungi/sanger/sanger.mock7.rRNA.ITSx \
 &> $LFR/Source/REF/fungi/sanger/sanger.mock7.rRNA.ITSx.log
 awk -F " |\t|-" '{ \
 print $1"\t"$5-1"\t"$18"\tALL\t1\t+\t"$5"\t"$18"\t50,50,05";
 print $1"\t"$5"\t"$6"\t"$4"\t1\t+\t"$5"\t"$6"\t50,50,05"; \
 print $1"\t"$8"\t"$9"\t"$7"\t1\t+\t"$8"\t"$9"\t50,50,05"; \
 print $1"\t"$11"\t"$12"\t"$10"\t1\t+\t"$11"\t"$12"\t50,50,05"; \
 print $1"\t"$14"\t"$15"\t"$13"\t1\t+\t"$14"\t"$15"\t50,50,05"; \
 print $1"\t"$17"\t"$18"\t"$16"\t1\t+\t"$17"\t"$18"\t50,50,05"; \
}' $LFR/Source/REF/fungi/sanger/sanger.mock7.rRNA.ITSx.positions.txt |sed 's/://' \
> $LFR/Source/REF/fungi/sanger/sanger.mock7.rRNA.ITSx.bed

```

**#Get subunits of zymo without potential deletion region
```bash
awk '{if(pn!=$1&&b!=""){print b"-"pp;b=$1":"$2}else if(b==""&&$3/(pc+0.1)>10){b=$1":"$2}else if(b!=""&&($3+0.1)/(pc+0.1)<0.1){print b"-"$2-1;b=""};pn=$1;pp=$2;pc=$3}' SAM/Z3/CLIP/clip2zymo.ref.cov > SAM/Z3/CLIP/clip2zymo.ref.region.txt
```
get (with modified):
```
Bacillus_subtilis:9-4516
Enterococcus_faecalis:1-4260
Escherichia_coli:7-1607
Escherichia_coli:1709-4748
Lactobacillus_fermentum:313-4869
Listeria_monocytogenes:1-4528
Pseudomonas_aeruginosa:402-5152
Salmonella_enterica:8-2436
Salmonella_enterica:2550-4764
Staphylococcus_aureus:313-5118
```
```bash
samtools faidx -r SAM/Z3/CLIP/clip2zymo.ref.region.txt $LFR/Source/REF/zymo/D6305.rRNA.fa > $LFR/Source/REF/zymo/D6305.fix.rRNA.fa
#join segments together
barrnap --threads 4 --kingdom bac < $LFR/Source/REF/zymo/D6305.fix.rRNA.fa > $LFR/Source/REF/zymo/D6305.fix.rRNA.barrnap.gff 2> $LFR/Source/REF/zymo/D6305.fix.rRNA.barrnap.log
awk -F "\t|;" '($0!~/^#/){ gsub("Name=", "", $9); \
if($1==arr[1] && $7==arr[7] && (($9 == "16S_rRNA" && arr[9] == "23S_rRNA")|| ($9 == "23S_rRNA" && arr[9] == "16S_rRNA"))){ print $1"\t"arr[5]+1"\t"$4-1"\tITS_rRNA\t1\t"$7"\t"arr[5]+1"\t"$4-1"\t50,50,05";}	\
print $1"\t"$4"\t"$5"\t"$9"\t1\t"$7"\t"$4"\t"$5"\t50,50,05"; \
arr[1]=$1;arr[4]=$4;arr[5]=$5;arr[9]=$9;arr[7]=$7;}' \
< $LFR/Source/REF/zymo/D6305.fix.rRNA.barrnap.gff > $LFR/Source/REF/zymo/D6305.fix.rRNA.barrnap.bed


```

```R
library(dplyr)
library(plyr)
m6 <- read.table("SAM/Z1/summary.BI.megahit.clip2zymo.1.m6")
uc <- read.table("SAM/Z1/CLIP/preclust.uc")
uc.c <- uc%>%filter(V1=="C")
uc.h <- uc%>%filter(V1=="H")
uc.c$V10 <- uc.c$V9
uc.df <- rbind(uc.c,uc.h)
uc.m6 <- merge(uc.df,m6,by.x="V9",by.y="V1")

uc.consist <- ddply(uc.m6,"V10.x",summarise,refs=length(unique(V2.y)),count=length(V2.y))
```


### test checkBeadRegion
```bash
metabbq checkBeadRegion -i SAM/Z1/summary.BI.megahit.clip2zymo.m6 -r $LFR/Source/REF/zymo/D6305.rRNA.barrnap.bed -o SAM/Z1/CLIP/clip2zymo.bead.anno -c SAM/Z1/CLIP/clip2zymo.clip.anno -v

perl -d /ldfssz1/ST_META/F16HQSB1SY2636_TumorMeta_ZH/fangchao/metaSeq/util/checkBeadRegion -i SAM/Z1/summary.BI.megahit.clip2zymo.m6 -r $LFR/Source/REF/zymo/D6305.rRNA.barrnap.bed -o SAM/Z1/CLIP/clip2zymo.bead.anno -c SAM/Z1/CLIP/clip2zymo.clip.anno -g SAM/Z1/CLIP/clip2zymo.ref.cov -v

#for fix reference
makeblastdb -dbtype nucl -in $LFR/Source/REF/zymo/D6305.fix.rRNA.fa
blastn -num_threads 8 -db $LFR/Source/REF/zymo/D6305.fix.rRNA.fa -query SAM/Z1/summary.BI.megahit.clip.fasta -outfmt '6 std qlen' -out SAM/Z1/summary.BI.megahit.clip2zymo.fix.m6
metabbq checkBeadRegion -i SAM/Z1/summary.BI.megahit.clip2zymo.fix.m6 -r $LFR/Source/REF/zymo/D6305.fix.rRNA.barrnap.bed -o SAM/Z1/CLIP/clip2zymo.fix.bead.anno -c SAM/Z1/CLIP/clip2zymo.fix.clip.anno -v

for i in `sort -k2,2nr SAM/Z1/CLIP/Ecoli.bead.anno|head -500|cut -f1`;do grep -A 1 $i SAM/Z1/summary.BI.megahit.clip.fasta;
done > SAM/Z1/CLIP/fix.top500.clip.fasta
for i in `sort -k2,2nr SAM/Z1/CLIP/Ecoli.bead.anno|head -500|cut -f1`;
do cat SAM/Z1/Assemble_BI/$i/sort.1.fq|awk 'FNR%4==1{gsub("/1$","");gsub("/","/1\tBX:Z:")}{print}';
done > SAM/Z1/CLIP/fix.top500.1.fq &
for i in `sort -k2,2nr SAM/Z1/CLIP/Ecoli.bead.anno|head -500|cut -f1`;
do cat SAM/Z1/Assemble_BI/$i/sort.2.fq|awk 'FNR%4==1{gsub("/2$","");gsub("/","/2\tBX:Z:")}{print}';
done > SAM/Z1/CLIP/fix.top500.2.fq

quast.py -t 16 SAM/Z1/CLIP/fix.top500.clip.fasta -R $LFR/Source/REF/zymo/D6305.fix.rRNA.fa  -o SAM/Z1/CLIP/clip2zymo.fix_500_quast \
-1 SAM/Z1/CLIP/fix.top500.1.fq -2 SAM/Z1/CLIP/fix.top500.2.fq
```

### batch generation
```bash
metabbq smk -j -npk SAM/Z{1,2,3}/CLIP/clip2zymo.bead.anno

mkdir -p STAT/CLIP
for i in Z{1,2,3};do
	awk -v ID=$i '{print ID"\t"$0}' SAM/$i/CLIP/clip2zymo.bead.anno;
done > STAT/CLIP/clip2zymo.bead.anno
for i in Z{1,2,3};do
	awk -v ID=$i '{print ID"\t"$0}' SAM/$i/CLIP/clip2zymo.ref.cov;
done > STAT/CLIP/clip2zymo.ref.cov

metabbq smk --configfile F.config.yaml -j -npk SAM/M{4,5,6}/CLIP/clip2CloseRef.bead.anno
for i in M{4,5,6};do awk -v ID=$i '{print ID"\t"$0}' SAM/$i/CLIP/clip2CloseRef.bead.anno;
done > STAT/CLIP/F.M.clip2CloseRef.bead.anno
for i in M{4,5,6};do awk -v ID=$i '{print ID"\t"$0}' SAM/$i/CLIP/clip2CloseRef.ref.cov;
done > STAT/CLIP/F.M.clip2CloseRef.ref.cov

#for fix ref
metabbq smk -j -npk SAM/Z{1,2,3}/CLIP/clip2zymo.fix.bead.anno
for i in Z{1,2,3};do
	awk -v ID=$i '{print ID"\t"$0}' SAM/$i/CLIP/clip2CloseRef.fix.bead.anno;
done > STAT/CLIP/clip2CloseRef.fix.bead.anno
for i in Z{1,2,3};do
	awk -v ID=$i '{print ID"\t"$0}' SAM/$i/CLIP/clip2CloseRef.fix.ref.cov;
done > STAT/CLIP/clip2CloseRef.fix.ref.cov

```

### cal beads' contigs' kmers
```bash
cat SAM/Z1/Assemble_BI/*/sort.1.fq > SAM/Z1/mash/Asm.sort.1.fq &
cat SAM/Z1/Assemble_BI/*/sort.2.fq > SAM/Z1/mash/Asm.sort.2.fq
mash sketch -p 8 -k 32 -s 2000 -r -B SAM/Z1/mash/Asm.sort.{1,2}.fq -o SAM/Z1/mash/Asm
mash dist -p 16 -d 0.2 SAM/Z1/mash/Asm.msh SAM/Z1/mash/Asm.msh|gzip > SAM/Z1/mash/Asm.dist.gz
```


### Align reads to ref:
```bash
bwa index $LFR/Source/REF/zymo/D6305.rRNA.fa
bwa mem -t 32 -db $LFR/Source/REF/zymo/D6305.rRNA.fa SAM/Z1/clean/fastp.sort.{1,2}.fq.gz \
| samtools sort -@ 3 - | samtools view -@ 3 -b - > SAM/Z1/READ/fastp2CloseRef.bam
samtools depth SAM/Z1/READ/fastp2CloseRef.bam -a > SAM/Z1/READ/fastp2CloseRef.bam.depth

for i in Z{2,3,8,9};do
	bwa mem -t 8 -db $LFR/Source/REF/zymo/D6305.rRNA.fa SAM/$i/clean/fastp.sort.{1,2}.fq.gz \
| samtools sort -T /dev/shm/samtools.$i -@ 3 - | samtools view -@ 3 -b - > SAM/$i/READ/fastp2CloseRef.bam && \
samtools depth SAM/$i/READ/fastp2CloseRef.bam -a > SAM/$i/READ/fastp2CloseRef.bam.depth &
done


# shotgun zymo
/hwfssz1/ST_META/F16ZQSB1SY3016_InfantMeta_LJH/User/yangfangming/tmp/D6300.Q.raw.path
bwa mem -t 16 -db $LFR/Source/REF/zymo/D6305.rRNA.fa /zfssz6/CNGB_DATA/BGISEQ01/BGISEQ500/P17H10200N0283_Temp/CL200105917_L02/CL200105917_L02_563_{1,2}.fq.gz | samtools sort -@ 3 - | samtools view -@ 3 -b - \
> TEST/shotgun/read2CloseRef.bam && samtools depth TEST/shotgun/read2CloseRef.bam -a > TEST/shotgun/read2CloseRef.bam.depth
```


# Rescaffolding test
```bash
#get long
awk -F "_" '(NF==9){if($9>999){p=1}else{p=0}} (p==1){print}' SAM/Z1/summary.BI.megahit.clip.fasta \
> SAM/Z1/CLIP/gt1kb.fasta
awk -F "_" '(NF==9){if($9<1000){p=1}else{p=0}} (p==1){print}' SAM/Z1/summary.BI.megahit.clip.fasta \
> SAM/Z1/CLIP/le1kb.fasta
bwa index SAM/Z1/CLIP/gt1kb.fasta
bwa bwasw -t 16 SAM/Z1/CLIP/gt1kb.fasta SAM/Z1/CLIP/le1kb.fasta -f SAM/Z1/CLIP/bwasw.1kb.sam

#
awk '($7=="Escherichia_coli"&&$8>99){print}' SAM/Z1/CLIP/clip2zymo.bead.anno > SAM/Z1/CLIP/Ecoli.bead.anno
for i in `sort -k2,2nr SAM/Z1/CLIP/Ecoli.bead.anno|cut -f1|head -30`;do grep -A 1 $i SAM/Z1/summary.BI.megahit.clip.fasta;
done > SAM/Z1/CLIP/Ecoli.clip.fasta &
for i in `sort -k2,2nr SAM/Z1/CLIP/Ecoli.bead.anno|cut -f1`;do cat SAM/Z1/Assemble_BI/$i/sort.1.fq;
done > SAM/Z1/CLIP/Ecoli.read.1.fasta &
for i in `sort -k2,2nr SAM/Z1/CLIP/Ecoli.bead.anno|cut -f1`;do cat SAM/Z1/Assemble_BI/$i/sort.2.fq;
done > SAM/Z1/CLIP/Ecoli.read.2.fasta

bwa index SAM/Z1/CLIP/Ecoli.clip.fasta
#bwa bwasw -M -t 16 SAM/Z1/CLIP/Ecoli.fasta SAM/Z1/CLIP/Ecoli.fasta -f SAM/Z1/CLIP/Ecoli.bwasw.sam
bwa mem -a -t16  SAM/Z1/CLIP/Ecoli.clip.fasta SAM/Z1/CLIP/Ecoli.clip.fasta > SAM/Z1/CLIP/Ecoli.clip.bwamem.sam
samtools view SAM/Z1/CLIP/Ecoli.bwamem.sam|awk '$1!=$3'|column -t|les
abyss-longseqdist -k14  SAM/Z1/CLIP/Ecoli.clip.bwamem.sam | grep -v "l=" > SAM/Z1/CLIP/Ecoli.bwamem.dist.dot
#reads to clip
bwa mem -a -t16 -S -P -k40 SAM/Z1/CLIP/Ecoli.fasta SAM/Z1/CLIP/Ecoli.read.fasta > SAM/Z1/CLIP/Ecoli.read.bwamem.sam
samtools view SAM/Z1/CLIP/Ecoli.read.bwamem.sam|les

#blast
makeblastdb -dbtype nucl -in SAM/Z1/CLIP/Ecoli.clip.fasta
blastn -num_threads 8 -perc_identity 99 -db SAM/Z1/CLIP/Ecoli.clip.fasta -query SAM/Z1/CLIP/Ecoli.clip.fasta -outfmt '6 std qlen slen' -out SAM/Z1/CLIP/Ecoli.clip.blastn.m6

samtools view SAM/Z1/CLIP/Ecoli.blastn.sam|awk '$1!=$3'|column -t|les
quast.py -t 16 SAM/Z1/CLIP/Ecoli.clip.fasta -R $LFR/Source/REF/zymo/D6305.fix.rRNA.fa -o SAM/Z1/CLIP/Ecoli.clip_quast
AdjList -k45 -m 100 --dot  SAM/Z1/CLIP/Ecoli.clip.fasta > SAM/Z1/CLIP/Ecoli.clip.dot

```

#AbySS test
```bash
mkdir SAM/Z1/ABYSS/
for i in `sort -k2,2nr SAM/Z1/CLIP/Ecoli.bead.anno|cut -f1|head -10`;
do cat SAM/Z1/Assemble_BI/$i/sort.1.fq|awk 'FNR%4==1{gsub("/1$","");gsub("/","/1\tBX:Z:")}{print}';
done > SAM/Z1/ABYSS/t.sort.1.fq &
for i in `sort -k2,2nr SAM/Z1/CLIP/Ecoli.bead.anno|cut -f1|head -10`;
do cat SAM/Z1/Assemble_BI/$i/sort.2.fq|awk 'FNR%4==1{gsub("/2$","");gsub("/","/2\tBX:Z:")}{print}';
done > SAM/Z1/ABYSS/t.sort.2.fq


abyss-pe k=48 name=SAM/Z1/ABYSS/Z1 \
	pe="pea" pea='SAM/Z1/ABYSS/t.sort.1.fq SAM/Z1/ABYSS/t.sort.2.fq' \
	mp='mpa' mpa='SAM/Z1/ABYSS/t.sort.1.fq SAM/Z1/ABYSS/t.sort.2.fq' \
	lr='lra' lra='SAM/Z1/ABYSS/t.sort.1.fq SAM/Z1/ABYSS/t.sort.2.fq'

cd SAM/Z1/ABYSS/

abyss-pe k=48 name=Z1 in='t.sort.1.fq t.sort.2.fq'
quast.py -t 16 Z1-scaffolds.fa -R $LFR/Source/REF/zymo/D6305.rRNA.fa -o quast -m 100 -1 t.sort.1.fq -2 t.sort.2.fq


abyss-pe k=48 name=Z1 \
 lib='pea' pea='t.sort.1.fq t.sort.2.fq' \
 lra='lra' lra='t.sort.1.fq t.sort.2.fq' \
 long='longa' longa='../summary.BI.megahit.clip.1kb.fasta' \
 c=5 j=16

megahit -t 16 --k-step 6 --prune-level 0 --min-count 1 -1 t.sort.1.fq -2 t.sort.2.fq -o megahit |tee megahit.log
quast.py -t 16 megahit/final.contigs.fa -R $LFR/Source/REF/zymo/D6305.rRNA.fa -o megahit/quast -m 100 -1 t.sort.1.fq -2 t.sort.2.fq
```

# adjget test
```bash
# STEP1  cluster
vsearch --threads 16 --cluster_fast  SAM/Z1/summary.BI.megahit.clip.fasta --strand both --fasta_width 0 --centroids SAM/Z1/CLIP/all.clip.clust.fasta --uc SAM/Z1/CLIP/all.clip.uc --id .99 -iddef 0 &> SAM/Z1/CLIP/all.clip.clust.log

#rename clip.fasta
perl filter.pl < SAM/Z1/summary.BI.megahit.clip.fasta > SAM/Z1/CLIP/all.clip.rename.fasta
#cluster
vsearch --threads 16 --cluster_fast  SAM/Z1/CLIP/all.clip.rename.fasta --strand both --fasta_width 0 --centroids SAM/Z1/CLIP/all.clip.rename.clust.fa \
--uc SAM/Z1/CLIP/all.clip.rename.clust.uc --id .99 -iddef 0 &> SAM/Z1/CLIP/all.clip.rename.clust.log

#makedb
makeblastdb -dbtype nucl -in SAM/Z1/CLIP/all.clip.rename.clust.fa
#blast
blastn -num_threads 12 -perc_identity 99 -db SAM/Z1/CLIP/all.clip.rename.clust.fa -query SAM/Z1/CLIP/all.clip.rename.clust.fa -outfmt '6 std qlen slen' -out SAM/Z1/CLIP/all.clip.rename.clust.blastn.m6
#scaffold
#sed 's/_/./g' SAM/Z1/CLIP/all.clip.rename.clust.uc > SAM/Z1/CLIP/all.clip.rename.clust.uc2
metabbq AdjGet -i SAM/Z1/CLIP/all.clip.rename.clust.fa -u SAM/Z1/CLIP/all.clip.rename.clust.uc2 -b SAM/Z1/CLIP/all.clip.rename.clust.blastn.m6 -m 5 -o SAM/Z1/CLIP/all.clip.rename.clust -v -n 5 &> SAM/Z1/CLIP/all.clip.rename.clust.adj.log

awk 'FNR%2==1{if($2>300){p=1}else{p=0}}(p==1){print}' SAM/Z1/CLIP/all.clip.rename.clust.merge.fa > SAM/Z1/CLIP/all.clip.rename.clust.merge.3k.fa
blastn -num_threads 8 -perc_identity 99 -db $LFR/Source/REF/zymo/D6305.fix.rRNA.fa -query SAM/Z1/CLIP/all.clip.rename.clust.merge.3k.fa -outfmt '6 std qlen slen qcovs' -out SAM/Z1/CLIP/all.clip.rename.clust.merge.3k.ref.m6


sed 's/ /_/g' SAM/Z1/CLIP/all.clip.rename.clust.merge.fa > SAM/Z1/CLIP/all.clip.rename.clust.merge.tmp.fa

quast.py -m 3000 -t 16 SAM/Z1/CLIP/all.clip.rename.clust.merge.tmp.fa -R $LFR/Source/REF/zymo/D6305.fix.rRNA.fa -o SAM/Z1/CLIP/all.clip.rename.clust.merge_quast
#

#
vsearch --threads 16 --cluster_fast SAM/Z1/CLIP/all.clip.rename.clust.merge.tmp.fa --strand both --fasta_width 0 --minsl .9 \
--centroids SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.fa --uc SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.uc --id .99 -iddef 0

#filter
perl -e 'open I, "<SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.uc";while(<I>){ next unless $_ =~ /^C/;
  @s=split;$HS{$s[8]} = $s[2] if $s[2] > 8};close I;$/=">";while(<>){
    chomp;next unless $_;@s=split "\n";if(exists $HS{$s[0]}){print ">$HS{$s[0]}\_$_"}}' < SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.fa \
    > SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.fa


quast.py -m 4000 -t 16 SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.fa -R $LFR/Source/REF/zymo/D6305.fix.rRNA.fa -o SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5_quast

#rRNA predict
barrnap --kingdom bac --threads 4 --reject 0.1 SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.fa \
--outseq SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.fa &> SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.log
grep -A1 "16S" SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap > SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.fa
vsearch --threads 16 --cluster_fast SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.fa --strand both --fasta_width 0 --sizeout \
--centroids SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.fa --uc SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.uc --id .98
vsearch --sortbysize SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.fa --fasta_width 0 --minsize 2 --sizeout --output SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.size2.fa
perl -e 'while(<>){if($_=~/^>16S_rRNA\:\:(\d+\_\d+\.\d+\_\d+\_\d+)_BI.*size=(\d+)$/){print ">$1_$2\n"}else{print}}' SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.size2.fa \
> SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.size2.rename.fa
cat $LFR/Source/REF/zymo/D6305.ssrRNA.fasta SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.size2.rename.fa > SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.size2.ref.fa
mafft --thread 16 SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.size2.ref.fa > SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.size2.ref.mafft.fa
sed 's/^.*::/>/;s/:/./g;s/;size=/_/' SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.size2.ref.mafft.fa > SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.size2.mafft.tmp.fa
trimal -gappyout -in SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.size2.mafft.tmp.fa -out SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.size2.mafft.trimal.fa


```

# build tree
```bash
#grep --no-group-separator -A1 "S_1" $LFR/Source/REF/zymo/D6305.ssrRNA.fasta > TEST/D6305.ssrRNA1.fasta
cat $LFR/Source/REF/zymo/D6305.ssrRNA.fasta SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.fa > SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.ref.fa
mafft --thread 16 SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.ref.fa > SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.ref.mafft.fa
sed 's/^.*::/>/;s/_BI.*$//' SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.ref.mafft.fa > SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.ref.mafft.tmp.fa
trimal -gappyout -in SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.ref.mafft.tmp.fa -out SAM/Z1/CLIP/all.clip.rename.clust.merge.clust.m5.barrnap.16S.clust.ref.trimal.fa

#16S
blastn -num_threads 8 -perc_identity 99 -db SAM/Z1/CLIP/all.16S.clust.ref.fa -query SAM/Z1/CLIP/all.16S.fasta -outfmt '6 std qlen' -out SAM/Z1/CLIP/all.16S.ref.m6
#use clips 16S
grep -A1 "16S_rRNA" SAM/Z1/summary.BI.megahit.rRNA.fasta > SAM/Z1/CLIP/all.16S.fasta
vsearch --threads 16 --fastx_filter SAM/Z1/CLIP/all.16S.fasta --fastq_minlen 1000 --fastaout SAM/Z1/CLIP/all.16S_1k.fasta
vsearch --threads 16 --cluster_fast SAM/Z1/CLIP/all.16S_1k.fasta --strand both --fasta_width 0 --sizeout --centroids SAM/Z1/CLIP/all.16S.clust.fa --uc SAM/Z1/CLIP/all.16S.clust.uc --id .99
vsearch --sortbysize SAM/Z1/CLIP/all.16S.clust.fa --minsize 50 --sizeout --output SAM/Z1/CLIP/all.16S.clust.size2.fa
perl -e 'while(<>){if($_=~/^>16S_rRNA\:\:(BI\d+)\_(.*)\_flag=(\d+)\_multi=(\d+).*\_len=(\d+).*size=(\d+)$/){print ">$1_$2_f$3_$4_$5\[$6\]\n"}else{print}}' \
SAM/Z1/CLIP/all.16S.clust.size2.fa > SAM/Z1/CLIP/all.16S.clust.size2.rename.fa
#awk 'FNR%2==1{match($0,/size=([0-9]+)/,a);if(a[1]>1){p=1}else{p=0}}(p==1){print}' SAM/Z1/CLIP/all.16S.clust.fa > SAM/Z1/CLIP/all.16S.clust.size2.fa
grep -A1 --no-group-separator "S_1" $LFR/Source/REF/zymo/D6305.ssrRNA.fasta > SAM/Z1/CLIP/all.16S.clust.ref.fa
  makeblastdb -dbtype nucl -in SAM/Z1/CLIP/all.16S.clust.ref.fa
  blastn -num_threads 8 -perc_identity 99 -db SAM/Z1/CLIP/all.16S.clust.ref.fa -query SAM/Z1/CLIP/all.16S.clust.size2.rename.fa -outfmt '6 std qlen' -out SAM/Z1/CLIP/all.16S.clust.size2.rename.m6
cat SAM/Z1/CLIP/all.16S.clust.size2.rename.fa >> SAM/Z1/CLIP/all.16S.clust.ref.fa
mafft --thread 16 SAM/Z1/CLIP/all.16S.clust.ref.fa > SAM/Z1/CLIP/all.16S.clust.ref.mafft.fa
sed 's/^.*::/>/;s/:/./g;s/;size=/_/' SAM/Z1/CLIP/all.16S.clust.ref.mafft.fa > SAM/Z1/CLIP/all.16S.clust.ref.mafft.tmp.fa
trimal -gappyout -in SAM/Z1/CLIP/all.16S.clust.ref.mafft.tmp.fa -out SAM/Z1/CLIP/all.16S.clust.ref.mafft.trimal.fa

#23S
blastn -num_threads 8 -perc_identity 98 -db SAM/Z1/CLIP/all.23S.clust.ref.fa -query SAM/Z1/CLIP/all.23S.fasta -outfmt '6 std qlen' -out SAM/Z1/CLIP/all.23S.ref.m6
#use clips 23S

grep -A1 "23S_rRNA" SAM/Z1/summary.BI.megahit.rRNA.fasta > SAM/Z1/CLIP/all.23S.fasta
vsearch --threads 16 --fastx_filter SAM/Z1/CLIP/all.23S.fasta --fastq_minlen 1000 --fastaout SAM/Z1/CLIP/all.23S_1k.fasta
vsearch --threads 16 --cluster_fast SAM/Z1/CLIP/all.23S_1k.fasta --strand both --fasta_width 0 --sizeout --centroids SAM/Z1/CLIP/all.23S.clust.fa --uc SAM/Z1/CLIP/all.23S.clust.uc --id .99
vsearch --sortbysize SAM/Z1/CLIP/all.23S.clust.fa --minsize 10 --sizeout --output SAM/Z1/CLIP/all.23S.clust.size2.fa
perl -e 'while(<>){if($_=~/^>23S_rRNA\:\:(BI\d+)\_(.*)\_flag=(\d+)\_multi=(\d+).*\_len=(\d+).*size=(\d+)$/){print ">$1_$2_f$3_$4_$5\[$6\]\n"}else{print}}' \
SAM/Z1/CLIP/all.23S.clust.size2.fa > SAM/Z1/CLIP/all.23S.clust.size2.rename.fa
#awk 'FNR%2==1{match($0,/size=([0-9]+)/,a);if(a[1]>1){p=1}else{p=0}}(p==1){print}' SAM/Z1/CLIP/all.23S.clust.fa > SAM/Z1/CLIP/all.23S.clust.size2.fa
awk '/>/{if($0~/23S/){p=1}else{p=0}}(p==1){print}' $LFR/Source/REF/zymo/D6305.markers.fa > SAM/Z1/CLIP/all.23S.clust.ref.fa
  makeblastdb -dbtype nucl -in SAM/Z1/CLIP/all.23S.clust.ref.fa
  blastn -num_threads 8 -perc_identity 99 -db SAM/Z1/CLIP/all.23S.clust.ref.fa -query SAM/Z1/CLIP/all.23S.clust.size2.rename.fa -outfmt '6 std qlen' -out SAM/Z1/CLIP/all.23S.clust.size2.rename.m6
cat SAM/Z1/CLIP/all.23S.clust.size2.rename.fa >> SAM/Z1/CLIP/all.23S.clust.ref.fa
mafft --thread 16 SAM/Z1/CLIP/all.23S.clust.ref.fa > SAM/Z1/CLIP/all.23S.clust.ref.mafft.fa
sed 's/^.*::/>/;s/:/./g;s/;size=/_/' SAM/Z1/CLIP/all.23S.clust.ref.mafft.fa > SAM/Z1/CLIP/all.23S.clust.ref.mafft.tmp.fa
trimal -gappyout -in SAM/Z1/CLIP/all.23S.clust.ref.mafft.tmp.fa -out SAM/Z1/CLIP/all.23S.clust.ref.mafft.trimal.fa


#overlap of 16S and 23S
grep -Eo "BI[0-9]+" SAM/Z1/CLIP/all.23S.fasta|sort|uniq > SAM/Z1/CLIP/all.23S.BID
grep -Eo "BI[0-9]+" SAM/Z1/CLIP/all.16S.fasta|sort|uniq > SAM/Z1/CLIP/all.16S.BID
```



#test LOTU phylotree
```bash
barrnap --threads 4 --kingdom bac < SAM/Z2/CLIP/all.clust.fa > SAM/Z2/CLIP/all.clust.fa.barrnap.gff 2> SAM/Z2/CLIP/all.clust.fa.barrnap.log
metabbq getAmpSeq -r B -g SAM/Z2/CLIP/all.clust.fa.barrnap.gff -i SAM/Z2/CLIP/all.clust.fa -o SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.fa
awk -F "=|:" 'FNR%2==1{if($2>=50){p=1}else{p=0}}(p==1){print}' SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.fa > SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.size50.fa
cat SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.size50.fa $LFR/Source/REF/zymo/D6305.rRNA.fa > SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.size50.addRef.fa
mafft --thread 16 SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.size50.addRef.fa > SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.size50.addRef.mafft.fa

trimal -gappyout -in SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.size50.addRef.mafft.fa -out SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.size50.addRef.mafft.trimal.fa

##
cut -f2 SAM/Z2/CLIP/pick.LID.txt|xargs -n1 grep -h -A1 - SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.fa \
> SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.cov3.fa 2>/dev/null


#Only 16S
metabbq getAmpSeq -r S -p 100 -g SAM/Z2/CLIP/all.clust.fa.barrnap.gff -i SAM/Z2/CLIP/all.clust.fa -o SAM/Z2/CLIP/all.clust.fa.barrnap.SSU.fa
perl -e 'open I,"<SAM/Z2/CLIP/pick.LID.txt";while(<I>){@s=split;$HS{$s[1]}=1};$/=">";while(<>){chomp;next unless $_;@s=split /;|\n/;if(exists $HS{$s[0]}){print ">".$_}}' < SAM/Z2/CLIP/all.clust.fa.barrnap.SSU.fa > SAM/Z2/CLIP/all.clust.fa.barrnap.SSU.cov3.fa
cat SAM/Z2/CLIP/all.clust.fa.barrnap.SSU.cov3.fa $LFR/Source/REF/zymo/D6305.ssrRNA.fasta > SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.cov3.addRef.fa
mafft --thread 16 SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.cov3.addRef.fa > SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.cov3.addRef.mafft.fa

trimal -gappyout -in SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.cov3.addRef.mafft.fa -out SAM/Z2/CLIP/all.clust.fa.barrnap.rRNA.cov3.addRef.mafft.trimal.fa


```

#test LOTU_270:
```bash
metabbq IO tree2fa -i SAM/SBR/PREDICT/barrnap.ssu.fa -r phylum -d CLADE0.72 -t SAM/SBR/KRAKEN/test/data/added.txt -o SAM/SBR/PREDICT/phylum.ssu.fa -v
mafft --thread 16 SAM/SBR/PREDICT/phylum.ssu.fa > SAM/SBR/PREDICT/phylum.ssu.mafft.fa
trimal -gappyout -in SAM/SBR/PREDICT/phylum.ssu.mafft.fa -out SAM/SBR/PREDICT/phylum.ssu.mafft.trimal.fa

raxmlHPC-HYBRID-AVX -f a -p 12345 -s SAM/SBR/PREDICT/phylum.ssu.mafft.trimal.fa -x 12345 -T 20 -# 20 -m GTRCAT  -w $PWD/SAM/SBR/PREDICT -n phylum.ssu.raxml

```

# test raxml
```bash
metabbq IO clade2tree -i SAM/SBR/CLIP/id90def4.clust.fa -a SAM/SBR/ANNO/CLIP.map.merge.clip.anno \
 -d $LFR/Source/REF/KRAKEN2/rDNA -s SAM/SBR/CLIP/id90def4.clade.uc -o SAM/SBR/KRAKEN/test -m bac -v

awk '/^>/{if($0~/23S|28S/){p=1}else{p=0}};(p==1){print}' $LFR/Source/REF/zymo/D6305.markers.fa > $LFR/Source/REF/zymo/D6305.markers.LSU.fa

awk '/^>/{if($0~/16S|18S|23S|25S|28S/){p=1}else{p=0}};(p==1){print}' $LFR/Source/REF/zymo/D6305.markers.fa > $LFR/Source/REF/zymo/D6305.markers.both.fa

#perl -e 'while(<>){@s=split;if($s[3]=~/16S|18S/ && $s[5] eq "+"){print "$s[0]:$s[6]-$s[7]\n"}}' < $LFR/Source/REF/zymo/D6305.rRNA.barrnap.bed > $LFR/Source/REF/zymo/D6305.rRNA.barrnap.16S+.region
#perl -e 'while(<>){@s=split;if($s[3]=~/16S|18S/ && $s[5] eq "-"){print "$s[0]:$s[6]-$s[7]\n"}}' < $LFR/Source/REF/zymo/D6305.rRNA.barrnap.bed > $LFR/Source/REF/zymo/D6305.rRNA.barrnap.16S-.region
samtools faidx -r $LFR/Source/REF/zymo/D6305.rRNA.barrnap.16S+.region $LFR/Source/REF/zymo/D6305.rRNA.fa > $LFR/Source/REF/zymo/D6305.select.both.ssu.rep.fa
samtools faidx -i -r $LFR/Source/REF/zymo/D6305.rRNA.barrnap.16S-.region $LFR/Source/REF/zymo/D6305.rRNA.fa >> $LFR/Source/REF/zymo/D6305.select.both.ssu.rep.fa
sed -i 's/:/ /' $LFR/Source/REF/zymo/D6305.select.both.ssu.rep.fa

#perl -e 'while(<>){@s=split;if($s[3]=~/23S|28S|25S/ && $s[5] eq "+"){print "$s[0]:$s[6]-$s[7]\n"}}' < $LFR/Source/REF/zymo/D6305.rRNA.barrnap.bed > $LFR/Source/REF/zymo/D6305.rRNA.barrnap.23S+.region
#perl -e 'while(<>){@s=split;if($s[3]=~/23S|28S|25S/ && $s[5] eq "-"){print "$s[0]:$s[6]-$s[7]\n"}}' < $LFR/Source/REF/zymo/D6305.rRNA.barrnap.bed > $LFR/Source/REF/zymo/D6305.rRNA.barrnap.23S-.region
samtools faidx -r $LFR/Source/REF/zymo/D6305.rRNA.barrnap.23S+.region $LFR/Source/REF/zymo/D6305.rRNA.fa > $LFR/Source/REF/zymo/D6305.select.both.lsu.rep.fa
samtools faidx -i -r $LFR/Source/REF/zymo/D6305.rRNA.barrnap.23S-.region $LFR/Source/REF/zymo/D6305.rRNA.fa >> $LFR/Source/REF/zymo/D6305.select.both.lsu.rep.fa
sed -i 's/:/ /' $LFR/Source/REF/zymo/D6305.select.both.lsu.rep.fa


```

For species:

**with anno*
```bash
metabbq IO tree2fa -s SAM/SBR/CLIP/all.SSU.fa -l SAM/SBR/CLIP/all.LSU.fa -r species -d CLADE0.99 -t SAM/SBR/KRAKEN/test/data/added.txt -o SAM/SBR/PREDICT/sp.a.fa -a -v
#SSU:497,LSU:131,BOTH
cat $LFR/Source/REF/zymo/D6305.ssrRNA.rep.fasta SAM/SBR/PREDICT/sp.a.SSU.fa > SAM/SBR/PREDICT/sp.a.ssu.fa

mafft --thread 16 SAM/SBR/PREDICT/sp.a.ssu.fa > SAM/SBR/PREDICT/sp.a.ssu.mafft.fa
trimal -gappyout -in SAM/SBR/PREDICT/sp.a.ssu.mafft.fa -out SAM/SBR/PREDICT/sp.a.ssu.mafft.trimal.fa

raxmlHPC-HYBRID-AVX -f a -x 12345 -p 12345 -s SAM/SBR/PREDICT/sp.a.ssu.mafft.trimal.fa -T 20 -# 40 -m GTRCAT  -w $PWD/SAM/SBR/PREDICT -n sp.a.ssu.raxml


```

```bash
metabbq IO tree2fa -s SAM/SBR/CLIP/all.SSU.fa -l SAM/SBR/CLIP/all.LSU.fa -r species -d CLADE0.99 -t SAM/SBR/KRAKEN/test/data/added.txt -o SAM/SBR/PREDICT/species.fa -p -v
```
Written 32307 clade seqs. 37697 records cannot found seqs.
SSU: 26572 (561+25475+536) | LSU: 7240(132+6158+950) | BOTH:1505(72+1433+1486)

```bash
#SSU:497,LSU:131,BOTH
perl -e 'open I,"<SAM/SBR/PREDICT/species.SSU.fa";$o=$/;$/=">";while(<I>){chomp;@s=split /\n/;$id=shift @s;$HS{$id}=join("",@s)};$/=$o;
open M,"<SAM/SBR/summary.BI.megahit.clip.metadata.tsv";while(<M>){@s=split;my %H=();while($s[7]=~s/(FWD|REV)(\+|\-)//){
  $H{$1}++;$H{$2}++};if($H{"FWD"} ==1 && $H{"REV"} == 0){$s[6]=~s/C//;$AD{"$s[0].$s[2].$s[6]"}++ }};
while(<>){@s=split /\t/;if($s[3]>1399&&$AD{$s[2]}){$count ++;if($HS{$s[0]}){print ">$s[0]\n$HS{$s[0]}\n"}else{print STDERR "$s[0]\n"}}};print STDERR $count."\n"' \
< SAM/SBR/PREDICT/species.SSU.inf > SAM/SBR/PREDICT/species.ssu.filter.fa 2> SAM/SBR/PREDICT/species.ssu.filter.log

cat $LFR/Source/REF/zymo/D6305.select.both.ssu.rep.fa SAM/SBR/PREDICT/species.ssu.filter.fa > SAM/SBR/PREDICT/species+ref.ssu.fa

mafft --thread 48 SAM/SBR/PREDICT/species+ref.ssu.fa > SAM/SBR/PREDICT/species+ref.ssu.mafft.fa
trimal -gappyout -in SAM/SBR/PREDICT/species+ref.ssu.mafft.fa -out SAM/SBR/PREDICT/species+ref.ssu.mafft.trimal.fa

raxmlHPC-HYBRID-AVX -f a -p 12345 -x 12345 -s SAM/SBR/PREDICT/species+ref.ssu.mafft.trimal.fa -T 40 -# 40 -m GTRCAT  -w $PWD/SAM/SBR/PREDICT -n species+ref.ssu.raxml

perl phylo.trans.pl < SAM/SBR/PREDICT/RAxML_bipartitionsBranchLabels.species+ref.ssu.raxml > SAM/SBR/PREDICT/RAxML_bipartitionsBranchLabels.species.a.ssu.raxml

#both
perl -e 'open I,"<SAM/SBR/PREDICT/species.both.SSU.fa";$o=$/;$/=">";while(<I>){chomp;@s=split /\n/;$id=shift @s;$HS{$id}=join("",@s)};$/=$o;
while(<>){@s=split /\t/;if($s[3]>999&&$s[4]>999){$count ++;if($HS{$s[0]}){print ">$s[0]\n$HS{$s[0]}\n"}else{print STDERR "$s[0]\n"}}};print STDERR $count."\n"' \
< SAM/SBR/PREDICT/species.both.inf > SAM/SBR/PREDICT/species.both.ssu.filter.fa 2> SAM/SBR/PREDICT/species.both.ssu.filter.log
perl -e 'open I,"<SAM/SBR/PREDICT/species.both.LSU.fa";$o=$/;$/=">";while(<I>){chomp;@s=split /\n/;$id=shift @s;$HS{$id}=join("",@s)};$/=$o;
while(<>){@s=split /\t/;if($s[3]>999&&$s[4]>999){$count ++;if($HS{$s[0]}){print ">$s[0]\n$HS{$s[0]}\n"}else{print STDERR "$s[0]\n"}}};print STDERR $count."\n"' \
< SAM/SBR/PREDICT/species.both.inf > SAM/SBR/PREDICT/species.both.lsu.filter.fa 2> SAM/SBR/PREDICT/species.both.lsu.filter.log

cat $LFR/Source/REF/zymo/D6305.select.both.ssu.rep.fa SAM/SBR/PREDICT/species.both.ssu.filter.fa > SAM/SBR/PREDICT/species.both.ssu+.filter.fa
cat $LFR/Source/REF/zymo/D6305.select.both.lsu.rep.fa SAM/SBR/PREDICT/species.both.lsu.filter.fa > SAM/SBR/PREDICT/species.both.lsu+.filter.fa

#ssu part
mafft --thread 16 SAM/SBR/PREDICT/species.both.ssu+.filter.fa > SAM/SBR/PREDICT/species.both.ssu.mafft.fa
trimal -gappyout -in SAM/SBR/PREDICT/species.both.ssu.mafft.fa -out SAM/SBR/PREDICT/species.both.ssu.trimal.fa
#lsu part
mafft --thread 16 SAM/SBR/PREDICT/species.both.lsu+.filter.fa > SAM/SBR/PREDICT/species.both.lsu.mafft.fa
trimal -gappyout -in SAM/SBR/PREDICT/species.both.lsu.mafft.fa -out SAM/SBR/PREDICT/species.both.lsu.trimal.fa
# merge
perl -e 'open I,"<SAM/SBR/PREDICT/species.both.ssu.trimal.fa";$/=">";while(<I>){chomp;next unless $_; @s=split;$id=shift @s;$id=~s/ .*//;$HS{$id}=join("",@s)};
while(<>){chomp;next unless $_; @s=split;$id=shift @s;$id=~s/ .*//;die "ssu missin for $id" unless exists $HS{$id};
$str = join("",$HS{$id},@s); @parts = $str =~ /(.{1,60})/g;print ">$id\n".join("\n",@parts)."\n"}' \
< SAM/SBR/PREDICT/species.both.lsu.trimal.fa > SAM/SBR/PREDICT/species.both.trimal.fa


raxmlHPC-HYBRID-AVX -f a -p 12345 -s SAM/SBR/PREDICT/species.both.trimal.fa -x 12345 -T 40 -# 40 -m GTRCAT  -w $PWD/SAM/SBR/PREDICT/RAXML -n species+ref.both.raxml
perl phylo.trans.pl SAM/SBR/PREDICT/species.both.inf < SAM/SBR/PREDICT/RAXML/RAxML_bipartitionsBranchLabels.species+ref.both.raxml > SAM/SBR/PREDICT/RAXML/RAxML_bipartitionsBranchLabels.species+ref.both.anno.raxml

```
> Written 49750 clade seqs. 20253 records cannot found seqs.
> SSU: 36902 (558+33561+2783) | LSU: 17640(128+15589+1923) | BOTH:4792(126+4666+4706)


For genus:
```bash
metabbq IO tree2fa -s SAM/SBR/CLIP/all.SSU.fa -l SAM/SBR/CLIP/all.LSU.fa -r genus -d CLADE0.97 -t SAM/SBR/KRAKEN/test/data/added.txt -o SAM/SBR/PREDICT/genus.fa -v
cat $LFR/Source/REF/zymo/D6305.ssrRNA.fasta SAM/SBR/PREDICT/genus.ssu.fa > SAM/SBR/PREDICT/genus+ref.ssu.fa
mafft --thread 16 SAM/SBR/PREDICT/genus+ref.ssu.fa > SAM/SBR/PREDICT/genus+ref.ssu.mafft.fa
trimal -gappyout -in SAM/SBR/PREDICT/genus+ref.ssu.mafft.fa -out SAM/SBR/PREDICT/genus+ref.ssu.mafft.trimal.fa
raxmlHPC-HYBRID-AVX -f a -p 12345 -s SAM/SBR/PREDICT/genus+ref.ssu.mafft.trimal.fa -x 12345 -T 20 -# 20 -m GTRCAT  -w $PWD/SAM/SBR/PREDICT -n genus+ref.ssu.raxml

```
> Written 398 seqs (SSU:335,LSU:288,BOTH:225). 8 records cannot found seqs
> SSU: 5691 (278+4946+467) | LSU: 3002(195+2252+555) | BOTH:1114(225+889+1022)

For family:
```bash
metabbq IO tree2fa -s SAM/SBR/CLIP/all.SSU.fa -l SAM/SBR/CLIP/all.LSU.fa -r family -d CLADE0.89 -t SAM/SBR/KRAKEN/test/data/added.txt -o SAM/SBR/PREDICT/family.fa -v
cat $LFR/Source/REF/zymo/D6305.ssrRNA.fasta SAM/SBR/PREDICT/family.SSU.fa > SAM/SBR/PREDICT/family+ref.ssu.fa
mafft --thread 16 SAM/SBR/PREDICT/family+ref.ssu.fa > SAM/SBR/PREDICT/family+ref.ssu.mafft.fa
trimal -gappyout -in SAM/SBR/PREDICT/family+ref.ssu.mafft.fa -out SAM/SBR/PREDICT/family+ref.ssu.mafft.trimal.fa
raxmlHPC-HYBRID-AVX -f a -p 12345 -s SAM/SBR/PREDICT/family+ref.ssu.mafft.trimal.fa -x 12345 -T 40 -# 40 -m GTRCAT  -w $PWD/SAM/SBR/PREDICT -n family+ref.ssu.raxml

```
> Written 560 clade seqs. 2606 records cannot found seqs.
> SSU: 530 (93+389+48) | LSU: 223(79+29+115) | BOTH:193(99+94+163)

For phylum:
```bash
metabbq IO tree2fa -s SAM/SBR/CLIP/all.SSU.fa -l SAM/SBR/CLIP/all.LSU.fa -r phylum -d CLADE0.72 -t SAM/SBR/KRAKEN/test/data/added.txt -o SAM/SBR/PREDICT/phylum.fa -v
```
  Written 34 clade seqs. 1310 records cannot found seqs.
  SSU: 32 (19+6+7) | LSU: 28(11+2+15) | BOTH:26(24+2+22)

  Written 27 clade seqs. 1317 records cannot found seqs.
  SSU: 25 (21+1+3) | LSU: 15(7++8) | BOTH:13(13++11)

```bash
#SSU
perl -e 'open I,"<SAM/SBR/PREDICT/phylum.SSU.fa";$o=$/;$/=">";while(<I>){chomp;@s=split /\n/;$id=shift @s;$HS{$id}=join("",@s)};$/=$o;
while(<>){@s=split /\t/;if($s[3]>999){$count ++;if($HS{$s[0]}){print ">$s[0]\n$HS{$s[0]}\n"}else{print STDERR "$s[0]\n"}}};print STDERR $count."\n"' \
< SAM/SBR/PREDICT/phylum.SSU.inf > SAM/SBR/PREDICT/phylum.ssu.filter.fa 2> SAM/SBR/PREDICT/species.ssu.filter.log

cat $LFR/Source/REF/zymo/D6305.select.both.ssu.rep.fa SAM/SBR/PREDICT/phylum.ssu.filter.fa > SAM/SBR/PREDICT/phylum+ref.ssu.fa
mafft --thread 16 SAM/SBR/PREDICT/phylum+ref.ssu.fa > SAM/SBR/PREDICT/phylum+ref.ssu.mafft.fa
trimal -gappyout -in SAM/SBR/PREDICT/phylum+ref.ssu.mafft.fa -out SAM/SBR/PREDICT/phylum+ref.ssu.mafft.trimal.fa
raxmlHPC-HYBRID-AVX -f a -p 12345 -s SAM/SBR/PREDICT/phylum+ref.ssu.mafft.trimal.fa -x 12345 -T 40 -# 40 -m GTRCAT  -w $PWD/SAM/SBR/PREDICT/RAXML -n phylum+ref.ssu.raxml
perl phylo.trans.pl SAM/SBR/PREDICT/phylum.SSU.inf < SAM/SBR/PREDICT/RAXML/RAxML_bipartitionsBranchLabels.phylum+ref.ssu.raxml > SAM/SBR/PREDICT/RAXML/RAxML_bipartitionsBranchLabels.phylum+ref.ssu.anno.raxml

#LSU
cat $LFR/Source/REF/zymo/D6305.select.both.lsu.rep.fa SAM/SBR/PREDICT/phylum.LSU.fa > SAM/SBR/PREDICT/phylum+ref.lsu.fa
mafft --thread 16 SAM/SBR/PREDICT/phylum+ref.lsu.fa > SAM/SBR/PREDICT/phylum+ref.lsu.mafft.fa
trimal -gappyout -in SAM/SBR/PREDICT/phylum+ref.lsu.mafft.fa -out SAM/SBR/PREDICT/phylum+ref.lsu.mafft.trimal.fa
raxmlHPC-HYBRID-AVX -f a -p 12345 -s SAM/SBR/PREDICT/phylum+ref.lsu.mafft.trimal.fa -x 12345 -T 40 -# 40 -m GTRCAT  -w $PWD/SAM/SBR/PREDICT/RAXML -n phylum+ref.lsu.raxml
perl phylo.trans.pl SAM/SBR/PREDICT/phylum.LSU.inf < SAM/SBR/PREDICT/RAXML/RAxML_bipartitionsBranchLabels.phylum+ref.lsu.raxml > SAM/SBR/PREDICT/RAXML/RAxML_bipartitionsBranchLabels.phylum+ref.lsu.anno.raxml


#both
perl -e 'open I,"<SAM/SBR/PREDICT/phylum.both.SSU.fa";$o=$/;$/=">";while(<I>){chomp;@s=split /\n/;$id=shift @s;$HS{$id}=join("",@s)};$/=$o;
while(<>){@s=split /\t/;if($s[3]>999&&$s[4]>999){$count ++;if($HS{$s[0]}){print ">$s[0]\n$HS{$s[0]}\n"}else{print STDERR "$s[0]\n"}}};print STDERR $count."\n"' \
< SAM/SBR/PREDICT/phylum.both.inf > SAM/SBR/PREDICT/phylum.both.ssu.filter.fa 2> SAM/SBR/PREDICT/species.both.ssu.filter.log
perl -e 'open I,"<SAM/SBR/PREDICT/phylum.both.LSU.fa";$o=$/;$/=">";while(<I>){chomp;@s=split /\n/;$id=shift @s;$HS{$id}=join("",@s)};$/=$o;
while(<>){@s=split /\t/;if($s[3]>999&&$s[4]>999){$count ++;if($HS{$s[0]}){print ">$s[0]\n$HS{$s[0]}\n"}else{print STDERR "$s[0]\n"}}};print STDERR $count."\n"' \
< SAM/SBR/PREDICT/phylum.both.inf > SAM/SBR/PREDICT/phylum.both.lsu.filter.fa 2> SAM/SBR/PREDICT/species.both.lsu.filter.log

cat $LFR/Source/REF/zymo/D6305.select.both.ssu.rep.fa SAM/SBR/PREDICT/phylum.both.ssu.filter.fa > SAM/SBR/PREDICT/phylum.both.ssu+.filter.fa
cat $LFR/Source/REF/zymo/D6305.select.both.lsu.rep.fa SAM/SBR/PREDICT/phylum.both.lsu.filter.fa > SAM/SBR/PREDICT/phylum.both.lsu+.filter.fa

#ssu part
mafft --thread 16 SAM/SBR/PREDICT/phylum.both.ssu+.filter.fa > SAM/SBR/PREDICT/phylum.both.ssu.mafft.fa
trimal -gappyout -in SAM/SBR/PREDICT/phylum.both.ssu.mafft.fa -out SAM/SBR/PREDICT/phylum.both.ssu.trimal.fa
#lsu part
mafft --thread 16 SAM/SBR/PREDICT/phylum.both.lsu+.filter.fa > SAM/SBR/PREDICT/phylum.both.lsu.mafft.fa
trimal -gappyout -in SAM/SBR/PREDICT/phylum.both.lsu.mafft.fa -out SAM/SBR/PREDICT/phylum.both.lsu.trimal.fa
# merge
perl -e 'open I,"<SAM/SBR/PREDICT/phylum.both.ssu.trimal.fa";$/=">";while(<I>){chomp;next unless $_; @s=split;$id=shift @s;$id=~s/ .*//;$HS{$id}=join("",@s)};
while(<>){chomp;next unless $_; @s=split;$id=shift @s;$id=~s/ .*//;die "ssu missin for $id" unless exists $HS{$id};
$str = join("",$HS{$id},@s); @parts = $str =~ /(.{1,60})/g;print ">$id\n".join("\n",@parts)."\n"}' \
< SAM/SBR/PREDICT/phylum.both.lsu.trimal.fa > SAM/SBR/PREDICT/phylum.both.trimal.fa


raxmlHPC-HYBRID-AVX -f a -p 12345 -s SAM/SBR/PREDICT/phylum.both.trimal.fa -x 12345 -T 40 -# 40 -m GTRCAT  -w $PWD/SAM/SBR/PREDICT/RAXML -n phylum+ref.both.raxml
perl phylo.trans.pl SAM/SBR/PREDICT/phylum.both.inf < SAM/SBR/PREDICT/RAXML/RAxML_bipartitionsBranchLabels.phylum+ref.both.raxml > SAM/SBR/PREDICT/RAXML/RAxML_bipartitionsBranchLabels.phylum+ref.both.anno.raxml

```


>("domain",0.66,"kingdom",0.7,"phylum",0.76,"class",0.77,"order",0.83,"family",0.89,"genus",0.97,"species",0.995);
