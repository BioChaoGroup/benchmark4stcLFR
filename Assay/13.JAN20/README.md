# Test data generated at Jan 2020

# Sample information

Above information were save in sheet [dataInfo.tsv](../../../Results/JAN20/dataInfo.tsv).

Experiments information under directory:
 `BGI Tech Solutions (Hongkong) Co., Ltd\stLFR4META - 文档\Experiments\`

# Analysis preparation

Organize directories. `ORI(GIN)` were organized by sequencing events; `SAM(PLE)` were organized by sample name.

```bash
mkdir -p ../../Results/JAN20/{ORI,SAM}
ln -s ../../Results/JAN20/ORI ./
```

Sample info are organized in excel file `stLFR4META - 文档\Experiments\2001范飞\20年测试样本信息对照表.xlsx`
Then generate directories for each sequence data:

```bash
awk '{system("mkdir -p ORI/"$4"/input");
system("ln -s "$5" ORI/"$4"/input/rawSeq_1.fq.gz");
sub("_1.fq","_2.fq",$5);
system("ln -s "$5" ORI/"$4"/input/rawSeq_2.fq.gz");
}' ../../Results/JAN20/dataInfo.tsv
```

Run QC step:

```bash
#for RCA bac
metabbq smk -j -np ORI/FF{1,2,3}_{1,3,4,5,7,8,17,22}/clean/BB.stat
#for RCA and non-RCA fungi
metabbq smk -j -np ORI/FF{4,5,6,A,B,C}_{9,10,11,12,13,20,15,16}/clean/BB.stat
#for Non-RCA bac & fungi
metabbq smk -j -np ORI/FF{7,8,9}_{1,3,17,4,5,22,7,8}/clean/BB.stat

#stat
for i in `ls -d ORI/FF*`;do metabbq stat basic -p $i -o $i/stat/basic.log & done
cat ORI/FF*/stat/basic.log|sort -nr|uniq > STAT/ORI.basic.log
```

After `fastp`, we can merge samples together:
**RCA bacteria**

```bash
ln -s ../../Results/JAN20/SAM ./
mkdir -p SAM/{U,S,O}{1,2,3}/{clean,.tmp.{1,2}}
for i in {1,2,3};do for j in {1,2};do
echo cat ORI/FF$i\_{1,17}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/U$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/U$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{3,4}/clean/fastp.sort.$j.fq  \| paste - - - - \| sort -T SAM/S$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/S$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{5,22}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/O$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/O$i/clean/fastp.sort.$j.fq \&
done;done
metabbq smk -j -np SAM/{U,S,O}{1,2,3}/clean/BB.stat
#
mkdir -p SAM/{M,Z}{1,2,3}
for i in {1,2,3};do
	ln -s ../../ORI/FF$i\_7/clean SAM/M$i/
	ln -s ../../ORI/FF$i\_8/clean SAM/Z$i/
done
```

**RCA fungi**

```bash
mkdir -p SAM/{S,O}{4,5,6}/{clean,.tmp.{1,2}}
for i in {4,5,6};do for j in {1,2};do
echo cat ORI/FF$i\_{9,10,11}/clean/fastp.sort.$j.fq  \| paste - - - - \| sort -T SAM/S$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/S$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{12,13,20}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/O$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/O$i/clean/fastp.sort.$j.fq \&
done;done
metabbq smk -j -np SAM/{S,O}{4,5,6}/clean/BB.stat
#
mkdir -p SAM/{M,Y}{4,5,6}
for i in {4,5,6};do
	ln -s ../../ORI/FF$i\_16/clean SAM/M$i/
	ln -s ../../ORI/FF$i\_15/clean SAM/Y$i/
done
```

**Non-RCA bacteria**

```bash
mkdir -p SAM/{S,O,Z}{7,8,9}/{clean,.tmp.{1,2}}
for i in {7,8,9};do for j in {1,2};do
echo cat ORI/FF$i\_{1,3,17}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/S$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/S$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{4,5,22}/clean/fastp.sort.$j.fq  \| paste - - - - \| sort -T SAM/O$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/O$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{7,8}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/Z$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/Z$i/clean/fastp.sort.$j.fq \&
done;done
metabbq smk -j -np SAM/{S,O,Z}{7,8,9}/clean/BB.stat
#
```

**Non-RCA fungi**

```bash
mkdir -p SAM/{S,O,M}{A,B,C}/{clean,.tmp.{1,2}}
for i in {A,B,C};do for j in {1,2};do
echo cat ORI/FF$i\_{9,10,11}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/S$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/S$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{12,13,20}/clean/fastp.sort.$j.fq  \| paste - - - - \| sort -T SAM/O$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/O$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{15,16}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/M$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/M$i/clean/fastp.sort.$j.fq \&
done;done
metabbq smk -j -np SAM/{S,O,M}{A,B,C}/clean/BB.stat
#
```

stat combined basic numbers:
```bash
for i in `ls -d SAM/??`;do mkdir -p $i/stat && metabbq stat basic -p $i -o $i/stat/basic.log & done
for i in `ls ORI|grep FF`;do awk -v i=$i '{print i"\t"$0}' ORI/$i/clean/BB.stat; done > STAT/ORI.BB.stat
for i in `ls SAM`;do awk -v i=$i '{print i"\t"$0}' SAM/$i/clean/BB.stat; done > STAT/SAM.BB.stat

cat SAM/{{M,Z}{1,2,3},{M,Y}{4,5,6}}/stat/basic.log|sort|uniq > STAT/SAM.basic.single.log
cat SAM/{{U,S,O}{1,2,3},{S,O}{4,5,6}}/stat/basic.log|sort|uniq > STAT/SAM.basic.RCA.log
cat SAM/{{S,O,Z}{7,8,9},{S,O,M}{A,B,C}}/stat/basic.log|sort|uniq > STAT/SAM.basic.NOR.log

```

Run following steps:

```bash
metabbq smk -j --resources mem=120 -np SAM/{{U,S,O,M,Z}{1,2,3},{S,O,M,Y}{4,5,6},{S,O,Z}{7,8,9},{S,O,M}{A,B,C}}/Assemble_BI/ID.lst
```

RCA isolate assemble:
```bash
# for bac
metabbq smk -j -np SAM/{U,S,O}{1,2,3}}/summary.BI.megahit.clip.all.fasta
# for fungi
metabbq smk --configfile config.F.yaml -j -np SAM/{S,O,M,Y}{4,5,6}/summary.BI.megahit.clip.all.fasta

```

# get fungal mock ref
```bash
#from sanger
cd $LFR/Source/REF/fungi/sanger
sh work.sh
vsearch --allpairs_global sanger.mock7.rRNA.fa --acceptall --blast6out sanger.mock7.rRNA.vsearch.allpairs.m6 --iddef 4
ITSx -t F --cpu 2 -i sanger.mock7.rRNA.fa -o sanger.mock7.ITSx &> sanger.mock7.ITSx.log
perl onetime.cutSubunits.pl
vsearch --allpairs_global sanger.mock7.ITSx.SSU.fa --acceptall --blast6out sanger.mock7.SSU.vsearch.allpairs.m6 --iddef 4
vsearch --allpairs_global sanger.mock7.ITSx.ITS.fa --acceptall --blast6out sanger.mock7.ITS.vsearch.allpairs.m6 --iddef 4
vsearch --allpairs_global sanger.mock7.ITSx.LSU.fa --acceptall --blast6out sanger.mock7.LSU.vsearch.allpairs.m6 --iddef 4
```
From above alignment we know the mock similarity could be as high as 99.5%! Thus we need a higher threashold for LOTU cluster.

Find following codes from below scripts:
**upstream analysis**
> Assess.JAN20.Rmd

**Mock samples assembly results**
> Assess.JAN20.MOCK.BAC.figures.Rmd
> Assess.JAN20.MOCK.FUNGI.figures.Rmd

**Real soil samples assembly results**
> Assess.JAN20.SOIL.BAC.figures.Rmd
> Assess.JAN20.SOIL.FUNGI.figures.Rmd



# Appends: Test codes
### test: range of rpb or kmers for assemble
```bash
mkln ../../Results/JAN20/TEST
mkdir TEST/Z1
ln -s ../../SAM/Z1/clean TEST/Z1/

metabbq smk --config p_cluster_minR=10 p_cluster_maxR=99 p_cluster_ranP=100 -j -npk TEST/Z1/mash/BI.1.fq
metabbq smk --config p_asm_minKmers=0 p_asm_minKmerCoverage=1 -j -npk TEST/Z1/batch.assemble.BI.sh
 metabbq smk --force --config p_asm_min=200 -j -pk TEST/Z1/summary.BI.megahit.clip.fasta

for i in Z{1,2,3};do
	mkdir -p SAM/$i/CLIP && vsearch --threads 8 --cluster_fast SAM/$i/summary.BI.megahit.clip.fasta --iddef $def --id 0.99 --strand both --fasta_width 0 --centroids SAM/$i/CLIP/preclust$def.fa -uc SAM/$i/CLIP/preclust$def.uc --sizeout &> SAM/$i/CLIP/preclust$def.log \
	&& awk -F "_" '/>/{if($9>1499&&$10>1){p=1}else{p=0}} (p==1){print}' SAM/$i/CLIP/preclust$def.fa > SAM/$i/CLIP/preclust$def\_1.5k.fasta \
	&& quast.py -t 8 SAM/$i/summary.BI.megahit.clip_1.5k.fasta -R $LFR/Source/REF/zymo/D6305.rRNA.fa -o SAM/$i/summary.BI.megahit.clip_1.5k_quast &
done

def=4
for i in M{4,5,6};do
	vsearch --threads 8 --cluster_fast SAM/$i/summary.BI.megahit.clip.fasta --iddef $def --id 0.99 --strand both --fasta_width 0 --centroids SAM/$i/CLIP/preclust$def.fa -uc SAM/$i/CLIP/preclust$def.uc --sizeout &> SAM/$i/CLIP/preclust$def.log &&
	awk -F "_|;size=" '/>/{if($9>1999&&$10>1){p=1}else{p=0}} (p==1){print}' SAM/$i/CLIP/preclust$def.fa > SAM/$i/CLIP/preclust$def\_2k.fasta &&  quast.py -t 8 SAM/$i/CLIP/preclust$def\_2k.fasta -R $LFR/Source/REF/fungi/merge/contig_8.fa -o SAM/$i/CLIP/preclust$def\_2k_quast &
done

 quast.py -t 32 TEST/Z1/summary.BI.megahit.clip.fasta -R $LFR/Source/REF/zymo/D6305.rRNA.fa -o TEST/Z1/summary.BI.megahit.clip.fasta_quast
```

**Check fungal reference ranges:
```bash
for i in Z{1,2,3};do
	blastn -num_threads 8 -db $LFR/Source/REF/zymo/D6305.rRNA.fa -query SAM/$i/summary.BI.megahit.clip.fasta -outfmt '6 std qlen qcov' -out SAM/$i/summary.BI.megahit.clip2zymo.m6 &
done
for i in M{4,5,6};do
	blastn -num_threads 8 -db $LFR/Source/REF/fungi/merge/contig_8.fa -query SAM/$i/summary.BI.megahit.clip.fasta -outfmt '6 std qlen qcov' -out SAM/$i/summary.BI.megahit.clip2ref8.m6 &
done

blastn -num_threads 32 -db $LFR/Source/REF/fungi/merge/contig_14.fa -query SAM/M5/summary.BI.megahit.clip.fasta -outfmt '6 std qlen' -out SAM/M5/summary.BI.megahit.clip2F14.m6
```
