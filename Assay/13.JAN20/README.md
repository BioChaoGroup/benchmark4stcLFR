# Test data generated at Jan 2020

# Sample information

Above information were save in sheet `dataInfo.tsv`.

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

Run following steps:

```bash
metabbq smk -j --resources mem=120 -np SAM/{{U,S,O,M,Z}{1,2,3},{S,O,M,Y}{4,5,6},{S,O,Z}{7,8,9},{S,O,M}{A,B,C}}/Assemble_BI/ID.lst
```

RCA isolate assemble:
```bash
metabbq smk -j -np SAM/{{U,S,O,M,Z}{1,2,3},{S,O,M,Y}{4,5,6}}/summary.BI.megahit.rRNA.fasta

metabbq smk -j -npk SAM/{{U,O,M,Z}{1,2,3},S{1,2},{S,O,M,Y}{4,5,6}}/VSEARCH/barrnap.RSUs.fasta

```

**fungi** needs specific parameters for RNA units prediction:
```bash
rename rRNA bac SAM/{S,O,M,Y}{4,5,6}/summary.BI.megahit.rRNA.fasta{,.barrnap}
metabbq smk --config sampleType=F -j -npk SAM/{S,O,M,Y}{4,5,6}/summary.BI.megahit.rRNA.fasta

```

**stat barrnap RSU seq info**
```bash
metabbq stat RSU -l SAM/S1/VSEARCH/barrnap.LSU.fasta -s SAM/S1/VSEARCH/barrnap.SSU.fasta -o SAM/S1/stat/barrnap.RSU.stat
metabbq smk -j -npk SAM/{{U,O,S,M,Z}{1,2,3},{S,O,M,Y}{4,5,6}}/stat/barrnap.RSU.stat
```
