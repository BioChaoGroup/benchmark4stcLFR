# Test data generated at Jan 2020

# Sample information

Above information were save in sheet `dataInfo.tsv`.

# Analysis preparation
Organize directories. `ORI(GIN)` were organized by sequencing events; `SAM(PLE)` were organized by sample name.
```bash
mkdir -p ../../Results/JAN20/{ORI,SAM}
ln -s ../../Results/JAN20/ORI ./
```
Then generate directories for each sequence data:
```bash
awk '{system("mkdir -p ORI/"$4"/input");
system("ln -s "$5" ORI/"$4"/input/rawSeq_1.fq.gz");
sub("_1.fq","_2.fq",$5);
system("ln -s "$5" ORI/"$4"/input/rawSeq_2.fq.gz");
}' dataInfo.tsv
```

Run QC step:
```bash
#for RCA bac
metabbq smk -j -np ORI/FF{1,2}_{1,3,4,5,7,8,17,22}/clean/BB.stat
#for RCA fungi
metabbq smk -j -np ORI/FF{5,6}_{9,10,11,12,13,20,15,16}/clean/BB.stat
#for Non-RCA bac & fungi
metabbq smk -j -np ORI/FF{9_{1,3,17,4,5,22,7,8},A_{9,10,11,12,13,20,15,16}}/clean/BB.stat

```

After `fastp`, we can merge samples together:
**RCA bacteria**
```bash
ln -s ../../Results/JAN20/SAM ./
mkdir -p SAM/{U,S,O}{1,2}/{clean,.tmp.{1,2}}
for i in {1,2};do for j in {1,2};do
echo cat ORI/FF$i\_{1,17}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/U$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/U$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{3,4}/clean/fastp.sort.$j.fq  \| paste - - - - \| sort -T SAM/S$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/S$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{5,22}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/O$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/O$i/clean/fastp.sort.$j.fq \&
done;done
metabbq smk -j -np SAM/{U,S,O}{1,2}/clean/BB.stat
#
for i in {1,2};do
	ln -s ../../ORI/FF$i\_7/clean SAM/M$i/
	ln -s ../../ORI/FF$i\_8/clean SAM/Z$i/
done
```
**RCA fungi**
```bash
mkdir -p SAM/{S,O,M,Y}{5,6}/{clean,.tmp.{1,2}}
for i in {5,6};do for j in {1,2};do
echo cat ORI/FF$i\_{9,10,11}/clean/fastp.sort.$j.fq  \| paste - - - - \| sort -T SAM/S$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/S$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{12,13,20}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/O$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/O$i/clean/fastp.sort.$j.fq \&
done;done
metabbq smk -j -np SAM/{S,O}{5,6}/clean/BB.stat
#
for i in {5,6};do
	ln -s ../../ORI/FF$i\_16/clean SAM/M$i/
	ln -s ../../ORI/FF$i\_15/clean SAM/Y$i/
done
```
**Non-RCA bacteria**
```bash
mkdir -p SAM/{S,O,Z}9/{clean,.tmp.{1,2}}
for i in {9};do for j in {1,2};do
echo cat ORI/FF$i\_{1,3,17}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/S$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/S$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{4,5,22}/clean/fastp.sort.$j.fq  \| paste - - - - \| sort -T SAM/O$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/O$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{7,8}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/Z$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/Z$i/clean/fastp.sort.$j.fq \&
done;done
metabbq smk -j -np SAM/{S,O,Z}9/clean/BB.stat
#
```
**Non-RCA fungi**
```bash
mkdir -p SAM/{S,O,M}A/{clean,.tmp.{1,2}}
for i in {A};do for j in {1,2};do
echo cat ORI/FF$i\_{9,10,11}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/S$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/S$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{12,13,20}/clean/fastp.sort.$j.fq  \| paste - - - - \| sort -T SAM/O$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/O$i/clean/fastp.sort.$j.fq \&
echo cat ORI/FF$i\_{15,16}/clean/fastp.sort.$j.fq \| paste - - - - \| sort -T SAM/M$i/.tmp.$j -k2,2 -t \"/\" \| tr \"\\t\" \"\\n\" \> SAM/M$i/clean/fastp.sort.$j.fq \&
done;done
metabbq smk -j -np SAM/{S,O,M}A/clean/BB.stat
#
```
Run following steps:
```bash

```
