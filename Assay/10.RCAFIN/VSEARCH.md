Fungal mock strains:

| ID   | Prefix | Tag  | Taxonomy                    | Alternative                |
| ---- | ------ | ---- | --------------------------- | -------------------------- |
| 158  | JOMC   | Au   | Aspergillus ustus           | Aspergillus puniceus       |
| 159  | BCGH   | Tk   | Trichoderma koningii        | Trichoderma pseudokoningii |
| 160  | JQFZ   | Pe   | Penicillium expansum        |                            |
| 161  | AACD   | An   | Aspergillus nidulans        |                            |
| 162  | JMSF   | Pc   | Penicillium chrysogenum     |                            |
| 163  | MBDJ   | Tl   | Trichoderma longibrachiatum |                            |
| 164  | AAIL   | Tr   | Trichoderma viride          | Trichoderma reesei         |
| 165  | ABDG   | Ta   | Trichoderma atroviride      |                            |

> 159:160:161=100:10:1

```bash
SAM="RCA197M01_4"
threads=32
input=$SAM/summary.BI.megahit.contig.CLIP4.fasta
pfx=$SAM/VSEARCH/summary.BI.megahit.contig.CLIP4

mkdir -p $SAM/ITSx
awk -F '_' '(FNR%2==1){if($8>1999){t=1}else{t=0}};(t==1){print}' $input > $pfx.2k.fa
ITSx --cpu 4 -t F -i $pfx.2k.fa -o $SAM/ITSx/contig.CLIP.2k.ITSx

awk -F '_' '(FNR%2==1){if($8>499){t=1}else{t=0}};(t==1){print}' $input > $pfx.500.fa
ITSx --cpu 4 -t F -i $pfx.500.fa -o $SAM/ITSx/contig.CLIP.500.ITSx


mkdir -p $SAM/VSEARCH
echo "Dereplicate (non-singleton)"
vsearch --threads $threads --derep_fulllength $input  \
--output $pfx.derep.full.fasta \
-uc $pfx.derep.full.uc --fasta_width 0

# Precluster at 99.9% before chimera detection
vsearch --threads $threads --cluster_fast $pfx.derep.full.fasta \
--id 0.999 --strand both --fasta_width 0 \
--uc $pfx.derep.preclustered.uc --minuniquesize 1 \
--centroids $pfx.derep.preclustered.fasta

# De novo chimera detection
vsearch --threads $threads --uchime_denovo $pfx.derep.preclustered.fasta \
--fasta_width 0 \
--nonchimeras $pfx.derep.denovo.nonchimeras.fasta

pct=99
# Cluster nonchimeras at $pct% and relabel with OTU_n, generate OTU table
vsearch --threads $threads --cluster_fast  $pfx.derep.denovo.nonchimeras.fasta \
--id 0.$pct --strand both --fasta_width 0 \
--relabel LFR_ --relabel_keep \
--centroids $pfx.derep.nonchimeras.LFRs.$pct.fasta \
--uc $pfx.derep.nonchimeras.LFRs.$pct.uc

grep ">" $pfx.derep.nonchimeras.LFRs.$pct.fasta|sed 's/>//;s/ /\t/g' \
> $pfx.derep.nonchimeras.LFRs.$pct.IDs

blastn -num_threads $threads -perc_identity 98 -qcov_hsp_perc 50 -db ../../Source/REF/MIX7.fa -query $pfx.derep.nonchimeras.LFRs.$pct.fasta \
-outfmt '6 std qlen' -out $pfx.derep.nonchimeras.LFRs.98.b6

blastn -num_threads $threads -perc_identity 98 -qcov_hsp_perc 50 -db ../../Source/REF/zymo/D6305.genomes -query $pfx.LFRs.fasta \
-outfmt '6 std qlen' -out $pfx.LFRs.zymo.genome.b6

blastn -num_threads $threads -perc_identity 98 -qcov_hsp_perc 50 -db ../../Source/REF/zymo/D6305.markers.fa -query $pfx.LFRs.fasta -outfmt '6 std qlen' -out $pfx.LFRs.markers.b6

blastn -num_threads $threads -perc_identity 98 -qcov_hsp_perc 50 -db ../../Source/REF/zymo/D6305.rRNA.fa -query $pfx.LFRs.fasta -outfmt '6 std qlen' -out $pfx.LFRs.rRNA.b6

metabbq findRegion2 -i $pfx.LFRs.zymo.genome.b6 -o $pfx.LFRs.zymo.genome.b6.more

#mapping all beads to the FLRs
bwa index $pfx.LFRs.fasta
bwa mem -t $threads $pfx.LFRs.fasta \
$SAM/clean/fastp.sort.{1,2}.fq > $pfx.LFRs.bwa

metabbq beadStat sam -i $pfx.LFRs.bwa \
-o $pfx.LFRs.b2r -v

metabbq beadStat sam2b -i $pfx.LFRs.b2r -o $pfx.LFRs.bb1
 
awk '$5>=0.5' $pfx.derep.nonchimeras.LFRs.$pct.bb1 > $pfx.derep.nonchimeras.LFRs.$pct.bb2

awk '$5>0.5{print $4}' $pfx.derep.nonchimeras.LFRs.$pct.bb1|sort|uniq -c \
> $pfx.derep.nonchimeras.LFRs.$pct.prop

perl -e 'while(<>){chomp;@a=split;
$HS{$a[3]}{$a[6]}{BB} ++; 
$HS{$a[3]}{$a[6]}{US} += $a[1]; $HS{$a[3]}{$a[6]}{MS} += $a[2];
$HS{$a[3]}{$a[6]}{U1} += $a[1] * $a[4]; 
$HS{$a[3]}{$a[6]}{M2} += $a[1] * $a[7] + $a[2] * $a[8];
};foreach my $h1(sort keys %HS){foreach my $h2 (sort keys %{$HS{$h1}}){
printf ("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\n",
$h1,$h2,$HS{$h1}{$h2}{BB},$HS{$h1}{$h2}{US},$HS{$h1}{$h2}{MS},
($HS{$h1}{$h2}{US}==0)?0:$HS{$h1}{$h2}{U1}/$HS{$h1}{$h2}{US},
($HS{$h1}{$h2}{MS}+$HS{$h1}{$h2}{US}==0)?0:$HS{$h1}{$h2}{M2}/($HS{$h1}{$h2}{MS}+$HS{$h1}{$h2}{US}))}}' \
< $pfx.derep.nonchimeras.LFRs.$pct.bb1 > $pfx.derep.nonchimeras.LFRs.$pct.bbs

#zymo
SAM=RCA197M01_2
vsearch --threads $threads --allpairs_global $SAM/VSEARCH/contig.LFRs.fasta --uc $SAM/VSEARCH/contig.LFRs.allpairs.uc --iddef 0 --id 0.95

perl -e 'while(<>){chomp;@a=split;$HS{$a[9]}++}; open IN,"< '$SAM'/VSEARCH/contig.LFRs.fasta"; while(<IN>){@a=split />| /;$seq=<IN>;unless(defined $HS{$a[1]}){print $_; print $seq}}' < $SAM/VSEARCH/contig.LFRs.allpairs.uc > $SAM/VSEARCH/contig.LFRs.filter.fasta


## backup ##########################################################
awk '($1!=p){print;p=$1}' $pfx.derep.denovo.nonchimeras.b6
awk '($3>99&&$4/$13>0.8){print}' $pfx.derep.denovo.nonchimeras.b6 > $pfx.derep.denovo.nonchimeras.filter.b6

ITSx --cpu 4 -t F -i $pfx.derep.preclustered.comp.fasta --detailed_results T \
-o $pfx.derep.preclustered.comp.ITSx
mafft  --auto --reorder $pfx.derep.denovo.nonchimeras.fasta > $pfx.derep.denovo.nonchimeras.mafft.fasta


#zymo
SAM=RCA199M03_1
SAM=RCA199M03_2
metabbq RCACLIPS -a AAGTCGGACGCTGATAAGGTCGCCATGCCTCTCAGTACTCCGACTT -fwd AGAGTTTGATCATGGCTCAG -rev CTTAGATGCCTTCAGC -i $SAM/summary.BI.megahit.contig.fasta -o $SAM/summary.BI.megahit.clip.bac.fasta -v 2> $SAM/summary.BI.megahit.clip.bac.fasta.log

metabbq RCACLIPS -a AAGTCGGACGCTGATAAGGTCGCCATGCCTCTCAGTACTCCGACTT -fwd AGAGTTTGATCATGGCTCAG -rev AAGGAGGTGATCCAGCCGCA -i $SAM/summary.BI.megahit.contig.fasta -o $SAM/summary.BI.megahit.clip.bac2.fasta -v 2> $SAM/summary.BI.megahit.clip.bac2.fasta.log


barrnap --threads 4 --reject 0.1 $SAM/VSEARCH/contig.LFRs.fasta --outseq $SAM/VSEARCH/contig.LFRs.barrnap

barrnap --threads 4 --reject 0.05 $SAM/summary.BI.megahit.clip.bac.fasta --outseq $SAM/summary.BI.megahit.clip.bac.barrnap



```

#zymo non-RCA

```bash
SAM=RCA197M01_8
SAM=RCA197M02_8
SAM=RCA199M03_8

vsearch --threads 8 --fastq_mergepairs $SAM/clean/fastp.sort.1.fq \
--reverse $SAM/clean/fastp.sort.2.fq --fasta_width 0 \
--fastqout $SAM/clean/fastp.sort.merged.fq \
--fastqout_notmerged_fwd $SAM/clean/fastp.sort.unmerged.F.fq \
--fastqout_notmerged_rev $SAM/clean/fastp.sort.unmerged.R.fq
### old
Merging reads 100%
  67679483  Pairs
  41246048  Merged (60.9%)
  26433435  Not merged (39.1%)

Pairs that failed merging due to various reasons:
  17886324  too few kmers found on same diagonal
     76401  multiple potential alignments
     61213  too many differences
   3700231  alignment score too low, or score drop to high
    160781  overlap too short
   4548485  staggered read pairs

Statistics of all reads:
     91.19  Mean read length

Statistics of merged reads:
    108.65  Mean fragment length
     43.40  Standard deviation of fragment length
      0.26  Mean expected error in forward sequences
      0.48  Mean expected error in reverse sequences
      0.22  Mean expected error in merged sequences
      0.19  Mean observed errors in merged region of forward sequences
      0.37  Mean observed errors in merged region of reverse sequences
      0.56  Mean observed errors in merged region
### new
Merging reads 100%
  67583804  Pairs
  23101371  Merged (34.2%)
  44482433  Not merged (65.8%)

Pairs that failed merging due to various reasons:
  17886324  too few kmers found on same diagonal
    166347  multiple potential alignments
     61213  too many differences
   3699867  alignment score too low, or score drop to high
    160781  overlap too short
  22507901  staggered read pairs

Statistics of all reads:
    100.00  Mean read length

Statistics of merged reads:
    141.25  Mean fragment length
     25.71  Standard deviation of fragment length
      0.31  Mean expected error in forward sequences
      0.62  Mean expected error in reverse sequences
      0.38  Mean expected error in merged sequences
      0.18  Mean observed errors in merged region of forward sequences
      0.45  Mean observed errors in merged region of reverse sequences
      0.62  Mean observed errors in merged region
###

metabbq beadsWrite3.pl -x --r1 $SAM/clean/fastp.sort.merged.fq -v 
awk '$1!~/0000/{print ($3-$2+1)/4}' $SAM/clean/fastp.sort.merged.fq.idx|sort|uniq -c|sort -k2,2nr|awk '{b=b+$1;print $0"\t"b}' > $SAM/clean/merged.BB.stat

mkdir -p $SAM/mash/split

metabbq binWrite fqpick -x $SAM/clean/fastp.sort.merged.fq.idx -c 10 -m 30 -b 0 -r 0 -i $SAM/clean/fastp.sort.merged.fq -o $SAM/mash/merged.bMin10.fq
metabbq beadsWrite3.pl -x --r1 $SAM/mash/merged.bMin10.fq

metabbq binWrite fqSplit -x $SAM/mash/merged.bMin10.fq.idx -p 100 -i $SAM/mash/merged.bMin10.fq -o $SAM/mash/split/merged.

for i in {00..99};do echo mash sketch -p 6 -r -B $SAM/mash/split/merged.$i $SAM/mash/split/merged.$i -o $SAM/mash/split/merged.$i;done > batch.sketch.sh

for i in {00..99};do for j in {00..99}; do echo "mash dist -p 8 -d 0.1 $SAM/mash/split/merged.$i.msh $SAM/mash/split/merged.$j.msh|gzip > $SAM/mash/split/merged.$i\_$j.dist.gz"; done ; done > batch.mash.sh

zcat $SAM/mash/split/merged.{00..99}_{00..99}.dist.gz > $SAM/mash/bMin2.raw.dist

zcat $SAM/mash/split/merged.{00..99}_{00..99}.dist.gz | awk '($3>0.02&&$3<0.04){print}'  > $SAM/mash/bMin2.clean.dist
```



