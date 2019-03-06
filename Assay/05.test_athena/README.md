# CODE NOTE

Athena guide: https://github.com/abishara/athena_meta

1.prepare input
```sh
mkdir -p APR842_00_00/input
ln -sf $mSeq2/instance/clean/APR842_00_00.fps.1.fq.gz APR842_00_00/input/fastp.sort.1.fq.gz
ln -sf $mSeq2/instance/clean/APR842_00_00.fps.2.fq.gz APR842_00_00/input/fastp.sort.2.fq.gz

cd APR842_00_00

#Transform barcode format to fit athena
perl -e 'open I1,"pigz -p 4 -dc $ARGV[0]|";open I2,"pigz -p 4 -dc $ARGV[1]|";
while(<I1>){
  @a = split("/",$_);@b = split("/",<I2>);
  $r1="$a[0]/$a[1] BC:Z:$a[1]-1\n".<I1>.<I1>.<I1>;
  $r2="$b[0]/$b[1] BC:Z:$b[1]-1\n".<I2>.<I2>.<I2>;
  unless($a[1] =~ /0000/){ print "$r1$r2"}
}' input/fastp.sort.1.fq.gz input/fastp.sort.2.fq.gz \
> input/sort_atn_sp.fastq
```

2.Assemble
```sh
spades.py -t 48 --meta -o spadesMeta --12 input/sort_atn_sp.fastq
```

3.Create a bwa index
```bash
mkdir -p index align
ln -sf ../spadesMeta/contigs.fasta index/contigs.fasta
bwa index index/contigs.fasta
bwa mem -t 60 -C -p index/contigs.fasta input/sort_atn_sp.fastq| \
samtools sort -@ 8 -o align/reads.2.metaspades-contigs.bam -
samtools index align/reads.2.metaspades-contigs.bam
```

4.Setup a configuration file
create a config.json and following content:
```
{
    "input_fqs": "input/sort_atn_sp.fastq",
    "ctgfasta_path": "index/contigs.fasta",
    "reads_ctg_bam_path": "align/reads.2.metaspades-contigs.bam"
}
```

5.Run Athena
```
source activate athena_env
athena-meta --config config.json --threads 64
```

# with athena.smk
```bash
ln -s ../../Results/DEC560_d10_00
#preview
snakemake -np DEC560_d10_00/athena/run.athena.sh
#run
snakemake -p DEC560_d10_00/athena/run.athena.sh
sh DEC560_d10_00/athena/run.athena.sh
```

#Assessment
```
quast.py DEC560_d10_00/flye-asm-1/scaffolds.fasta \
        -1 DEC560_d10_00/input/fastp.sort.1.fq.gz \
        -2 DEC560_d10_00/input/fastp.sort.2.fq.gz \
        -o DEC560_d10_00/quast_test_output \
        -t 24
```
