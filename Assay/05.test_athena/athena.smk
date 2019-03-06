#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Athena process
# Call from:         .
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

#  Temp Main
rule all:
    input:
        expand("beadPool/{sample}.B.filter.dist", sample="NULL")

rule ATHENA_1_barcodeTag:
    input:
        i1 = "{sample}/input/fastp.sort.1.fq.gz",
        i2 = "{sample}/input/fastp.sort.2.fq.gz"
    output:  "{sample}/input/sort_atn_sp.fastq"
    shell:
        """
        perl -e 'open I1,\"pigz -p 4 -dc $ARGV[0]|\";open I2,\"pigz -p 4 -dc $ARGV[1]|\";
        while(<I1>){{
          @a = split("/",$_);@b = split("/",<I2>);
          $r1="$a[0]/$a[1] BC:Z:$a[1]-1\\n".<I1>.<I1>.<I1>;
          $r2="$b[0]/$b[1] BC:Z:$b[1]-1\\n".<I2>.<I2>.<I2>;
          unless($a[1] =~ /0000/){{ print \"$r1$r2\"}}
        }}' {input.i1} {input.i2} > {output}
        """

rule ATHENA_2_assemble:
    input:  "{sample}/input/sort_atn_sp.fastq"
    output: "{sample}/spadesMeta/contigs.fasta"
    params: "{sample}/spadesMeta"
    shell:
        "spades.py -t 48 --meta -o {params} --12 {input}"

rule ATHENA_3_index:
    input:
        idx = "{sample}/spadesMeta/contigs.fasta",
        inf = "{sample}/input/sort_atn_sp.fastq"
    output: "{sample}/align/reads.2.metaspades-contigs.bam"
    params: "{sample}/spadesMeta"
    shell:
        "bwa index {input}\n"
        "bwa mem -t 60 -C -p {input.idx} {input.inf}| "
        "samtools sort -@ 8 -o {output} - \n"
        "samtools index {output}"

rule ATHENA_4_config:
    input:
        fqs = "{sample}/input/sort_atn_sp.fastq",
        ctg = "{sample}/spadesMeta/contigs.fasta",
        bam = "{sample}/align/reads.2.metaspades-contigs.bam"
    output:   "{sample}/config.json"
    shell:
        """
        echo {{
          "input_fqs": "{input.fqs}",
          "ctgfasta_path": "{input.ctg}",
          "reads_ctg_bam_path": "{input.bam}"
        }} > {output}
        """

rule ATHENA_5_FinalShell:
    input: "{sample}/config.json"
    output:"{sample}/run.athena.sh"
    shell:
        "source activate athena_env\n"
        "athena-meta --config {input} --threads 64\n"

rule ATHENA_6_annotation:
    input:
        r1  = "{sample}/input/fastp.sort.1.fq.gz",
        r2  = "{sample}/input/fastp.sort.2.fq.gz",
        scaf = "{sample}/results/olc/flye-asm-1/scaffolds.fasta",
        refSeq  = "",
        refID   = ""
    params:
        oDir = "{sample}/reAlign",
        b6   = "{sample}/reAlign/FlyScaf2ref.blast6"
    output:"{sample}/FlyScaf2ref.blast6.anno"
    shell:
        "metabbq reAlign.sh 24 {input.r1} {input.r2} {input.scaf} {input.ref} {params.oDir}\n"
        "metabbq anno.pl {input.refID} {params.b6} > {output}"

rule ATHENA_7_circos:
    input:
        r1  = "{sample}/input/fastp.sort.1.fq.gz",
        r2  = "{sample}/input/fastp.sort.2.fq.gz",
        scaf = "{sample}/results/olc/flye-asm-1/scaffolds.fasta",
        refSeq  = "$mSeq2/instance/REF/ncbi.5Ref.fasta",
        refID   = "$mSeq2/instance/REF/ncbi.5Ref.faID"
    params:
        oDir = "{sample}/reAlign",
        b6   = "{sample}/reAlign/FlyScaf2ref.blast6"
    output:"{sample}/FlyScaf2ref.blast6.anno"
    shell:
        "metabbq circos.sh 100 97 reAlign $mSeq2/instance/REF/ncbi.5Ref.faID circos"
