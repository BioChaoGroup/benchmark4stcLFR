#!/usr/bin/env snakemake
# (c) 2016 - 2019 Chao IN-HORSE SHARE ONLY
# ===================================================================
# Description:       Beads Cluster Generation
# Call from:         ../Snakefile
# Author:            Chao | fangchao@genomics.cn
# ===================================================================

# Init
## Read config
configfile: "config.yaml"

rule AS1M_0_sketch:
    input:
        id = "{sample}/clean/fastp.sort.1.fq.idx",
        x1 = "{sample}/clean/fastp.sort.1.fq",
        x2 = "{sample}/clean/fastp.sort.2.fq",
        bb = "{sample}/clean/BB.stat"
    output:  "{sample}/{dtag}/mash/bMin2.msh"
    params:
        pfx = "{sample}/{dtag}/mash/bMin2",
        minR= config['p_cluster_minR'],
        maxR= config['p_cluster_maxR'],
        topB= config['p_cluster_topB'],
        ranP= config['p_cluster_ranP']
    shell:
        "export maxC={params.maxR} \n export minC={params.minR}\n"
        "echo choose minc = $minC , maxc = $maxC \n"
        "metabbq binWrite fqpick -x {input.id} -c {params.minR} -m {params.maxR} -b {params.topB} -r {params.ranP} -i {input.x1} -o {params.pfx}.sort.1.fq & \n"
        "metabbq binWrite fqpick -x {input.id} -c {params.minR} -m {params.maxR} -b {params.topB} -r {params.ranP} -i {input.x2} -o {params.pfx}.sort.2.fq & \n"
        "wait && mash sketch -r -B {params.pfx}.sort.1.fq {params.pfx}.sort.2.fq -o {params.pfx}"


rule AS1M_1_distRaw:
    input:  "{sample}/{dtag}/mash/bMin2.msh"
    output: "{sample}/{dtag}/mash/bMin2.raw.dist"
    shell:  "mash dist -p 24 -d 0.2 {input} {input}  > {output}"

rule stat_1_beadAnno:
    input:
        rawDist = "{sample}/{dtag}/mash/bMin2.raw.dist",
        uniqAnn = "{sample}/{dtag}/mash/bMin2.raw.dist"
    output: "{sample}/{dtag}/mash/bMin2.raw.dist.anno.stat"
    shell:
        "perl beadAnno.pl {input.rawDist} {input.uniqAnn} > {output}"

rule AS1M_2a_distFilter:
    input:  "{sample}/{dtag}/mash/bMin2.raw.dist"
    output: "{sample}/{dtag}/mash/bMin2.clean.dist"
    params:
        min = config["p_dist_min"],
        max = config["p_dist_max"]
    shell:  "awk '($3>{params.min}&&$3<{params.max}){{print}}' {input}  > {output}"

rule AS1M_2b_distReNum:
    input:  "{sample}/{dtag}/mash/bMin2.clean.dist"
    output:
        d = "{sample}/{dtag}/mash/bMin2.bc.dist",
        m = "{sample}/{dtag}/mash/bMin2.bc.map"
    shell:  "metabbq filterMASHoutput.pl -i  {input} -o {output.d} -M  {output.m}"

rule AS1M_3_convert:
    input:   "{sample}/{dtag}/mash/bMin2.bc.dist"
    output:
        b  = "{sample}/{dtag}/mash/bMin2.bc.bin",
        w  = "{sample}/{dtag}/mash/bMin2.bc.w"
    shell:   "convert -i {input} -w {output.w} -o {output.b}"

rule AS1M_4_community:
    input:
        b  = "{sample}/{dtag}/mash/bMin2.bc.bin",
        w  = "{sample}/{dtag}/mash/bMin2.bc.w"
    output:   "{sample}/{dtag}/mash/bMin2.bc.tree"
    shell:   "community {input.b} -w {input.w} -l -1 -v > {output}"

if config["p_dist_lv"]:
    output5 = "{sample}/{dtag}/mash/lv" + str(config["p_dist_lv"]) + "/tree.cluster"
    outdir5 = "{sample}/{dtag}/mash/lv" + str(config["p_dist_lv"])
    rule AS1M_5_clusterMAP:
        input:
            t = "{sample}/{dtag}/mash/bMin2.bc.tree",
            m = "{sample}/{dtag}/mash/bMin2.bc.map"
        output: output5
        params:
            dir = outdir5,
            lv  = config["p_dist_lv"]
        shell:
            "hierarchy {input.t} -l {params.lv} > {params.dir}/tree.lst\n"
            "metabbq change.id.pl -n {params.dir}/tree.lst "
            "-m {input.m} -v -o {output}"

    c6 = "{sample}/{dtag}/mash/lv" + str(config["p_dist_lv"]) + "/tree.cluster"
    o6 = "{sample}/{dtag}/mash/lv" + str(config["p_dist_lv"]) + "/tree.cluster.count"
    rule AS1M_6_stat:
        input:
            i = "{sample}/clean/fastp.sort.1.fq.idx",
            c = c6
        output: o6
        shell:
            "perl -e 'open IN,\"<'{input.i}'\";while(<IN>){{chomp;@a=split;"
            "$num=($a[2]-$a[1]+1)/4;$HASH{{$a[0]}}=$num}};while(<>){{chomp;@a=split;"
            "$HB{{$a[0]}}{{R}}+=$HASH{{$a[1]}};$HB{{$a[0]}}{{C}}++}};"
            "foreach my $c (sort {{$a<=>$b}} keys %HB){{"
            "print \"$c\t$HB{{$c}}{{C}}\t$HB{{$c}}{{R}}\n\"}}' < {input.c} > {output}"

    oct7 = "{sample}/{dtag}/mash/lv" + str(config["p_dist_lv"]) + "/tree.cluster.count.main"
    ocl7 = "{sample}/{dtag}/mash/lv" + str(config["p_dist_lv"]) + "/tree.cluster.main"
    rule AS1M_7_main:
        input:
            count   = o6,
            cluster = c6
        output:
            count   = oct7,
            cluster = ocl7,
            linkc   = "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster.main",
            linkt   = "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster.count.main"
        params:
            rpc = config["p_rpc_min"],
            bpc = config["p_bpc_min"],
            lv = str(config["p_dist_lv"])
        shell:
            "metabbq clusterHelper main -a -i {input.cluster} -c {input.count} "
            "-o {output.cluster} -s {output.count} -v\n"
            "ln -sf lv{params.lv}/tree.cluster.main {output.linkc}\n"
            "ln -sf lv{params.lv}/tree.cluster.count.main {output.linkt}"

else:
    rule AS1M_5_clusterMAP:
        input:
            t = "{sample}/{dtag}/mash/bMin2.bc.tree",
            m = "{sample}/{dtag}/mash/bMin2.bc.map"
        output: "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster"
        shell:
            "export maxLv=$[ `hierarchy {input.t} | wc -l` - 2 ]\n"
            "hierarchy {input.t} -l $maxLv > {input.t}.Lv$maxLv.lst\n"
            "metabbq change.id.pl -n {input.t}.Lv$maxLv.lst "
            "-m {input.m} -v -o {input.t}.Lv$maxLv.cluster\n"

    rule AS1M_6_stat:
        input:
            i = "{sample}/clean/fastp.sort.1.fq.idx",
            c = "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster"
        output: "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster.count"
        shell:
            "perl -e 'open IN,\"<'{input.i}'\";while(<IN>){{chomp;@a=split;"
            "$num=($a[2]-$a[1]+1)/4;$HASH{{$a[0]}}=$num}};while(<>){{chomp;@a=split;"
            "$HB{{$a[0]}}{{R}}+=$HASH{{$a[1]}};$HB{{$a[0]}}{{C}}++}};"
            "foreach my $c (sort {{$a<=>$b}} keys %HB){{"
            "print \"$c\t$HB{{$c}}{{C}}\t$HB{{$c}}{{R}}\n\"}}' < {input.c} > {output}"

    rule AS1M_7_main:
        input:
            count   = "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster.count",
            cluster = "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster"
        output:
            count   = "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster.count.main",
            cluster = "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster.main"
        params:
            rpc = config["p_rpc_min"],
            bpc = config["p_bpc_min"]
        shell:
            "metabbq clusterHelper main -r {params.rpc} -b {params.bpc} "
            "-i {input.cluster} -c {input.count} "
            "-o {output.cluster} -s {output.count} -v\n"


rule AS1M_8_write:
    input:
        s1 = "{sample}/clean/fastp.sort.1.fq",
        s2 = "{sample}/clean/fastp.sort.2.fq",
        bc = "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster.main"
    output: "{sample}/{dtag}/Assemble_mashBC.log"
    params:
        outDir = "{sample}/{dtag}/Assemble_mashBC",
        threads = config["threads"]
    shell:
        "metabbq beadsWrite3.pl -b {input.bc} -f fq -t 4 -o {params.outDir} -v "
        "--r1 {input.s1}  --r2 {input.s2} > {output}\n"

rule AS1M_9_initASMsh:
    input:
        log = "{sample}/{dtag}/Assemble_mashBC.log",
        bc = "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster.count.main",
    output: "{sample}/{dtag}/batch.assemble.BC.sh"
    params:
        samDir = "{sample}/{dtag}",
        outDir = "{sample}/{dtag}/Assemble_mashBC",
        mode   = config["method"]["assemble"]["mode"],
        ref1   = config["REF_FA"],
        ref2   = config["REF_FA"]
    shell:
        "for i in `sort -k2,2nr {input.bc} | cut -f1`; do "
        "echo sh ./bcPost.template.sh {params.mode} {params.samDir} mashBC $i BC "
        "{params.ref1} {params.ref2} ;"
        "done > {params.samDir}/batch.assemble.BC.sh\n"

outXs = "{sample}/{dtag}/summary." + str(config["method"]["assemble"]["mode"]) + ".contig.tsv"
outXa= "{sample}/{dtag}/summary." + str(config["method"]["assemble"]["mode"]) + ".contig.fasta"
rule AS1M_X_summary:
    input: "{sample}/{dtag}/batch.assemble.BC.sh"
    output:
        stat = outXs,
        fa   = outXa
    params:
        outDir = "{sample}/{dtag}/Assemble_mashBC",
        mAsm   = config["p_asm_min"],
        mode   = config["method"]["assemble"]["mode"]
    shell:
        """
        for i in `ls {params.outDir}`;do
          metabbq clusterHelper asmPick -l {params.mAsm} -b $i -i {params.outDir}/$i/{params.mode}/final.contigs.fa
        done > {output.fa}
        for i in `ls {params.outDir}`;do
          awk -v bc=$i '/^>/{{print bc"\\t"$0}}' {params.outDir}/$i/{params.mode}/final.contigs.fa | sed 's/>//;s/ flag=/\\t/;s/ multi=/\\t/;s/ len=/\\t/'
        done > {output.stat}
        """


rule AS1M_ext_circos:
    input:
        log = "{sample}/{dtag}/Assemble_mashBC.log",
        bc = "{sample}/{dtag}/mash/bMin2.bc.tree.target.cluster.count.main",
    output: "{sample}/{dtag}/batch.circos.BC.sh"
    params:
        refDB  = config["REF_GEN"],
        refID  = config["REF_ID"],
        zoom   = config["p_cc_zoom"],
        cut    = config["p_cc_cut"],
        inDir  = "{sample}/{dtag}/Assemble_mashBC",
        outDir = "{sample}/{dtag}/circos"
    shell:
        """
        for i in `ls {params.inDir}`; do
          echo metabbq reAlign.sh 8 {params.inDir}/$i/sort.{{1,2}}.fq \
          {params.inDir}/$i/scaffolds.fasta {params.refDB} {params.outDir}/$i/reAlign
          echo metabbq circos.sh {params.zoom} {params.cut} {params.outDir}/$i/reAlign \
          {params.refID} {params.outDir}/$i/circos;
        done > {output}
        """
