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

rule TT1M_0_sketch:
    input:
        id = "{sample}/clean/fastp.sort.1.fq.idx",
        x1 = "{sample}/clean/fastp.sort.1.fq",
        x2 = "{sample}/clean/fastp.sort.2.fq",
        bb = "{sample}/clean/BB.stat"
    output:  "{sample}/TT1M/mash/bMin2.msh"
    params:
        pfx = "{sample}/TT1M/mash/bMin2",
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

rule TT1M_1_distRaw:
    input:  "{sample}/TT1M/mash/bMin2.msh"
    output: "{sample}/TT1M/mash/bMin2.raw.dist"
    shell:  "mash dist -p 24 -d 0.2 {input} {input}  > {output}"

rule TT1M_2a_distFilter:
    input:  "{sample}/TT1M/mash/bMin2.raw.dist"
    output: "{sample}/TT1M/mash/bMin2.clean.dist"
    params:
        min = config["p_dist_min"],
        max = config["p_dist_max"]
    shell:  "awk '($3>{params.min}&&$3<{params.max}){{print}}' {input}  > {output}"

rule TT1M_2b_distReNum:
    input:  "{sample}/TT1M/mash/bMin2.clean.dist"
    output:
        d = "{sample}/TT1M/mash/bMin2.bc.dist",
        m = "{sample}/TT1M/mash/bMin2.bc.map"
    shell:  "metabbq filterMASHoutput.pl -i  {input} -o {output.d} -M  {output.m}"

rule TT1M_3_convert:
    input:   "{sample}/TT1M/mash/bMin2.bc.dist"
    output:
        b  = "{sample}/TT1M/mash/bMin2.bc.bin",
        w  = "{sample}/TT1M/mash/bMin2.bc.w"
    shell:   "convert -i {input} -w {output.w} -o {output.b}"

rule TT1M_4_community:
    input:
        b  = "{sample}/TT1M/mash/bMin2.bc.bin",
        w  = "{sample}/TT1M/mash/bMin2.bc.w"
    output:   "{sample}/TT1M/mash/bMin2.bc.tree"
    shell:   "community {input.b} -w {input.w} -l -1 -v > {output}"

rule TT1M_5_clusterMAP:
    input:
        t = "{sample}/TT1M/mash/bMin2.bc.tree",
        m = "{sample}/TT1M/mash/bMin2.bc.map"
    output: "{sample}/TT1M/mash/bMin2.bc.tree.target.cluster"
    shell:
        "export maxLv=$[ `hierarchy {input.t} | wc -l` - 2 ]\n"
        "hierarchy {input.t} -l $maxLv > {input.t}.Lv$maxLv.lst\n"
        "metabbq change.id.pl -n {input.t}.Lv$maxLv.lst "
        "-m {input.m} -v -o {input.t}.Lv$maxLv.cluster\n"
        "ln -s bMin2.bc.tree.Lv$maxLv.cluster {output}"

rule TT1M_6_stat:
    input:
        i = "{sample}/clean/fastp.sort.1.fq.idx",
        c = "{sample}/TT1M/mash/bMin2.bc.tree.target.cluster"
    output: "{sample}/TT1M/mash/bMin2.bc.tree.target.cluster.count"
    shell:
        "perl -e 'open IN,\"<'{input.i}'\";while(<IN>){{chomp;@a=split;"
        "$num=($a[2]-$a[1]+1)/4;$HASH{{$a[0]}}=$num}};while(<>){{chomp;@a=split;"
        "$HB{{$a[0]}}{{R}}+=$HASH{{$a[1]}};$HB{{$a[0]}}{{C}}++}};"
        "foreach my $c (sort {{$a<=>$b}} keys %HB){{"
        "print \"$c\t$HB{{$c}}{{C}}\t$HB{{$c}}{{R}}\n\"}}' < {input.c} > {output}"

rule TT1M_7_main:
    input:
        count   = "{sample}/TT1M/mash/bMin2.bc.tree.target.cluster.count",
        cluster = "{sample}/TT1M/mash/bMin2.bc.tree.target.cluster"
    output:
        count   = "{sample}/TT1M/mash/bMin2.bc.tree.target.cluster.count.main",
        cluster = "{sample}/TT1M/mash/bMin2.bc.tree.target.cluster.main"
    params:
        rpc = config["p_rpc_min"],
        bpc = config["p_bpc_min"]
    shell:
        "awk '($2>={params.bpc}&&$3>={params.rpc}){{print}}' {input.count} | sort -k1,1n > {output.count}\n"
        "perl -e 'open IN,\"sort -k1,1n {output.count}|\";"
        "while(<IN>){{@a=split(/\\t/,$_);push @ids, $a[0]}}; $tag= shift @ids;"
        "while(<>){{@a=split /\\t/, $_; if($a[0] < $tag){{next}}elsif($a[0] == $tag){{"
         "print $_}}else{{$tag= shift @ids;print $_ if $a[0] == $tag }}}}' "
         "{input.cluster} |sort -k1,1n > {output.cluster}\n"



rule TT1M_8_write:
    input:
        s1 = "{sample}/clean/fastp.sort.1.fq",
        s2 = "{sample}/clean/fastp.sort.2.fq",
        bc = "{sample}/TT1M/mash/bMin2.bc.tree.target.cluster.main"
    output: "{sample}/TT1M/Assemble_mashBC.log"
    params:
        outDir = "{sample}/TT1M/Assemble_mashBC",
        threads = config["threads"]
    shell:
        "metabbq beadsWrite3.pl -b {input.bc} -f fq -t 4 -o {params.outDir} -v "
        "--r1 {input.s1}  --r2 {input.s2} > {output}\n"

rule TT1M_9_initASMsh:
    input:
        log = "{sample}/TT1M/Assemble_mashBC.log",
        bc = "{sample}/TT1M/mash/bMin2.bc.tree.target.cluster.count.main",
    output: "{sample}/TT1M/batch.assemble.BC.sh"
    params:
        sType  = config["sampleType"],
        samDir = "{sample}/TT1M",
        outDir = "{sample}/TT1M/Assemble_mashBC",
    shell:
        "for i in `sort -k2,2nr {input.bc} | cut -f1`; do "
        "echo metabbq bcPost.template.sh {params.sType} {params.samDir} mashBC $i BC;"
        "done > {params.samDir}/batch.assemble.BC.sh\n"

rule TT1M_10_circos:
    input:
        log = "{sample}/TT1M/Assemble_mashBC.log",
        bc = "{sample}/TT1M/mash/bMin2.bc.tree.target.cluster.count.main",
    output: "{sample}/TT1M/batch.circos.BC.sh"
    params:
        refDB  = config["REF_GEN"],
        refID  = config["REF_ID"],
        zoom   = config["p_cc_zoom"],
        cut    = config["p_cc_cut"],
        inDir  = "{sample}/TT1M/Assemble_mashBC",
        outDir = "{sample}/TT1M/circos"
    shell:
        """
        for i in `ls {params.inDir}`; do
          echo metabbq reAlign.sh 8 {params.inDir}/$i/sort.{{1,2}}.fq \
          {params.inDir}/$i/scaffolds.fasta {params.refDB} {params.outDir}/$i/reAlign
          echo metabbq circos.sh {params.zoom} {params.cut} {params.outDir}/$i/reAlign \
          {params.refID} {params.outDir}/$i/circos;
        done > {output}
        """
