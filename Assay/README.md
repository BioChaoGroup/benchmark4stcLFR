# Technical Note

### Landscape

![fig1. landscape of the pipeline](../Results/general/fig.pipeline.png)

**Fig.1** Landscape of pipeline using stLFR barcodes information

Figure 1 shows a pipeline of how we detecting barcodes and used it to improve the assembly results. The core process is to treat sequences in bead level and cluster them according to the similarity of reads into each other beads. The clustered subset is called Bead-barcodes Cluters. Sequence belongint to each cluter was assembled and annotated parallelled.

### Method

#### Sampling

For Fungus, we prepared a mock sample purched from ZYMO.

**Sequencing** was performed by DIPSEQ T1 (from MGI), PE100 mode.

**Quality control** was implented in a modified **[fastp](https://github.com/Scelta/fastp)**, along with the function to detect bead barcodes.

**VSEARCH** was used to compute duplicates and reads pre-clusters.

Codes and pipelines are availiable by accessing the project **metaseq**.

### Results

#### Bead Barcodes detection

The **bead barcodes** (BBs) are idetntified by 3 sequeces of 10 bp attached in the tail of each read 2.  

Generally, more than 85% of reads are successfully idtentiable with a barely offests (< 1% )  (Fig.2).

![fig.2](./01.stLFR_barcode_split/BBdetectionPct.png)

> **Fig.2 Reads distribution of bead barcodes detection**.
> Data tested from 1 millon randomly picked reads.

We totally identified Among the identified **68645** unique BBs. However, about **40%** of them are only attached by one pair-read (fig3. black).

![fig.3](./04.test_vsearch/BBwithReadsin.png)

>**Fig.3 Distribution of bead barcodes with different number of reads attached to it**.
>Black, all reads included; red, only dereplicated reads included

Since amplicon also induce duplicates, we also attempted to remove them. Among **1.94 m** merged paired-reads, **242k** duplicate centroids with identified by VSEARCH, **1.46 m** are removed as duplicates.

However, we **found reads with different BB belonging to a specific cnetroid**, which means they a "false-true duplicates". After split those centroids with unique BB,  **314086** duplicate centroids found and **71603** reads saved .  Besides, the dereplicated data shows that the proportion of singloton BB are increased to **55.7%** (fig3. red). More than half of the detected beads only got one single read with it.

```bash
# CODE NOTE:
cd instance/APR842_00/VSEARCH
# return duplicates without BB info
grep "^>" read.merge.fasta|wcl
1943671
grep "^S" read.merge.derep.uc|grep -v "0000" |wc -l
242483
# return duplicates seperated by BB
grep -v "0000" read.merge.derep.bead.uc.S|wc -l
314086
```

#### Bead Barcodes cluster

Those "false-true duplicates" are however useful for us to cluster the potential BBs from the same Amplicon DNA together. Though they belonging to different BB, they still have excatly the same sequence, which means they are more likely come from the same DNA. By using the distance matrix weighted by the number and length of reads 'duplicated' between each two BBs, cluters can be defined and used fro assembly. With less cluters to get, more double-BBs-cluters found (fig.4). While the reads distribution a relativly stable among different levels to cluter (fig.5).

![fig.4](04.test_vsearch/BCdistofBBs.png)

> Fig. 4

![fig.5](04.test_vsearch/BCdistofReads.png)

> fig.5

```bash
#CODE NOTE:
hierarchy read.merge.derep.2T.bc.graph.tree
Number of levels: 5
level 0: 20403 nodes
level 1: 4805 nodes   # I chose this one for test
level 2: 2255 nodes
level 3: 1889 nodes
level 4: 1867 nodes
```

Finally, we tested **4805** Beads Barcode Cluters (BBCs or BCs) identified by level 1 and processed assembly with reads belonging to each cluster.

#### Annotation summary

Assemble were processed by Spades with `--careful` and `--cov-cutoff` to avoid the mismatchs.

![fig.6](04.test_vsearch/annotationStat.png)

> **Table 1. Assemle and annotation statistics obtained by randomly picked BCs with beads number from least to highest.**

#### BC annotation consistency

Among **4805** BCs, **2432** of them can be assemble to at least one scaffolds. Those scaffolds were then aligned to UNITE database. **1564** BCs with belonging scaffolds got alignment of  **length >= 50bp and identity >= 97%** were kept as their identification.

![fig.7](04.test_vsearch/AllAnnoStat.png)

> fig. 6 Annotation stat

Most of the BCs obtained unified indentification with same taxonomy.

| Diverse | species | genus | family | order | class | phylum | kindom |
| ------- | ------- | ----- | ------ | ----- | ----- | ------ | ------ |
| 1       | 1076    | 1095  | 1182   | 1373  | 1406  | 1408   | 1564   |
| 2       | 342     | 329   | 265    | 166   | 144   | 147    | 0      |
| 3       | 89      | 95    | 89     | 23    | 13    | 9      | 0      |
| 4       | 41      | 34    | 21     | 2     | 1     | 0      | 0      |
| 5       | 9       | 5     | 5      | 0     | 0     | 0      | 0      |
| 6       | 3       | 5     | 2      | 0     | 0     | 0      | 0      |
| 7       | 4       | 1     | 0      | 0     | 0     | 0      | 0      |

| Diverse% | species | genus  | family | order  | class  | phylum | kindom  |
| -------- | ------- | ------ | ------ | ------ | ------ | ------ | ------- |
| 1        | 68.80%  | 70.01% | 75.58% | 87.79% | 89.90% | 90.03% | 100.00% |
| 2        | 21.87%  | 21.04% | 16.94% | 10.61% | 9.21%  | 9.40%  | 0.00%   |
| 3        | 5.69%   | 6.07%  | 5.69%  | 1.47%  | 0.83%  | 0.58%  | 0.00%   |
| 4        | 2.62%   | 2.17%  | 1.34%  | 0.13%  | 0.06%  | 0.00%  | 0.00%   |
| 5        | 0.58%   | 0.32%  | 0.32%  | 0.00%  | 0.00%  | 0.00%  | 0.00%   |
| 6        | 0.19%   | 0.32%  | 0.13%  | 0.00%  | 0.00%  | 0.00%  | 0.00%   |
| 7        | 0.26%   | 0.06%  | 0.00%  | 0.00%  | 0.00%  | 0.00%  | 0.00%   |

> **Table 2 Annotation diversity in each BC** （counts and ratio）

Maybe due to the converved region between our 5 fungies, **20%**  BCs are seems combined with 2 species. Less than 1% are identified more than 4 species ,which is pesudo identification since we only mixed 5 fungi (one of which at very low concentration).

#### Case of BC

##### The simple case

**BC00005** , which contained only **3** BCs and **740** reads:

The mapping relationship from reads to genmoe reference shows below:

![BC00005.all](circos/BC00005.circos.i97.x10.png)

> **Fig.S1.a  The alignment relationshiop between reads to scafflods and scaffolds to genome references**.
> Black segments in the outerside are assembled scaffolds in the BC. Colored segments are genome reference ( x0.1 scaleed).

A simple BC inludes fewer BBs and reads, they can be easliy determined, though the longet  caffolds are  only 345.

**BC00935** , which contained only **47** BCs and **5074** reads:

The mapping relationship from reads to genmoe reference shows below:

![BC00935.all](circos/BC00935.circos.i97.x10.png)

> **Fig.S1.a  The alignment relationshiop between reads to scafflods and scaffolds to genome references**.
> Black segments in the outerside are assembled scaffolds in the BC. Colored segments are genome reference ( x0.1 scaleed). Black links are paried reads with each read aligned to different scaffolds. While colored links are alignemnt region between scaffolds to genome references. reads alignments are tile inside the scaffolds. The same BC are aligned at the same height and colored by BC IDs.

We can see 16 scaffolds are mainly come from Pe and Le.  2 of them are uniquly come from Sc. Though #14 can also aligned to Fv, this segement is too concerved to comfrim. The black linkes could help connect scalffods with gaps. Without those gaps, the target alignement position is posible to extended to as long as 2.5kb.

![BC00935.zoom](circos/BC00935.zoom.jpg)

> **Fig.S1 b. Zoomed scale of BC00005.** Reads with 100% alignemnt identification are drawn.

A closer look of the scaffolds #1 ~ #4.  #2 can both aligned to Pe and Lv, while the others are uniquely aligned to Lv. Since #1 ~ #4 share most of the BBs, which means they should be physically come from the sam DNA source, we can surly confirm that #2 are also all come from Lv.

##### The most complexed case

**BC00491** , which contained **89** BCs and **15769** reads:

The mapping relationship from reads to genmoe reference shows below:

![BC00491.all](circos/BC00491.circos.i97.x10.png)

> **Fig.S1.a  The alignment relationshiop between reads to scafflods and scaffolds to genome references**.
> Black segments in the outerside are assembled scaffolds in the BC. Colored segments are genome reference ( x0.1 scaleed). Black links are paried reads with each read aligned to different scaffolds. While colored links are alignemnt region between scaffolds to genome references. reads alignments are tile inside the scaffolds. The same BC are aligned at the same height and colored by BC IDs.

With the identification >= 97, assembled  scafffolds are mainly mapped to all of Pe, Hm, Le, Le, Fv, and Sc at least one chromosome (Fig.S1.a).  This indcated that the BC with too much BBs are likely combiend by different species sources.

![BC00491.zoom](circos/BC00491.zoom.jpg)

> **Fig.S1 b. Zoomed scale of BC00491.** Reads with 100% alignemnt identification are drawn.

The scaffolds #1 and #6 are assembled by reads comming from variouse BBs, while #1 could mapping both to Pe and Le. While #6 all aligned to Pe. So we can say, those reads with specific BBs assembled to be #6 should be the one from Pe. Then the other reads belonging to the same BBs should also come from Pe, thouth their comritbuted scaffolds might be aligned to other species due the potential conserved sequences.

This strategy can help us split BBs in a single BC to their  taxonomy source. Then refined assemble could be possible to make the scaffolds more perfect.
