# Technical Note

### Landscape

![fig1. landscape of the pipeline](../Results/general/fig.pipeline.png)

**Fig.1** Landscape of pipeline using stLFR barcodes information

Figure 1 shows a pipeline of how we detecting barcodes and used it to improve the assembly results. The core process is to treat sequences in bead level and cluster them according to the similarity of reads into each other beads. The clustered subset is called Bead-barcodes Cluters. Sequence belongint to each cluter was assembled and annotated parallelled.

### Method

#### Sampling

For Fungus, we prepared a mock sample, which is evenly mixtured by 

- *Saccharomyces cerevisiae* - Sc (Baker's yeast), 
- *Pleurotus eryngii* - Pe (Xing Bao Gu), 
- *Lentinula edodes* - Le (Xiang Gu, Shiitake), 
- *Flammulina velutipes* - Fv (Jin Zhen Gu), 
- *Hypsizygus marmoreus* - Hm (Xie Wei Gu, this fungus has very low PCR efficiency).  

For Bacteria, a soil sample collected from XXX were tested.

**Sequencing** was performed by BGISEQ500(?)

**Quality control** was implented in a modified **[fastp](https://github.com/Scelta/fastp)**, along with the function to detect bead barcodes.

**VSEARCH** was used to compute duplicates and reads pre-clusters.

Codes and pipelines are availiable by accessing the project **metaseq**.

### Results

