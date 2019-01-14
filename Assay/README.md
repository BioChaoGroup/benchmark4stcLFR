# Outline and Technical Note

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

For Bacteria, 