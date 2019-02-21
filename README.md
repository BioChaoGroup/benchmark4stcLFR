# MetaSeq Stages Testing Project

## Why

1. `metaSeq`流程开发涉及多个独立开发的项目，无法在一个项目内记录完全；
2. 测试过程涉及真实样本数据和路径，不便过早公开；
3. 测试环节存在较多版本、参数迭代的过程，过于繁琐，不便于在流程本身的项目中进行。

## How

本项目会对流程中的每个环节进行测试，测试内容包括：

1. 平行整合能够实现该过程的、来自不同语言或不同项目的工具；
2. 对该过程中涉及到的参数进行评估和优化；
3. 尝试对该过程的输出、指标进行标准化；

## What

### Requirment

**Environment**: `python 3.6` `perl 5.5` 

**Developing projects**：[metaSeq](https://github.com/ZeweiSong/metaSeq) ([biogit](https://biogit.cn/Fangchao/metaSeq))、[cOMG](https://biogit.cn/Fangchao/Omics_pipeline)、[fastp](https://github.com/OpenGene/fastp) ([biogit](https://biogit.cn/PUB/fastp))

**Third party program:**

- **Snakemake** - a pythonic workflow system ([bitbucket](https://bitbucket.org/snakemake/snakemake))

### Init this repo

```shell
export D_fastp="/path/to/fastp/reporsitory"
export D_metaSeq="/path/to/fastp/repository"
git clone git@biogit.cn:Fangchao/benchmark4stlfr.git
cd benchmark4stlfr
# Find and modifiy the DATA and DB directory according to your repo location.
ln -s ../../../DATA
ln -s ../../../DB
```



### Test sample

All raw data saved under: /hwfssz1/ST_META/P18Z10200N0059_stLFR_SZW/DATA/
Data is organized as the date sequenced (e.g. test20180301)

`/hwfssz1/ST_META/P18Z10200N0059_stLFR_SZW/DATA/`

- test20180301/   # First round stLFR test on yeast (1.9k and 2.5k amplicon)
  CL100063954_L01/        # 2.5 kb amplicon
  CL100063954_L02/        # 1.9 kb amplicon

- test20180410/   # Second round stLFR test using five or four fungi mixture (2.5k amplicon)

  > Saccharomyces cerevisiae - Sc (Baker's yeast)
  > Pleurotus eryngii - Pe (Xing Bao Gu)
  > Lentinula edodes - Le (Xiang Gu, Shiitake)
  > Flammulina velutipes - Fv (Jin Zhen Gu)
  > Hypsizygus marmoreus - Hm (Xie Wei Gu, this fungus has very low PCR efficiency)

**Note**: 

1. Even mixture is the equal mix of **PCR product** of the five fungi (but **Hm has very low signals** in the data)
2. Stagger mixture is mixed with **1000Sc : 100Pe : 10Fv : 1Le**

| Tag    | Chip_lane ID     | content                       |
| ------ | ---------------- | ----------------------------- |
| APR831 | CL100076083_L01/ | Even mixture, 10 pg load      |
| APR832 | CL100076083_L02/ | Staggered mixture, 20 pg load |
| APR841 | CL100076084_L01/ | Staggered mixture, 10 pg load |
| APR842 | CL100076084_L02/ | Even mixture, 20 pg load      |

- test20180904/   # Third round stLFR test on adding umi adaptors using Le

**Note**:

:warning: These are NOT raw data but split data by Wang Ou
⚠️  Sequences are labelded as @CL200051307L2C001R001_14#715_449_646/1
⚠️ "#" instead of "/" was used as the first delimiter, our split script only use "/";

```
q1/     # Copy number of PCR product and umi is 1:1, this data is converted to JSON format
q5/     # Copy number of PCR product and umi is 1:5
```

#### Module

**stLFR barcodes split** - [README.md](./)

#### Technical Note

[REAMDME.md](./Assay/README.md)

