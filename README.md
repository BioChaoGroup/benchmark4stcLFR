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



#### Module

**stLFR barcodes split** - [README.md](./)

