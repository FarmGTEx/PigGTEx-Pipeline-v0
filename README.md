# Analysis pipelines for the pilot phase (V0) of PigGTEx 

## 1. Introduction

As an essential part of the ongoing Farm animal Genotype-Tissue Expression (**FarmGTEx**, https://www.farmgtex.org/) project, the **PigGTEx** (http://piggtex.farmgtex.org/) aims to build a comprehensive public catalogue of genetic regulatory variants in distinct biological contexts (e.g., tissues, cells, development, sex, environmental exposures, and genetic backgrounds) in pigs. The PigGTEx resources will be updated every 2-3 years together with Functional Annotation of Animal Genomes (FAANG) project (FAANG/FarmGTEx-TF). The PigGTEx will serve as a valuable resource for basic biological discovery, pig genetics, breeding, biotechnology, domestication, veterinary and human biomedicine research.



![PigGTEx-pipeline-V1.1](https://raw.githubusercontent.com/FarmGTEx/PigGTEx-Pipeline-v0/master/PigGTEx-pipeline-V1.1.jpg)



## 2. Analysis pipelines

This repository contains analysis pipelines used by the PigGTEx Consortium, including:

- Whole-genome sequence alignment, SNP calling, and quality control
- RNA-seq alignment, quantification, SNP calling, genotype imputation, and quality control
- Tissue-specific gene expression and co-expression analysis
- Bioinformatics analysis of WGBS data
- Single-cell RNA-Seq analysis
- Breed prediction for RNA-Seq samples
- Cis-heritability estimation for gene expression
- Molecular QTL (molQTL) mapping, meta-analysis, effect size measurement, and evolutionary conservation
- Allele-specific expression (ASE) analysis
- Functional annotation (snpEff and chromatin state) of molQTL
- Interpreting GWAS loci with molQTL: enrichment, expression-mediated heritability, transcriptome-wide association study, Mendelian Randomization, and co-localization
- Comparative analysis between pig and human



Pipeline components are public available in this repository and execution scripts are provided in snakemake workflow (.smk). 

## 3. Citation
#### **A compendium of genetic regulatory effects across pig tissues**

[https://www.nature.com/articles/s41588-023-01585-7](https://www.nature.com/articles/s41588-023-01585-7)

Jinyan Teng†, Yahui Gao†, Hongwei Yin†, Zhonghao Bai†, Shuli Liu†, Haonan Zeng†, The PigGTEx Consortium, Lijing Bai, Zexi Cai, Bingru Zhao, Xiujin Li, Zhiting Xu, Qing Lin, Zhangyuan Pan, Wenjing Yang, Xiaoshan Yu, Dailu Guan, Yali Hou, Brittney N. Keel, Gary A. Rohrer, Amanda K. Lindholm-Perry, William T. Oliver, Maria Ballester, Daniel Crespo-Piazuelo, Raquel Quintanilla, Oriol Canela-Xandri, Konrad Rawlik, Charley Xia, Yuelin Yao, Qianyi Zhao, Wenye Yao, Liu Yang, Houcheng Li, Huicong Zhang, Wang Liao, Tianshuo Chen, Peter Karlskov-Mortensen, Merete Fredholm, Marcel Amills, Alex Clop, Elisabetta Giuffra, Jun Wu, Xiaodian Cai, Shuqi Diao, Xiangchun Pan, Chen Wei, Jinghui Li, Hao Cheng, Sheng Wang, Guosheng Su, Goutam Sahana, Mogens Sandø Lund, Jack C.M. Dekkers, Luke Kramer, Christopher K. Tuggle, Ryan Corbett, Martien A.M. Groenen, Ole Madsen, Marta Gòdia, Dominique Rocha, Mathieu Charles, Cong-jun Li, Hubert Pausch, Xiaoxiang Hu, Laurent Frantz, Yonglun Luo, Lin Lin, Zhongyin Zhou, Zhe Zhang, Zitao Chen, Leilei Cui, Ruidong Xiang, Xia Shen, Pinghua Li, Ruihua Huang, Guoqing Tang, Mingzhou Li, Yunxiang Zhao, Guoqiang Yi, Zhonglin Tang, Jicai Jiang, Fuping Zhao, Xiaolong Yuan, Xiaohong Liu, Yaosheng Chen, Xuewen Xu, Shuhong Zhao, Pengju Zhao, Chris Haley, Huaijun Zhou, Qishan Wang, Yuchun Pan, Xiangdong Ding, Li Ma, Jiaqi Li, Pau Navarro, Qin Zhang, Bingjie Li, Albert Tenesa*, Kui Li*, George E. Liu*, Zhe Zhang*, Lingzhao Fang*

Nature Genetics, 2024, DOI: [10.1038/s41588-023-01585-7](https://doi.org/10.1038/s41588-023-01585-7)
