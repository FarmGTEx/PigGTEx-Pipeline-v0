# Analysis pipelines for the pilot phase (V0) of PigGTEx 

## 1. Introduction

As an essential part of Farm animal Genotype-Tissue Expression (**FarmGTEx**) project, the **PigGTEx** aims to build a comprehensive public catalogue of genetic regulatory variants in distinct biological contexts in pigs (e.g., tissues, cells, development, sex, environmental exposures and genetic background). The PigGTEx will serve as a valuable resource for basic biological discovery, animal breeding biotechnology, and human biomedicine research.



![PigGTEx-pipeline-V1.1](https://raw.githubusercontent.com/FarmGTEx/PigGTEx-Pipeline-v0/master/PigGTEx-pipeline-V1.1.svg)



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
