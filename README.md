---
title: "Analysis of human differentially expressed genes during the EMT process"
author: "Théo Roncalli"
date: "24 décembre 2021"
---

# New Generation Sequencing (NGS)

## Dependencies

The pipeline runs on [bash](https://www.nextflow.io/).
Some package are required for launching some commands such as fastqc, trimmomatic and featureCounts

```bash
sudo apt-get install -y fastqc # For using fastqc
sudo apt-get install -y subread # For using featureCounts
conda install -c bioconda trimmomatic # For using trimmomatic
```

## Hardware requirements

A machine with at least 16 GB of **FREE** RAM (to create the index and the mapping on the chromosome 18 of the reference genome).

## Executing The Pipeline

1. Clone the Github repository to your machine
```bash
git clone https://github.com/Theo-Roncalli/RNAseq-EMT.git
cd RNAseq-EMT
```

2. Importation of reads and reference genome
```bash
bash import.sh
```

3. Creation of the counting file which contains, for each HUGO codes in Chromosome 18, the numbers of reads per gene and per observation.
```bash
bash counting.sh
```
