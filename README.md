# New Generation Sequencing (NGS)

This work focuses on the study of [Yang et al. (2016)](http://www.ncbi.nlm.nih.gov/pubmed/?term=27044866) who were interested in the epthithelium-mesenchymal transition (EMT) process. In their work, the EMT has been induced by ectopic expression of Zeb1 in a cell lung cancer cell line (H358). The authors have studied RNAseq data over 7 days, starting from uninduced cells.

The initial data are available on the [NCBI site](http://www.ncbi.nlm.nih.gov/sra?term=SRP066794). In order to reduce time computation, we used only 0.5% of the total RNAseq data at the following address: http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz

## Dependencies

The pipeline runs on bash.
Some package are required for launching some commands such as fastqc, trimmomatic and featureCounts.

```bash
sudo apt-get install -y fastqc # For using fastqc
sudo apt-get install -y subread # For using featureCounts
conda install -c bioconda trimmomatic # For using trimmomatic
```

## Hardware requirements

A machine with at least 16 GB of **FREE** RAM (to create the index and the mapping on the chromosome 18 of the reference genome).

## Executing The Pipeline

The pipeline is used to create a file named "hugo-counts.txt" to which is associated, for each gene, the HUGO identifier and the number of reads aligned for each observation. This file is available in the repository _Data/Counts_. The steps are the followings.

1. Clone the Github repository to your machine
```bash
git clone https://github.com/Theo-Roncalli/RNAseq-EMT.git
cd RNAseq-EMT
```

2. Importation of reads and reference genome
```bash
bash import.sh
```

3. Creation of the counting file which contains, for each HUGO code in Chromosome 18, the numbers of reads per gene and per observation.
```bash
bash counting.sh
```

## Cleaning Repository

For cleaning the repository (i.e. delete _Data_ and _Figures_ folders), please type:
```bash
bash clean.sh
```
