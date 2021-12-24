#!/usr/bin/env bash

# Directory parameters

reads=Data/Reads
genome=Data/Genome

# Url parameters

reads_url=http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz
genome_url=http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz
annotation_url=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz


# Step 1: Download reads

mkdir -p ${reads}
echo "Downloading reads..."
wget ${reads_url} -P ${reads}
tar -zxvf ${reads}/TPrnaseq.tar.gz -C ${reads}

echo "--------Number of sequences per file--------"
for file in ${reads}/*.fastq
do
	grep ^+$ ${file} | echo "${file}  $(wc -l)";
done

# Step 2: download reference genome

mkdir ${genome} -p
echo "Downloading genome..."
wget ${genome_url} -P ${genome} -q
echo "Downloading annotations..."
wget ${annotation_url} -P ${genome} -q
gunzip ${genome}/*.gz
