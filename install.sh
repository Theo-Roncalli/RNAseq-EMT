#!/usr/bin/env bash

# Directory parameters

reads=Data/Reads
trimming=Data/Trimming
genome=Data/Genome
figures_reads=Figures/Reads
figures_trimming=Figures/Trimming

# Url parameters

reads_url=http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz
genome_url=http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz
annotation_url=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz

# Color parameters

RED='\033[1;31m'
GREEN='\033[1;32m'
BLUE='\033[1;34m'
NC='\033[0m'

# Step 1: Download reads

mkdir -p ${reads}
echo -e "\n${BLUE}Downloading reads...${NC}"
wget ${reads_url} -P ${reads}
echo -e "${GREEN}Done.${NC}"
echo -e "\n${BLUE}Unarchiving reads...${NC}"
tar -zxf ${reads}/TPrnaseq.tar.gz -C ${reads}
echo -e "${GREEN}Done.${NC}"


echo -e "\n${BLUE}---------------Number of reads per file---------------${NC}"
for file in ${reads}/*.fastq
do
	grep ^+$ ${file} | echo "${file}: $(wc -l) reads";
done

# Step 2: Quality control + Reads cleaning

mkdir -p ${figures_reads}
echo -e "\n${BLUE}Creation of the fastqc files on raw reads...${NC}"
fastqc -o ${figures_reads} -f fastq ${reads}/*.fastq -q
echo -e "${GREEN}Done.${NC}"

# Trimming procedure (Elimination of low quality sequences at the end of reads)
# conda install -c bioconda trimmomatic

mkdir -p ${trimming}

for read1_file in ${reads}/*.R1.fastq
do
	paired_file_with_path=${read1_file%.R1.fastq};
	paired_file_without_path=${paired_file_with_path#${reads}/};
	echo -e "\n${BLUE}Trimming ${paired_file_without_path%.sampled}...${NC}";
	trimmomatic PE ${paired_file_with_path}.R1.fastq ${paired_file_with_path}.R2.fastq -baseout ${trimming}/${paired_file_without_path}.fastq LEADING:20 TRAILING:20 MINLEN:50 -quiet
echo -e "${GREEN}Done.${NC}"
done

# Removal of the bases from the extremity with a quality lower than 20. If the final read is smaller than 50, it is discarded. file with U => discard. file with P => no discard.
# Remark: files 1U and 2U returns a small number of sequences (around 10 000) while files 1P and 2P returns a large number of sequences (a little smaller than the reads without cleaning)

mkdir -p ${figures_trimming}
echo -e "\n${BLUE}Creation of the fastqc files on trimmed reads...${NC}"
fastqc -o ${figures_trimming} -f fastq ${trimming}/*.fastq -q
echo -e "${GREEN}Done.${NC}"

# Step 2: download reference genome

mkdir ${genome} -p
echo -e "\n${BLUE}Downloading genome...${NC}"
wget ${genome_url} -P ${genome} -q
echo -e "${GREEN}Done.${NC}"
echo -e "\n${BLUE}Downloading annotations...${NC}"
wget ${annotation_url} -P ${genome} -q
echo -e "${GREEN}Done.${NC}"
gunzip ${genome}/*.gz
