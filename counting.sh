#!/usr/bin/env bash

# Directory parameters

trimming=Data/Trimming
genome=Data/Genome
index=Data/Index
mapping=Data/Mapping
counts=Data/Counts

# Performance parameters

nb_cpus_indexing=7
nb_cpus_mapping=7
nb_cpus_counting=7

# Color parameters

RED='\033[1;31m'
GREEN='\033[1;32m'
BLUE='\033[1;34m'
NC='\033[0m'

# Step 1: Mapping

#### Index (STAR) ####

mkdir ${index} -p
echo -e "\n${BLUE}Creating the index on ${genome}/chr18.fa...${NC}"
STAR --runMode genomeGenerate --runThreadN ${nb_cpus_indexing} \
	--genomeSAindexNbases 12 \
	--genomeDir ${index} \
	--genomeFastaFiles ${genome}/chr18.fa \
	--sjdbGTFfile ${genome}/gencode.v24lift37.basic.annotation.gtf
echo -e "${GREEN}Done.${NC}\n"

#### Mapping (STAR) ####

mkdir ${mapping} -p
for read1_file in ${trimming}/*1P.fastq
do
	paired_file_with_path=${read1_file%_1P.fastq};
	paired_file_without_path=${paired_file_with_path#${trimming}/};
	echo -e "${BLUE}Downloading BAM file with ${paired_file_without_path}...${NC}";
	STAR --runThreadN ${nb_cpus_mapping} --outFilterMultimapNmax 1 \
	--genomeDir ${index} \
	--outSAMattributes All --outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${mapping}/${paired_file_without_path}_ \
	--readFilesIn ${paired_file_with_path}_1P.fastq ${paired_file_with_path}_2P.fastq;
    echo -e "${GREEN}Done.${NC}\n";
done

# Be careful : we cannot use the files 1U and 2U because the number of returned sequences is not the same.
# Consequently, we use the files 1P and 2P since the number of returned sequences is exactly the same.
# Indeed, the mapping is not achievable with the files 1U and 2U.
# Message when using 1U and 2U: EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length

# Step 3: Index BAM (samtools index)

for bam_file in ${mapping}/*.bam
do
	bam_file_without_path=${bam_file#${mapping}/};
    echo -e "${BLUE}Indexing ${bam_file_without_path}...${NC}";
	samtools index ${bam_file};
    echo -e "${GREEN}Done.${NC}\n";
done

# For more information about the indexation of a BAM file, please type: samtools stats ${mapping}/BAM_file | less

# samtools stats ${mapping}/Day_0_1_chr18.sampled_Aligned.sortedByCoord.out.bam | less
# SN      reads MQ0:      0       # mapped and MQ=0
# SN      reads QC failed:        0
# There is no issue since we have prune with trimmomatic.
# SN      inward oriented pairs:  314156
# SN      outward oriented pairs: 10
# There is a problem for 10 reads. Indeed, there have been a rearrangement since 10 pairs are in the same reading frame.

# Step 4: Counting (featureCounts)

# If not already install, please type: apt-get install subread

mkdir -p ${counts}
echo -e "${BLUE}Running featureCounts...${NC}"
featureCounts -p -T ${nb_cpus_counting} -t gene -g gene_id -s 0 -a ${genome}/*.gtf -o ${counts}/counts.txt ${mapping}/*.bam
echo -e "${GREEN}Done.${NC}\n"

# Create a file with pairs between ENCODE and HUGO identifiers
echo -e "${BLUE}Creating file with pairs between ENCODE and HUGO identifiers...${NC}";
perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' \
	${genome}/*.gtf | sort | uniq > ${counts}/encode-to-hugo.tab
echo -e "${GREEN}Done.${NC}\n"

echo -e "${BLUE}Sorting ${counts}/counts.txt file...${NC}"
sort ${counts}/counts.txt > ${counts}/sort_counts.txt
echo -e "${GREEN}Done.${NC}\n"

# Before joining the two files encode-to-hugo.tab and sort_counts.txt, we remove the lines which does not contain counting
# Indeed, the counts.txt file contains a header wich is the footer of the sort_counts.txt file.
# Indeed, there are these two lines:
# # Program:featureCounts v1.6.0; Command:"featureCounts" "-p" "-T" "7" "-t" "gene" "-g" "gene_id" "-s" "0" "-a" "Data/Genome/gencode.v24lift37.basic.annotation.gtf" "-o" "Data/Counts/counts.txt" "Data/Mapping/Day_0_1_chr18.sampled_Aligned.sortedByCoord.out.bam" "Data/Mapping/Day_0_2_chr18.sampled_Aligned.sortedByCoord.out.bam" "Data/Mapping/Day_0_3_chr18.sampled_Aligned.sortedByCoord.out.bam" "Data/Mapping/Day_7_1_chr18.sampled_Aligned.sortedByCoord.out.bam" "Data/Mapping/Day_7_2_chr18.sampled_Aligned.sortedByCoord.out.bam" "Data/Mapping/Day_7_3_chr18.sampled_Aligned.sortedByCoord.out.bam"
# Geneid	Chr	Start	End	Strand	Length	Data/Mapping/Day_0_1_chr18.sampled_Aligned.sortedByCoord.out.bam	Data/Mapping/Day_0_2_chr18.sampled_Aligned.sortedByCoord.out.bam	Data/Mapping/Day_0_3_chr18.sampled_Aligned.sortedByCoord.out.bam	Data/Mapping/Day_7_1_chr18.sampled_Aligned.sortedByCoord.out.bam	Data/Mapping/Day_7_2_chr18.sampled_Aligned.sortedByCoord.out.bam	Data/Mapping/Day_7_3_chr18.sampled_Aligned.sortedByCoord.out.bam
sed -i '/^[#|Geneid]/d' ${counts}/sort_counts.txt
# For verifying that the footer is effectively removed, please type: tail ${counts}/sort_counts.txt

# Creation of the hugo-counts.txt file which contains, for each HUGO code in Chromosome 18, the numbers of reads per gene and per observation.
if [ $(cat ${counts}/encode-to-hugo.tab | wc -l) == $(cat ${counts}/sort_counts.txt | wc -l) ]
then
	echo -e "${BLUE}Creation of the hugo-counts.txt file...${NC}"
	join ${counts}/encode-to-hugo.tab ${counts}/sort_counts.txt | grep "chr18" > ${counts}/paired_counts.txt
	awk '{print $2 " " $8 " " $9 " " $10 " " $11 " " $12 " " $13}' ${counts}/paired_counts.txt > ${counts}/hugo-counts.txt
	echo -e "${GREEN}Done.${NC}\n";
	echo -e "\n${GREEN}The final file to use is ${counts}/hugo-counts.txt."
	echo -e "It contains, for each HUGO codes in Chromosome 18, the numbers of reads per gene and per observation.${NC}\n"
	exit 0
else
    echo -e "\n${RED}The creation of a file containing the HUGO codes"
	echo -e "and numbers of reads per gene and per observation"
	echo -e "is not available since the number of genes is not the same"
	echo -e "between the encode-to-hugo.tab and sort_counts.txt files.${NC}\n"
	exit 1
fi
