#!/usr/bin/env bash

# Directory parameters

reads=Data/Reads
trimming=Data/Trimming
genome=Data/Genome
index=Data/Index
mapping=Data/Mapping
counts=Data/Counts
figures_reads=Figures/Reads
figures_trimming=Figures/Trimming

# Url parameters

reads_url=http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz
genome_url=http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz
annotation_url=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz

# Performance parameters

nb_cpus_indexing=7
nb_cpus_mapping=7

# Step 1: Download reads

mkdir -p ${reads}
wget ${reads_url} -P ${reads}
tar -zxvf ${reads}/TPrnaseq.tar.gz -C ${reads}

echo "--------Number of sequences per file--------"
for file in ${reads}/*.fastq
do
	grep ^+$ ${file} | echo "${file}  $(wc -l)";
done

# Step 2: Quality control + Reads cleaning

mkdir -p ${figures_reads}
fastqc -o ${figures_reads} -f fastq ${reads}/*.fastq -q

# Trimming procedure (Elimination of low quality sequences at the end of reads)
# conda install -c bioconda trimmomatic

mkdir -p ${trimming}

for read1_file in ${reads}/*.R1.fastq
do
	paired_file_with_path=${read1_file%.R1.fastq};
	paired_file_without_path=${paired_file_with_path#${reads}/};
	echo "Trimming ${paired_file_without_path%.sampled}...";
	trimmomatic PE ${paired_file_with_path}.R1.fastq ${paired_file_with_path}.R2.fastq -baseout ${trimming}/${paired_file_without_path}.fastq LEADING:20 TRAILING:20 MINLEN:50 -quiet
	echo "Done."
done

# Removal of the bases from the extremity with a quality lower than 20. If the final read is smaller than 50, it is discarded. file with U => discard. file with P => no discard.
# Remark: files 1U and 2U returns a small number of sequences (around 10 000) while files 1P and 2P returns a large number of sequences (a little smaller than the reads without cleaning)

mkdir -p ${figures_trimming}
fastqc -o ${figures_trimming} -f fastq ${trimming}/*.fastq -q

# Step 3: Mapping

#### Genome and annotation recovery ####

mkdir ${genome} -p
echo "Downloading genome..."
wget ${genome_url} -P ${genome} -q
echo "Downloading annotations..."
wget ${annotation_url} -P ${genome} -q
gunzip ${genome}/*.gz

#### Index (STAR) ####

mkdir ${index} -p
STAR --runMode genomeGenerate --runThreadN ${nb_cpus_indexing} \
	--genomeSAindexNbases 12 \
	--genomeDir ${index} \
	--genomeFastaFiles ${genome}/chr18.fa \
	--sjdbGTFfile ${genome}/gencode.v24lift37.basic.annotation.gtf

#### Mapping (STAR) ####

mkdir ${mapping} -p
for read1_file in ${trimming}/*1P.fastq
do
	paired_file_with_path=${read1_file%_1P.fastq};
	paired_file_without_path=${paired_file_with_path#${trimming}/};
	echo "Downloading BAM file with ${paired_file_without_path}..."; 
	STAR --runThreadN ${nb_cpus_mapping} --outFilterMultimapNmax 1 \
	--genomeDir ${index} \
	--outSAMattributes All --outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${mapping}/${paired_file_without_path}_ \
	--readFilesIn ${paired_file_with_path}_1P.fastq ${paired_file_with_path}_2P.fastq;
done

# Be careful : we cannot use the files 1U and 2U because the number of returned sequences is not the same.
# Consequently, we use the files 1P and 2P since the number of returned sequences is exactly the same.
# Indeed, the mapping is not achievable with the files 1U and 2U.
# Message when using 1U and 2U: EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length

# Step 4: Index BAM (samtools index)

samtools index ${mapping}.*bam

samtools stats ${mapping}.*bam



# samtools stats ${BAM} | less
# SN      reads MQ0:      0       # mapped and MQ=0
# SN      reads QC failed:        0
# On n'a pas de problèmes car on a élaguer avec trimmomatic.
# SN      inward oriented pairs:  273211
# SN      outward oriented pairs: 12
# On a un (petit) problème pour 12 reads. Il y a eu un réarrangement. En effet, 12 paires sont dans le même sens.

# BAM file indexation

for file in ${mapping}/*.bam
do
	paired_file_without_path=${paired_file_with_path#${trimming}/};
	echo "Indexation of ${file}"
	samtools index ${file}
done

# Creation of mapping matrix

mkdir -p ${counts}
featureCounts -p -T 7 -t gene -g gene_id -s 0 -a ${genome}/*.gtf -o ${counts}/counts.txt ${mapping}/*.bam

perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' \
	${genome}/*.gtf | sort | uniq > temp

sort counts.txt > temp1
sort temp > temp2
join temp1 temp2 |grep "chr18" > temp3


# featureCounts -p -T 7 -t gene -g gene_id -s 0 -a Homo_sapiens.GRCh38.99.gtf -o counts.txt SRR628582.bam SRR628584.bam SRR628585.bam SRR628583.bam SRR628587.bam SRR628589.bam SRR628586.bam SRR628588.bam
# featureCounts -p -T 7 -t gene -g gene_id -s 0 -a Data/Genome/gencode.v24lift37.basic.annotation.gtf -o counts.txt Data/Mapping/Day_0_1_chr18.sampled_Aligned.sortedByCoord.out.bam Data/Mapping/Day_0_2_chr18.sampled_Aligned.sortedByCoord.out.bam
