#!/usr/bin/env bash

# Directory parameters

reads=Data/Reads
trimming=Data/Trimming
genome=Data/Genome
index=Data/Index
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

# Récupération des fichiers fastq

mkdir -p ${reads}
wget ${reads_url} -P ${reads}
tar -zxvf ${reads}/TPrnaseq.tar.gz -C ${reads}

echo "--------Number of sequences per file--------"
for file in ${reads}/*.fastq
do
	grep ^+$ ${file} | echo "${file}  $(wc -l)";
done

# Quality control

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

mkdir -p ${figures_trimming}
fastqc -o ${figures_trimming} -f fastq ${trimming}/*.fastq -q

# 5/ Mapping

# 5.1 Genome and annotation recovery

mkdir ${genome} -p
echo "Downloading genome..."
wget ${genome_url} -P ${genome} -q
echo "Downloading annotations..."
wget ${annotation_url} -P ${genome} -q
gunzip ${genome}/*.gz

# 5.2 Indexation with STAR

input_genome=${genome}

mkdir ${index} -p
STAR --runMode genomeGenerate --runThreadN ${nb_cpus_indexing} \
	--genomeSAindexNbases 12 \
	--genomeDir ${index} \
	--genomeFastaFiles ${input_genome}/chr18.fa \
	--sjdbGTFfile ${input_genome}/gencode.v24lift37.basic.annotation.gtf

# 5.3 Mapping with STAR

mkdir ${counts} -p
for read1_file in ${trimming}/*1P.fastq
do
	paired_file_with_path=${read1_file%_1P.fastq};
	paired_file_without_path=${paired_file_with_path#${trimming}/};
	echo "Downaloading BAM file with ${paired_file_without_path}..."; 
	STAR --runThreadN ${nb_cpus_mapping} --outFilterMultimapNmax 1 \
	--genomeDir ${index} \
	--outSAMattributes All --outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${counts}/${paired_file_without_path}_ \
	--readFilesIn ${paired_file_with_path}_1P.fastq ${paired_file_with_path}_2P.fastq;
done

# ATTENTION : on ne peut pas utiliser les trimming qui ont supprimé les petites séquences. En effet, les deux fichiers 1U et 2U n'ont pas le même nombre de lignes. Comment pourrait-on régler ce problème ?
# Message affiché mais portant à confusion : EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length


