#!/usr/bin/env bash

output_fastq=Data/Fastq
output_trimming=Data/Trimming
output_figures_fastq=Figures/Fastq
output_figures_trimming=Figures/Trimming
output_genome=Data/Genome
output_index=Data/Index
output_counts=Data/Counts

nb_cpus_indexing=7
nb_cpus_mapping=7

# Récupération des fichiers fastq

mkdir -p ${output_fastq}
wget http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz -P ${output_fastq} -q
tar -zxvf ${output_fastq}/TPrnaseq.tar.gz -C ${output_fastq}

echo "--------Number of sequences per file--------"
for file in ${output_fastq}/*.fastq
do
	grep + $file | echo "$file  $(wc -l)";
done

# Contrôle qualité

mkdir -p ${output_figures_fastq}
fastqc -o ${output_figures_fastq} -f fastq ${output_fastq}/*.fastq -q
# Trimming procedure (Elimination of low quality sequences at the end of reads)
# conda install -c bioconda trimmomatic

mkdir -p ${output_trimming}

for first_read_file in ${output_fastq}/*.R1.fastq
do
	paired_file_with_path=${first_read_file%.R1.fastq};
	paired_file_without_path=${paired_file_with_path#${output_fastq}/};
	echo "Trimming ${paired_file_without_path%.sampled}...";
	trimmomatic PE ${paired_file_with_path}.R1.fastq ${paired_file_with_path}.R2.fastq -baseout ${output_trimming}/${paired_file_without_path}.fastq LEADING:20 TRAILING:20 MINLEN:50 -quiet
	echo "Done."
done

# Removal of the bases from the extremity with a quality lower than 20. If the final read is smaller than 50, it is discarded. file with U => discard. file with P => no discard.

mkdir -p ${output_figures_trimming}
fastqc -o ${output_figures_trimming} -f fastq ${output_trimming}/*.fastq -q

# 5/ Mapping

# 5.1 Genome and annotation recovery

mkdir ${output_genome} -p
echo "Downloading genome..."
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz -P ${output_genome} -q
echo "Downloading annotations..."
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz -P ${output_genome} -q
gunzip ${output_genome}/*.gz

# 5.2 Indexation with STAR

input_genome=${output_genome}

mkdir ${output_index} -p
STAR --runMode genomeGenerate --runThreadN ${nb_cpus_indexing} \
	--genomeSAindexNbases 12 \
	--genomeDir ${output_index} \
	--genomeFastaFiles ${input_genome}/chr18.fa \
	--sjdbGTFfile ${input_genome}/gencode.v24lift37.basic.annotation.gtf

# 5.3 Mapping with STAR

mkdir ${output_counts} -p
for first_read_file in ${output_trimming}/*1P.fastq
do
	paired_file_with_path=${first_read_file%_1P.fastq};
	echo "Downaloading BAM file with ${paired_file_with_path#${output_trimming}/}..."; 
	STAR --runThreadN ${nb_cpus_mapping} --outFilterMultimapNmax 1 \
	--genomeDir ${output_index} \
	--outSAMattributes All --outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${output_counts}/${paired_file_with_path}_ \
	--readFilesIn ${paired_file_with_path}_1P.fastq ${paired_file_with_path}_2P.fastq;
done

# ATTENTION : on ne peut pas utiliser les trimming qui ont supprimé les petites séquences. En effet, les deux fichiers 1U et 2U n'ont pas le même nombre de lignes. Comment pourrait-on régler ce problème ?
# Message affiché mais portant à confusion : EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length


