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

# Performance parameters

nb_cpus_indexing=7
nb_cpus_mapping=7
nb_cpus_counting=7

# Step 1: Quality control + Reads cleaning

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

# Step 2: Mapping

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

# Step 3: Index BAM (samtools index)

for bam_file in ${mapping}/*.bam
do
	bam_file_without_path=${bam_file#${mapping}/};
	echo "Indexing BAM file with ${bam_file_without_path}";
	samtools index ${bam_file};
	samtools stats ${bam_file};
done

# For more information about the indexation of a BAM file, please type: samtools stats ${mapping}/BAM_file | less

for file in ${mapping}/*.bam
do
	bam_file_without_path=${bam_file#${mapping}/};
	samtools index ${bam_file};
done

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
featureCounts -p -T ${nb_cpus_counting} -t gene -g gene_id -s 0 -a ${genome}/*.gtf -o ${counts}/counts.txt ${mapping}/*.bam

# Create a file with pairs between ENCODE and HUGO identifiers
perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' \
	${genome}/*.gtf | sort | uniq > encode-to-hugo.tab

sort ${counts}/counts.txt > ${counts}/sort_counts.txt
join ${counts}/sort_counts.txt encode-to-hugo.tab | grep "chr18" > ${counts}/paired_counts.txt
rm encode-to-hugo.tab ${counts}/sort_counts.txt
