#!/bin/bash

##PBS configure
#PBS -N read_mapping
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
uname -sa  ## Information about the operating system
set -ex  ## Log everything,and quit if there is any error
START=$(date +%s.%N)

## Source module environment and load tools
export MODULEPATH=$MODULEPATH:/public/home/shipy3/modulefiles
source /etc/profile.d/modules.sh
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.21.1

## Set up directories
input_dir=~/chinmo_project/input/fastq/
fqpath=~/chinmo_project/input/output/fastp/
refgen=~/chinmo_project/input/fasta/UCSC/dm6.fa
outputpath=~/chinmo_project/output/bam/

for i in ${input_dir}sampleID.txt;
do
	# Map reads to reference genome,using BWA mem
	bwa mem -t 5 -R "@RG\tID:L01\tSM:${i}\tPL:UNKNOWN\tLB:${i}" \
	${refgen} ${fqpath}trimmed_${i}_R1.fastq.gz ${fqpath}trimmed_${i}_R2.fastq.gz >${fqpath}${i}.sam ;

	# Compress sam files to bam files
	samtools view -bS ${fqpath}${i}.sam > ${fqpath}${i}.bam;
	rm ${fqpath}${i}.sam 

	# Sort BAM files
	samtools sort -O bam -o ${outputpath}${i}_sorted.bam ${fqpath}${i}.bam
	rm ${fqpath}${i}.bam

	# Mark duplicates by using picard
	java -jar picard.jar MarkDuplicates \
    I=${outputpath}${i}_sorted.bam \
    O=${outputpath}${i}_sorted_markup.bam \
    M=${outputpath}${i}_sorted_markup_metrics.txt

    rm ${outputpath}${i}_sorted.bam

    # Index the bam file
    samtools index ${outputpath}${i}_sorted_markup.bam
done

# unload tools
module unload bwa/0.7.17
module unload samtools/1.9
module unload picard/2.21.1

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration


