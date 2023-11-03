#!/bin/bash

##PBS configure
#PBS -N fastp
#PBS -j oe
#PSB -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=4

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
uname -sa  ## Information about the operating system
set -ex  ## Log everything,and quit if there is any error
START=$(date +%s.%N)

# set up directories
input_dir=~/chinmo_project/input/fastq/
output_dir=~/chinmo_project/output/fastp/

## Quality control and read trimming by fastp
for samplename in ${input_dir}sampleID.txt
do
	/md01/shipy3/software/fastp -i ${input_dir}${samplename}_R1.fastq.gz \
	-I ${input_dir}${samplename}_R2.fastq.gz \
	-o ${output_dir}trimmed_${samplename}_R1.fastq.gz \
    -O ${output_dir}trimmed_${samplename}_R2.fastq.gz \
    -h ${output_dir}${samplename}_fastp.html \
    -j ${output_dir}${samplename}_fastp.json
done

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration






