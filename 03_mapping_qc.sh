#!/bin/bash

## PBS configure
#PBS -N mapping_qc
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=10

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
uname -sa  ## Information about the operating system
set -ex  ## Log everything,and quit if there is any error
START=$(date +%s.%N)

## Source module environment and load tools
export MODULEPATH=$MODULEPATH:/public/home/shipy3/modulefiles
source /etc/profile.d/modules.sh
module load picard/2.21.1

## Set variables
input_dir=~/chinmo_project/input/fastq/
refgen=~/chinmo_project/input/fasta/UCSC/dm6.fa
bam_dir=~/chinmo_project/output/bam/
output_dir=~/chinmo_project/output/bam_qc/

## Collect metrics about coverage and performance of WGS by using picard
for i in ${input_dir}sampleID.txt;
do
        java -jar ~/software/picard/2.21.1/picard.jar CollectWgsMetrics \
        I=${bam_dir}${i}_sorted_markup.bam \
        O=${output_dir}${i}_collect_wgs_metrics.txt \
        R=${refgen} \
        INCLUDE_BQ_HISTOGRAM=true ;
done

# unload tools
module unload picard/2.21.1

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration
