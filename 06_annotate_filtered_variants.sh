#!/bin/bash

## PBS configure
#PBS -N annotate_filtered_variants
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=20

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
uname -sa  ## Information about the operating system
set -ex  ## Log everything,and quit if there is any error
START=$(date +%s.%N)

## Source module envrionment and load tools
export MODULEPATH=$MODULEPATH:/public/home/shipy3/modulefiles
source /etc/profile.d/modules.sh
module load snpeff/4.3

## Set variables
input_dir=~/chinmo_project/input/fastq/
filtered_vcfs_dir=~/chinmo_project/output/haplotypecaller/filtered_vcfs/

## Add annotations to filtered variations
for i in ${input_dir}sampleID.txt;
do
	# Annotate variations by snpEff
    java -Xmx4g -jar ~/software/snpeff/snpEff/snpEff.jar BDGP6.86 \
        -v -strict ${filtered_vcfs_dir}${i}.filtered.vcf >${filtered_vcfs_dir}${i}.filtered.ann.vcf
done

# Unload tools
module unload snpeff/4.3

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration

