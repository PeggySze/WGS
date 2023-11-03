#!/bin/bash

## PBS configure
#PBS -N HaplotypeCaller
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
module load picard/2.21.1
module load gatk/4.1.4.0

## Set variables
input_dir=~/chinmo_project/input/fastq/
refgen=~/chinmo_project/input/fasta/dm6.fa
bam_dir=~/chinmo_project/output/bam/
vcf_dir=~/chinmo_project/output/haplotypecaller/

## Generate reference genome dictionary
java -jar ~/software/picard/2.21.1/picard.jar CreateSequenceDictionary \
        REFERENCE=${refgen} \
        OUTPUT=~/chinmo_project/input/fasta/dm6.dict


for i in ${input_dir}sampleID.txt;
do
	## Call variants per sample
	~/software/gatk/4.1.4.0/gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' HaplotypeCaller \
    -R ${refgen} \
    -I ${bam_dir}${i}_sorted_markup.bam \
    -O ${vcf_dir}${i}.g.vcf \
    -ERC GVCF;

    ## Perform joint genotyping on a singular sample by GenotypeGVCFs
    ~/software/gatk/4.1.4.0/gatk --java-options "-Xmx4g" GenotypeGVCFs\
    -R ${refgen} \
    -V ${vcf_dir}${i}.g.vcf \
    -O ${vcf_dir}${i}.HC.vcf 
done

# unload tools
module unload gatk/4.1.4.0
module unload picard/2.21.1


END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration



