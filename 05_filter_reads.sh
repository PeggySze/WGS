#!/bin/bash

## PBS configure
#PBS -N filter_variants_and_counts
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
module load gatk/4.1.4.0
module load vcftools/v0.1.16

## Set variables
input_dir=~/chinmo_project/input/fastq/
vcf_dir=~/chinmo_project/output/haplotypecaller/
unfiltered_vcf_dir=~/chinmo_project/output/haplotypecaller/unfiltered_vcfs/
filtered_vcf_dir=~/chinmo_project/output/haplotypecaller/filtered_vcfs/
filtered_var_counts_dir=~/chinmo_project/output/haplotypecaller/filtered_vcfs_metrics/


## Mark and remove variants below quality thresholds
for i in ${input_dir}sampleID.txt;
do
	# Select SNPs from VCF files
    vcftools --vcf ${vcf_dir}${i}.HC.vcf --remove-indels --recode --recode-INFO-all \
    --chr X --chr 2L --chr 2R --chr 3L --chr 3R --chr Y --chr 4 \
    --stdout >${unfiltered_vcf_dir}${i}.SNPs.vcf

    # Mark SNPs below quality thresholds
    ~/software/gatk/4.1.4.0/gatk VariantFiltration\
    -V ${unfiltered_vcf_dir}${i}.SNPs.vcf  \
    --filter-expression "QD<2.00" --filter-name "QD" \
    --filter-expression "FS>60.00" --filter-name "FS" \
    --filter-expression "SOR>3.00" --filter-name "SOR" \
    --filter-expression "MQ<40.00" --filter-name "MQ" \
    --filter-expression "MQRankSum<-12.50" --filter-name "MQRankSum" \
    --filter-expression "ReadPosRankSum<-5.00" --filter-name "ReadPosRankSum" \
    -O ${filtered_vcf_dir}${i}.marked.SNPs.vcf

    # Remove SNPs below quality thresholds
    ~/software/gatk/4.1.4.0/gatk SelectVariants \
    -V ${filtered_vcf_dir}${i}.marked.SNPs.vcf \
    --exclude-filtered \
    -O ${filtered_vcf_dir}${i}.filtered.SNPs.vcf

    # Select indels from VCF files
    vcftools --vcf ${vcf_dir}${i}.HC.vcf --keep-only-indels --recode --recode-INFO-all \
    --chr X --chr 2L --chr 2R --chr 3L --chr 3R --chr Y --chr 4 \
    --stdout >${unfiltered_vcf_dir}${i}.indels.vcf

    # Mark indels below quality thresholds
    ~/software/gatk/4.1.4.0/gatk VariantFiltration\
    -V ${unfiltered_vcf_dir}${i}.indels.vcf  \
    --filter-expression "QD<2.00" --filter-name "QD" \
    --filter-expression "FS>60.00" --filter-name "FS" \
    --filter-expression "SOR>3.00" --filter-name "SOR" \
    --filter-expression "MQ<40.00" --filter-name "MQ" \
    --filter-expression "MQRankSum<-12.50" --filter-name "MQRankSum" \
    --filter-expression "ReadPosRankSum<-8.00" --filter-name "ReadPosRankSum" \
    -O ${filtered_vcf_dir}${i}.marked.indels.vcf

    # Remove indels below quality thresholds
    ~/software/gatk/4.1.4.0/gatk SelectVariants \
    -V ${filtered_vcf_dir}${i}.marked.indels.vcf  \
    --exclude-filtered \
    -O ${filtered_vcf_dir}${i}.filtered.indels.vcf

    # Merge filtered SNPs and indels
    ~/software/gatk/4.1.4.0/gatk MergeVcfs \
    -I ${filtered_vcf_dir}${i}.filtered.SNPs.vcf \
    -I ${filtered_vcf_dir}${i}.filtered.indels.vcf \
    -O ${filtered_vcf_dir}${i}.filtered.vcf
done

## Split filtered biallelic and multi SNPs/Indels
for i in ${input_dir}sampleID.txt;
do
    # Generate filtered_bi_snps.vcf
    vcftools --vcf ${filtered_vcf_dir}${i}.filtered.SNPs.vcf --recode --recode-INFO-all \
    --min-alleles 2 --max-alleles 2 --stdout >${filtered_vcf_dir}${i}_filtered_bi_SNPs.vcf

    # Generate filtered_multi_snps.vcf
    vcftools --vcf ${filtered_vcf_dir}${i}.filtered.SNPs.vcf --recode --recode-INFO-all \
    --min-alleles 3 --stdout >${filtered_vcf_dir}${i}_filtered_multi_SNPs.vcf

    # Generate filtered_bi_indels.vcf
    vcftools --vcf ${filtered_vcf_dir}${i}.filtered.indels.vcf --recode --recode-INFO-all \
    --min-alleles 2 --max-alleles 2 --stdout >${filtered_vcf_dir}${i}_filtered_bi_indels.vcf

    # Generate filtered_multi_indels.vcf
    vcftools --vcf ${filtered_vcf_dir}${i}.filtered.indels.vcf --recode --recode-INFO-all \
    --min-alleles 3 --stdout >${filtered_vcf_dir}${i}_filtered_multi_indels.vcf ;
done

## Count the number of filtered SNPs and Indels
for i in ${input_dir}sampleID.txt;
do
    # Count total filtered variations
    grep -c "^[^#]" ${filtered_vcf_dir}${i}.filtered.vcf  >>${filtered_var_counts_dir}total_filtered_var_counts

    # Count the number of filtered biallelic SNPs
    grep -c "^[^#]" ${filtered_vcf_dir}${i}_filtered_bi_SNPs.vcf >>${filtered_var_counts_dir}filtered_bi_SNPs_co

    # Count the number of filtered multiallelic SNPs
    grep -c "^[^#]" ${filtered_vcf_dir}${i}_filtered_multi_SNPs.vcf >>${filtered_var_counts_dir}filtered_multi_s

    # Count the number of filtered biallelic indels
    grep -c "^[^#]" ${filtered_vcf_dir}${i}_filtered_bi_indels.vcf >>${filtered_var_counts_dir}filtered_bi_indel

    # Count the number of filtered multiallelic indels
    grep -c "^[^#]" ${filtered_vcf_dir}${i}_filtered_multi_indels.vcf >>${filtered_var_counts_dir}filtered_multi
done

# Unload tools
module unload gatk/4.1.4.0
module unload vcftools/v0.1.16

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration

