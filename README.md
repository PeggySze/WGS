This repository contains code to call germline SNPs and indels using [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

<img src="https://github.com/PeggySze/WGS/blob/main/GATK_pipeline.png" style="display: block; margin: auto;" />

### Step1: Quality control and trimming of raw reads by fastp
- scirpt: 01_fastp.sh

### Step2: Mapping reads to reference genome by BWA mem
- script: 02_read_mapping.sh

### Step3: Examination of mapping quality 
- metrics:
	- the percentage of mapped reads
	- avergae sequencing coverage 
- scirpt: 03_mapping_qc.sh

### Step4: Call SNPs by GATK HaplotypeCaller
- scirpt: 04_HaplotypeCaller.sh

### Step5: Filter variants
- scirpt: 05_filter_reads.sh

### Step6: Annotate filtered variants
- scirpt: 06_annotate_filtered_variants.sh

