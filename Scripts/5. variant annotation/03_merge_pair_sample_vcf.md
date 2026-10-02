# merge the pair of samples' VCF files
* For each of the three most divergent farm pairs, the corresponding SNP VCF files from the two farms were merged
* resulting in a vcf file containing both overlapping and sample-specific SNPs




## SLURM bash script

### Create the bash script:

```bash
01_merge_vcf.sh

```

### Script content
```bash

#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 10:00:00
#SBATCH -p shared
#SBATCH -J meg_vcf
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se


module load bcftools/1.20


####### arguments
SAMPLE_ID_1=$1
SAMPLE_ID_2=$2

SAMPLE_ID_3=$3  # optional



####### input
input_base='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-03_normalise_vcf'

input_dir_1="$input_base/$SAMPLE_ID_1"

input_dir_2="$input_base/$SAMPLE_ID_2"

sample_snpbcf_1="$input_dir_1/${SAMPLE_ID_1}_filtered_snp.norm.pureSNP.bcf"

sample_snpbcf_2="$input_dir_2/${SAMPLE_ID_2}_filtered_snp.norm.pureSNP.bcf"


## if there is a 3rd sample

if [ -n "$SAMPLE_ID_3" ]; then   # if we have 3rd sample, -n = string length > 0
	input_dir_3="$input_base/$SAMPLE_ID_3"
	sample_snpbcf_3="$input_dir_3/${SAMPLE_ID_3}_filtered_snp.norm.pureSNP.bcf"
fi





#######  output
output_dir="/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-04_merge_vcf"

mkdir -p "$output_dir"

if [ -n "$SAMPLE_ID_3" ]; then
	output_merged_bcf="$output_dir/merged_${SAMPLE_ID_1}_${SAMPLE_ID_2}_${SAMPLE_ID_3}-filtered_snp.norm.pureSNP.bcf"
	output_merged_vcf="$output_dir/merged_${SAMPLE_ID_1}_${SAMPLE_ID_2}_${SAMPLE_ID_3}-filtered_snp.norm.pureSNP.vcf"
else
	output_merged_bcf="$output_dir/merged_${SAMPLE_ID_1}_${SAMPLE_ID_2}-filtered_snp.norm.pureSNP.bcf"
	output_merged_vcf="$output_dir/merged_${SAMPLE_ID_1}_${SAMPLE_ID_2}-filtered_snp.norm.pureSNP.vcf"
fi




####### merge a pair of samples vcf/bcf files of pure SNPs

if [ -n "$SAMPLE_ID_3" ]; then
	bcftools merge --threads 8 -Ob -m none $sample_snpbcf_1 $sample_snpbcf_2 $sample_snpbcf_3 > $output_merged_bcf

else
	bcftools merge --threads 8 -Ob -m none $sample_snpbcf_1 $sample_snpbcf_2 > $output_merged_bcf
fi



####### for reviewing and annotation convert BCF → VCF
bcftools view -Ov $output_merged_bcf -o $output_merged_vcf


## -O: Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v)

## -m none: merge method, keep all alternate alleles exactly as they are. Do not collapse or merge different ALT alleles. Safest option.

```


### Submit the SLURM job

run and provide the sample IDs you want to merge:
```bash
### merging a pair of farm samples: s15, s4 
# NOTE: the order of two samples name as arguments, should be same for all commands
sbatch 01_merge_vcf.sh Sample_XJ-4126-15 Sample_XJ-4126-4-ampure-1

### merging a pair of farm samples: s3, s4 
sbatch 01_merge_vcf.sh Sample_XJ-4126-3 Sample_XJ-4126-4-ampure-1


### merging a pair of farm samples: s4, G (replicates s7 + s9)
## as farm G has 2 replicates samples s7 and s9, so we just merged all three samples
sbatch 01_merge_vcf.sh Sample_XJ-4126-4-ampure-1 Sample_XJ-4126-7 Sample_XJ-4126-9  

```


## Running time:
~2 mins


## Input
one example:
```bash
# a pair of farms
025_snpEFF_vcf_annota/025-03_normalise_vcf/Sample_XJ-4126-15/

Sample_XJ-4126-15_filtered_snp.norm.pureSNP.bcf
Sample_XJ-4126-4-ampure-1_filtered_snp.norm.pureSNP.bcf

```


## Output
one example:

```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/
025_snpEFF_vcf_annota/025-04_merge_vcf> du -sh *

755M	merged_Sample_XJ-4126-15_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.bcf
4.0G	merged_Sample_XJ-4126-15_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.vcf

###### e.g. content of the merged vcf:
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample_XJ-4126-15       Sample_XJ-4126-4-ampure-1

```

