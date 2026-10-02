# VCF Annotation with SnpEff


## 1. Annotate the merged SNP VCFs for the three most divergent pairs

```bash
01_snpEFF_anno_merged.sh

```

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -n 1
#SBATCH -c 8          # 16g/2g = 8 cores
#SBATCH -t 01:00:00
#SBATCH -p shared
#SBATCH -J anno_meg
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se


module load bioinfo-tools
module load snpEff


# arguments
SAMPLE_ID_1=$1
SAMPLE_ID_2=$2

SAMPLE_ID_3=$3 # optional

## input
input_dir="/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-04_merge_vcf"

if [ -n "$SAMPLE_ID_3" ]; then
	input_merged_vcf="$input_dir/merged_${SAMPLE_ID_1}_${SAMPLE_ID_2}_${SAMPLE_ID_3}-filtered_snp.norm.pureSNP.vcf"
else
	input_merged_vcf="$input_dir/merged_${SAMPLE_ID_1}_${SAMPLE_ID_2}-filtered_snp.norm.pureSNP.vcf"
fi



# output
output_base="/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-05_snpEFF_annotation_merged_vcf"

if [ -n "$SAMPLE_ID_3" ]; then
	output_dir="$output_base/merged_${SAMPLE_ID_1}_${SAMPLE_ID_2}_${SAMPLE_ID_3}"
	out_filebase="merged_${SAMPLE_ID_1}_${SAMPLE_ID_2}_${SAMPLE_ID_3}"
else
	output_dir="$output_base/merged_${SAMPLE_ID_1}_${SAMPLE_ID_2}"
	out_filebase="merged_${SAMPLE_ID_1}_${SAMPLE_ID_2}"
fi

mkdir -p "$output_dir"





## parameters for running annotation by snpEFF :
ann_in=$input_merged_vcf  # input vcf
updown_stream_dist=2000 # set the up and down stream range: 2000bp
snpeff_config="/cfs/klemming/projects/snic/naiss2025-23-66/software/snpEFF/snpEff.config"

ann_out="$output_dir/${out_filebase}-filtered_snp.norm.pureSNP.snpeff.ann.vcf" #output annotation results
ann_csvStats_out="$output_dir/${out_filebase}-filtered_snp.norm.pureSNP.snpeff.ann.stats.csv" # show statistics results of annotation in csv
ann_htmlStats_out="$output_dir/${out_filebase}-filtered_snp.norm.pureSNP.snpeff.ann.stats.html"  # show statistics results of annotation in web html


## running annotation by snpEFF 
java -Xmx16g -jar /sw/bioinfo/snpEff/5.2/rackham/snpEff.jar \
-c $snpeff_config \
-ud $updown_stream_dist \
-csvStats $ann_csvStats_out \
-htmlStats $ann_htmlStats_out \
-o vcf \
HCON_2020 \
$ann_in \
> $ann_out


echo "Completed snpeff annotation for merged samples SNP "

```

### Submit the SLURM job
```bash
# annotate the merged vcfs from a pair of samples: 
# NOTE: the order of two samples name as arguments, should be same for all commands

#a pair of samples: s15, s4 
sbatch 01_snpEFF_anno_merged.sh Sample_XJ-4126-15 Sample_XJ-4126-4-ampure-1

# a pair of samples: s3, s4 
sbatch 01_snpEFF_anno_merged.sh Sample_XJ-4126-3 Sample_XJ-4126-4-ampure-1

### a pair of farm samples: s4, G (replicates s7 + s9)
## as farm G has 2 replicates samples s7 and s9,
sbatch 01_snpEFF_anno_merged.sh Sample_XJ-4126-4-ampure-1 Sample_XJ-4126-7 Sample_XJ-4126-9  

```

### running time:
5-6mins

### input
```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/
025_snpEFF_vcf_annota/025-04_merge_vcf/

merged_Sample_XJ-4126-15_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.vcf

merged_Sample_XJ-4126-3_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.vcf

merged_Sample_XJ-4126-4-ampure-1_Sample_XJ-4126-7_Sample_XJ-4126-9-filtered_snp.norm.pureSNP.vcf
```



### output 
one example for merged vcf from s15 vs s4:
```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-05_snpEFF_annotation_merged_vcf/
merged_Sample_XJ-4126-15_Sample_XJ-4126-4-ampure-1> du -sh *

104K	merged_Sample_XJ-4126-15_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.snpeff.ann.stats.csv # # statistics csv, the content is same to html
2.7M	merged_Sample_XJ-4126-15_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.snpeff.ann.stats.genes.txt # only include the annotation of variants on genes
388K	merged_Sample_XJ-4126-15_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.snpeff.ann.stats.html  # statistics html
5.8G	merged_Sample_XJ-4126-15_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.snpeff.ann.vcf # annotation info


## content example:  ann.vcf is the annotation info file
## ANN= is the annotation information...
ANN=T|intergenic_region|MODIFIER|CHR_START-HCON_00095840|CHR_START-gene:HCON_00095840|intergenic_region|CHR_START-gen
```


## 2. Annotate each SNP VCF for each sample

1. parameterized script:
```bash
01_snpEFF_annotation.sh

```

``` bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -n 1
#SBATCH -c 10          # 16g/2g = 8 cores
#SBATCH -t 10:00:00
#SBATCH -p shared
#SBATCH -J eff_anno
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se


module load bioinfo-tools
module load snpEff

# Use the sample ID passed as the first argument
SAMPLE_ID=$1
select_snp_indel=$2   # you can select `snp`, or `indel` file

# Define paths (adjust these if necessary)
input_dir='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/024_freebayes_vcf_filter_2'

output_base='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-02_snpEFF_annotation'

# Construct file paths:
# Sample_XJ-4126-1-ampure_filtered_indel.vcf, Sample_XJ-4126-1-ampure_filtered_snp.vcf, etc.
sample_snp_vcf="${input_dir}/${SAMPLE_ID}_filtered_${select_snp_indel}.vcf"
output_dir="${output_base}/${SAMPLE_ID}"
mkdir -p "$output_dir"


## parameters for running annotation
ann_in=$sample_snp_vcf  # input vcf
ann_out="$output_dir/${SAMPLE_ID}_filtered_${select_snp_indel}.snpeff.ann.vcf" #output annotation results
ann_csvStats_out="$output_dir/${SAMPLE_ID}_filtered_${select_snp_indel}.snpeff.ann.stats.csv" # show statistics results of annotation in csv
ann_htmlStats_out="$output_dir/${SAMPLE_ID}_filtered_${select_snp_indel}.snpeff.ann.stats.html"  # show statistics results of annotation in web html
updown_stream_dist=2000 # set the up and down stream range: 2000bp
snpeff_config="/cfs/klemming/projects/snic/naiss2025-23-66/software/snpEFF/snpEff.config"

## running annotation
java -Xmx16g -jar /sw/bioinfo/snpEff/5.2/rackham/snpEff.jar \
-c $snpeff_config \
-ud $updown_stream_dist \
-csvStats $ann_csvStats_out \
-htmlStats $ann_htmlStats_out \
-o vcf \
HCON_2020 \
$ann_in \
> $ann_out


echo "Completed snpeff annotation for $SAMPLE_ID ($select_snp_indel) "
```

2. wrapper script:
```bash
01_snpEFF_annotation_batch_run.sh
```

```bash
#!/bin/bash -l

# List the sample IDs
samples=(
  "Sample_XJ-4126-15"
  "Sample_XJ-4126-4-ampure-1"
  "Sample_XJ-4126-1-ampure"
  "Sample_XJ-4126-2-ampure"
  "Sample_XJ-4126-3"
  "Sample_XJ-4126-5"
  "Sample_XJ-4126-6"
  "Sample_XJ-4126-7"
  "Sample_XJ-4126-8"
  "Sample_XJ-4126-9"
  "Sample_XJ-4126-10"
  "Sample_XJ-4126-11"
  "Sample_XJ-4126-12"
  "Sample_XJ-4126-13"
  "Sample_XJ-4126-14"
)



# Loop over each sample ID and submit the SLURM job
for sample in "${samples[@]}"; do
	for v in snp indel; do         # do both snp and indel file 
    sbatch 01_snpEFF_annotation.sh $sample $v
    echo "Submitting job for $sample $v"
	done
done

```

### submit the slurm job

Run the wrapper to submit 15 jobs for all 15 samples at once:
```bash
bash 01_snpEFF_annotation.sh
```


### running time:
1-10mins


### input
one sample example:
```bash
analysis/variant_explore/illumina_HC_wgs_XJ4126/024_freebayes_vcf_filter_2

301M	Sample_XJ-4126-10_filtered_indel.vcf
3.0G	Sample_XJ-4126-10_filtered_snp.vcf
...
```



### output 
for each input vcf, there are 4 output files
```bash

/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-02_snpEFF_annotation

|-- Sample_XJ-4126-10
|   |-- Sample_XJ-4126-10_filtered_indel.snpeff.ann.stats.csv
|   |-- Sample_XJ-4126-10_filtered_indel.snpeff.ann.stats.genes.txt
|   |-- Sample_XJ-4126-10_filtered_indel.snpeff.ann.stats.html
|   |-- Sample_XJ-4126-10_filtered_indel.snpeff.ann.vcf

|   |-- Sample_XJ-4126-10_filtered_snp.snpeff.ann.stats.csv
|   |-- Sample_XJ-4126-10_filtered_snp.snpeff.ann.stats.genes.txt
|   |-- Sample_XJ-4126-10_filtered_snp.snpeff.ann.stats.html
|   -- Sample_XJ-4126-10_filtered_snp.snpeff.ann.vcf

```