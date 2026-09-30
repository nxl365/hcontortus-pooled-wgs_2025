# filter variants based on QUAL and DP thresholds
* removed sample 16 and 17 already, so now 15 samples
* QUAL: 20 , DP:   10  

## Tools
* bcftools
* vcflib



## SLURM bash script

### Create the bash script:

```bash
024_freebayes_vcf_filter_2.sh
024_freebayes_vcf_filter_2_batch_run.sh # then only need to submit this
```

### Script content

1. Parameterized script: `024_freebayes_vcf_filter_2.sh`

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -n 1
#SBATCH -c 30     # 15 farms, 30 files in total, so use 30 cores here
#SBATCH -t 00:30:00  ## add this
#SBATCH -p shared
#SBATCH -J filt_vcf2
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se

module load bioinfo-tools vcflib/2017-04-04
module load bcftools/1.20

INPUT_DIR=".../analysis/variant_explore/illumina_HC_wgs_XJ4126/019_freebayes_vc_merged_separate"

OUTPUT_DIR=".../analysis/variant_explore/illumina_HC_wgs_XJ4126/024_freebayes_vcf_filter_2"

mkdir -p "$OUTPUT_DIR"


sample_ID=$1

input_sample_indel="$INPUT_DIR/$sample_ID.indel.vcf"
input_sample_snp="$INPUT_DIR/$sample_ID.snp.vcf"

# keep variants with both QUAL>=20 AND DP>=5 in any sample
bcftools filter -i 'QUAL>=20 & INFO/DP>=10' $input_sample_indel > "$OUTPUT_DIR/${sample_ID}_filtered_indel.vcf"
bcftools filter -i 'QUAL>=20 & INFO/DP>=10' $input_sample_snp > "$OUTPUT_DIR/${sample_ID}_filtered_snp.vcf"


# FMT/DP or INFO/DP should be same, we only have one sample in one vcf file

# remember to check bcftool `&& vs & and || vs |` to make sure you use them correctly
```

2. Wrapper script for submitting all 15 samples at once:`024_freebayes_vcf_filter_2_batch_run.sh`

```bash
#!/bin/bash -l

# List the sample IDs based on your BAM file names (exclude files like metrics or index files)
samples=(
  "Sample_XJ-4126-1-ampure"
  "Sample_XJ-4126-2-ampure"
  "Sample_XJ-4126-3"
  "Sample_XJ-4126-4-ampure-1"
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
  "Sample_XJ-4126-15"
)

# Loop over each sample ID and submit the SLURM job
for sample in "${samples[@]}"; do
    echo "Submitting job for sample: $sample"
    sbatch 024_freebayes_vcf_filter_2.sh "$sample"
done


```

### Submit the SLURM job

1. if you want to run a single sample
```bash
# $1 can be given as a parameter <sample name> when you run it:
sbatch 024_freebayes_vcf_filter_2.sh [Sample1_name]
```

2. Run the wrapper to submit 15 jobs for all 15 samples at once:
```bash
bash 024_freebayes_vcf_filter_2.sh
```

## Running time:
15 samples from 1 — 2 mins


## Input 

for each sample, one `.indel.vcf`, one `.snp.vcf`:

```bash
analysis/variant_explore/illumina_HC_wgs_XJ4126/019_freebayes_vc_merged_separate/

443M	Sample_XJ-4126-10.indel.vcf
5.1G	Sample_XJ-4126-10.snp.vcf
...
```


## Output

* for each sample, 2 outputs

```bash
analysis/variant_explore/illumina_HC_wgs_XJ4126/024_freebayes_vcf_filter_2

301M	Sample_XJ-4126-10_filtered_indel.vcf
3.0G	Sample_XJ-4126-10_filtered_snp.vcf
...
```