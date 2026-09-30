# Extract VCF information for selecting filtering thresholds
we will mainly focus on `QUAL`, `DP`

## Tools
* bcftools
* vcflib



## SLURM bash script

### Create the bash script:

```bash
020_freebayes_extracted_vcf_info.sh
```

### Script content

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -n 1
#SBATCH -c 34     # 34 file in total, so use 34 cores here
#SBATCH -t 01:00:00  ## add this
#SBATCH -p shared
#SBATCH -J et_info
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se

module load bioinfo-tools vcflib/2017-04-04
module load bcftools/1.20

INPUT_DIR=".../analysis/variant_explore/illumina_HC_wgs_XJ4126/019_freebayes_vc_merged_separate"

OUTPUT_DIR=".../analysis/variant_explore/illumina_HC_wgs_XJ4126/020_freebayes_extracted_vcf_info"

mkdir -p "$OUTPUT_DIR"


# header
HEADER="CHROM\tPOS\tREF\tALT\tQUAL\tAF\tDP\tSRP\tSAP\tEPP\tGT\tRO\tAO"

# Loop over both SNP and INDEL VCFs
for vcf in "$INPUT_DIR"/*.vcf; do
  ( echo "processing $vcf"

    # Get the base name without the .vcf extension
    fname=$(basename "$vcf")
    sample_name="${fname%%.vcf}"
    
    output_file="$OUTPUT_DIR/${sample_name}_extracted_info.tsv"
 
    # write header to this output file
    echo -e $HEADER > $output_file

    # we did not use this: this can use more 30 mins
    # Randomly sample 2% of variants (change -r 0.02 to -r 0.01 if you need 1%)
    # bcftools view "$vcf" | vcfrandomsample -r 0.02 > "$OUTPUT_DIR/${sample_name}.subset.vcf" # vcfrandomsample from the vcflib package

    # Extract key fields: CHROM, POS, REF, ALT, QUAL, INFO/AF, INFO/DP, INFO/SRP, INFO/SAP, INFO/EPP, FORMAT/GT, FORMAT/RO, FORMAT/AO 
    bcftools query \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%AF\t%DP\t%SRP\t%SAP\t%EPP\t[%GT\t%RO\t%AO]\n' \
      $vcf \
      >> $output_file      #>> : Append to the file
  ) &
done

wait
echo "All completed"

## we extract more fields just in case
```

### Submit the SLURM job
```bash
bash 020_freebayes_extracted_vcf_info.sh
```


## Running time:
for all 17 samples: ~2min


## Input 

all 17 samples, e.g. for each sample, one `.indel.vcf`, one `.snp.vcf`

```bash
019_freebayes_vc_merged_separate/> ll

443M	Sample_XJ-4126-10.indel.vcf
5.1G	Sample_XJ-4126-10.snp.vcf 




## Output

for each sample, output one `.indel_extracted_info.tsv`, one `.snp_extracted_info.tsv`:

```bash
analysis/variant_explore/illumina_HC_wgs_XJ4126/020_freebayes_extracted_vcf_info

97M	Sample_XJ-4126-10.indel_extracted_info.tsv
1.1G	Sample_XJ-4126-10.snp_extracted_info.tsv

82M	Sample_XJ-4126-11.indel_extracted_info.tsv
935M	Sample_XJ-4126-11.snp_extracted_info.tsv

92M	Sample_XJ-4126-12.indel_extracted_info.tsv
1.1G	Sample_XJ-4126-12.snp_extracted_info.tsv

103M	Sample_XJ-4126-13.indel_extracted_info.tsv
1.1G	Sample_XJ-4126-13.snp_extracted_info.tsv

103M	Sample_XJ-4126-14.indel_extracted_info.tsv
1.1G	Sample_XJ-4126-14.snp_extracted_info.tsv

81M	Sample_XJ-4126-15.indel_extracted_info.tsv
875M	Sample_XJ-4126-15.snp_extracted_info.tsv

3.2M	Sample_XJ-4126-16.indel_extracted_info.tsv
102M	Sample_XJ-4126-16.snp_extracted_info.tsv

29M	Sample_XJ-4126-17.indel_extracted_info.tsv
440M	Sample_XJ-4126-17.snp_extracted_info.tsv

105M	Sample_XJ-4126-1-ampure.indel_extracted_info.tsv
1.1G	Sample_XJ-4126-1-ampure.snp_extracted_info.tsv

102M	Sample_XJ-4126-2-ampure.indel_extracted_info.tsv
1.1G	Sample_XJ-4126-2-ampure.snp_extracted_info.tsv

83M	Sample_XJ-4126-3.indel_extracted_info.tsv
885M	Sample_XJ-4126-3.snp_extracted_info.tsv

78M	Sample_XJ-4126-4-ampure-1.indel_extracted_info.tsv
793M	Sample_XJ-4126-4-ampure-1.snp_extracted_info.tsv

100M	Sample_XJ-4126-5.indel_extracted_info.tsv
1.1G	Sample_XJ-4126-5.snp_extracted_info.tsv

104M	Sample_XJ-4126-6.indel_extracted_info.tsv
1.1G	Sample_XJ-4126-6.snp_extracted_info.tsv

90M	Sample_XJ-4126-7.indel_extracted_info.tsv
941M	Sample_XJ-4126-7.snp_extracted_info.tsv

89M	Sample_XJ-4126-8.indel_extracted_info.tsv
944M	Sample_XJ-4126-8.snp_extracted_info.tsv

76M	Sample_XJ-4126-9.indel_extracted_info.tsv
777M	Sample_XJ-4126-9.snp_extracted_info.tsv
```


