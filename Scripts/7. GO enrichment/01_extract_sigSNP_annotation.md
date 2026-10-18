# extract significantly differentiated SNPs with annotations
For the three most divergent farm pairs (s4 VS s15,  s3 VS s4,  G (has 2 replicates s7, s9 ) VS s4), significantly differentiated SNPs with annotations were extracted from the three merged and annotated VCF files 

* Important: Make sure the POS column is consistent across all files (e.g., we have 1-based, where the first base of each chromosome is position 1).

## Bash script

1. parameterized script:

```bash

mkdir -p '/cfs/klemming/projects/snic/naiss2025-23-66/scripts/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-07_sig_SNP'

01_extract_sig_SNP_from_annoVCF.sh

```

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -J sig_snp
#SBATCH -c 8
#SBATCH -t 01:00:00 
#SBATCH -p shared
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se


module load bcftools

### arguments
SAMPLE_ID_1=$1
SAMPLE_ID_2=$2

windowFst_cut=$3
everysnpFst_cut=$4

## input1
input_anno_vcf=$5  

### input2
input_sig_snp_pos="/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-06_poola_fet_by_chr_everySNP/sig_snp_pos_list/${SAMPLE_ID_1}_${SAMPLE_ID_2}_sig_snp_pos_list-windowFst_10000_Cut${windowFst_cut}-everysnpFst_cut${everysnpFst_cut}.txt"

### output
output_dir="/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-07_sig_SNP/${SAMPLE_ID_1}_VS_${SAMPLE_ID_2}/windowFst_10000_Cut${windowFst_cut}-everysnpFst_cut${everysnpFst_cut}"

mkdir -p $output_dir

output_sig_snp_vcf="$output_dir/${SAMPLE_ID_1}_${SAMPLE_ID_2}_sig_snp_anno_cut${windowFst_cut}_${everysnpFst_cut}.vcf"


# Compress + index VCF
bgzip -c "$input_anno_vcf" > "${input_anno_vcf}.gz"
bcftools index --tbi "${input_anno_vcf}.gz" # index vcf just in case

# Get the intersecting sig SNPs in annotated vcf 
bcftools view --threads 8 -R "$input_sig_snp_pos" -Ov -o "$output_sig_snp_vcf" "${input_anno_vcf}.gz"


```


2. Wrapper script:

```bash
01_extract_sig_SNP_from_annoVCF-batch_run.sh

```

```bash
#!/bin/bash -l

##### for different sample pairs


## for s4 s15
sbatch 01_extract_sig_SNP_from_annoVCF.sh \
s4 \
s15 \
1 \
1 \
"/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-05_snpEFF_annotation_merged_vcf/merged_Sample_XJ-4126-15_Sample_XJ-4126-4-ampure-1/merged_Sample_XJ-4126-15_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.snpeff.ann.vcf"


## for s3 s4
sbatch 01_extract_sig_SNP_from_annoVCF.sh \
s3 \
s4 \
1 \
1 \
"/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-05_snpEFF_annotation_merged_vcf/merged_Sample_XJ-4126-3_Sample_XJ-4126-4-ampure-1/merged_Sample_XJ-4126-3_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.snpeff.ann.vcf"


## for G s4
sbatch 01_extract_sig_SNP_from_annoVCF.sh \
G \
s4 \
1 \
1 \
"/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-05_snpEFF_annotation_merged_vcf/merged_Sample_XJ-4126-4-ampure-1_Sample_XJ-4126-7_Sample_XJ-4126-9/merged_Sample_XJ-4126-4-ampure-1_Sample_XJ-4126-7_Sample_XJ-4126-9-filtered_snp.norm.pureSNP.snpeff.ann.vcf"

```

### submit the slurm job
Run the wrapper to run 3 jobs for 3 pairs of samples
```bash
bash 01_extract_sig_SNP_from_annoVCF-batch_run.sh
```



### running time
2 mins


## input

```bash

merged_Sample_XJ-4126-15_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.snpeff.ann.vcf

merged_Sample_XJ-4126-3_Sample_XJ-4126-4-ampure-1-filtered_snp.norm.pureSNP.snpeff.ann.vcf

merged_Sample_XJ-4126-4-ampure-1_Sample_XJ-4126-7_Sample_XJ-4126-9-filtered_snp.norm.pureSNP.snpeff.ann.vcf


G_s4_sig_snp_pos_list-windowFst_10000_Cut1-everysnpFst_cut1.txt

s3_s4_sig_snp_pos_list-windowFst_10000_Cut1-everysnpFst_cut1.txt

s4_s15_sig_snp_pos_list-windowFst_10000_Cut1-everysnpFst_cut1.txt

```

## output
uploaded these outputs to `Analysis_results/7.GO enrichment`
```bash

# analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-07_sig_SNP/s4_VS_s15/windowFst_10000_Cut1-everysnpFst_cut1

s4_s15_sig_snp_anno_cut1_1.vcf


# analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-07_sig_SNP/s3_VS_s4/windowFst_10000_Cut1-everysnpFst_cut1

s3_s4_sig_snp_anno_cut1_1.vcf


# analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-07_sig_SNP/G_VS_s4/windowFst_10000_Cut1-everysnpFst_cut1

G_s4_sig_snp_anno_cut1_1.vcf

```
