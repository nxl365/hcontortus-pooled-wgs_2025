# Normalise vcf file

We did 3 steps:
1. normalise VCFs
  Normalising a VCF usually means splitting multi-allelic calls, left-aligning indels, and ensuring that the alleles specified as REF in your VCF are also present in the reference genome.

2. extract pure SNP from snp.vcf with multi-allelic sites
* freebayes can call MNP, Indels, SNP, complex variants type...
* Therefore, even after extracting SNPs with `bcftools view --types snps` previously, there are still some multiallelic records containing other variant types if at least one ALT allele is a SNP.
* e.g.  
```bash
POS    REF  ALT
2982   CG   TA,TG
```

3. index bcf or vcf.gz


## SLURM bash script

### Create the bash script:

```bash
01_norm_pureSNP_vcf.sh
01_norm_pureSNP_vcf_batch_run.sh

```

### Script content

1. parameterized script: `01_norm_pureSNP_vcf.sh`
```bash
#!/bin/bash -l
#SBATCH -A naiss2026-3-101
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 00:30:00 # 30min
#SBATCH -p shared
#SBATCH -J norm_vcf
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se


module load bcftools/1.20

# Use the sample ID passed as the first argument
SAMPLE_ID=$1
select_snp_indel=$2   # `snp`, or `indel`


refer_genome='/cfs/klemming/projects/snic/naiss2025-23-66/variants_explore_2025-2/raw_data/Hc_reference_genome/haemonchus_contortus.PRJEB506.WBPS19.genomic.fa'

input_dir='/cfs/klemming/projects/snic/naiss2025-23-66/variants_explore_2025-2/analysis/variant_explore/illumina_HC_wgs_XJ4126/024_freebayes_vcf_filter_2'

output_base='/cfs/klemming/projects/snic/naiss2025-23-66/variants_explore_2025-2/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-03_normalise_vcf'


# Construct file paths:
# Sample_XJ-4126-1-ampure_filtered_indel.vcf, Sample_XJ-4126-1-ampure_filtered_snp.vcf, etc.
sample_snp_vcf="${input_dir}/${SAMPLE_ID}_filtered_${select_snp_indel}.vcf"
output_dir="${output_base}/${SAMPLE_ID}"
mkdir -p "$output_dir"


####### STEP1 normalising:  #######
#bcftools norm 
# 1. ´-m -both´: splits multi-allelic calls into separate variant calls (´-m -any´ can work as well, these two parameters only differ when merging, not when splitting, Here we are splitting, so no matter)

# 2. ´--check-ref w´: check the REF allele against the FASTA reference genome, and warn on mismatches.
# ´-f´: use this FASTA for checking REF alleles and for left-aligning/normalizing indels.
## but this is for indels and we will only focus on SNP later, so no matter now

# 3. ´-Ob´: output bcf, since I want to be faster for next merging process

bcftools norm --threads 8  -m-both $sample_snp_vcf \
| bcftools norm -Ob --check-ref w -f $refer_genome \
> "$output_dir/${SAMPLE_ID}_filtered_${select_snp_indel}.norm.bcf"


# just for viewing: convert BCF → compressed VCF
bcftools view -Oz "$output_dir/${SAMPLE_ID}_filtered_${select_snp_indel}.norm.bcf" \
-o "$output_dir/${SAMPLE_ID}_filtered_${select_snp_indel}.norm.vcf.gz"



####### STEP2  extract pure SNPs:  #######
input_bcf="$output_dir/${SAMPLE_ID}_filtered_${select_snp_indel}.norm.bcf"
output_bcf= "$output_dir/${SAMPLE_ID}_filtered_${select_snp_indel}.norm.pureSNP.bcf"
output_vcfgz= "$output_dir/${SAMPLE_ID}_filtered_${select_snp_indel}.norm.pureSNP.vcf.gz"

bcftools view --threads 8 -v snps -Ob -o $output_bcf $input_bcf

# just for viewing: convert BCF → compressed VCF
bcftools view -Oz $output_bcf -o $output_vcfgz



####### STEP3  index bcf  ####### for next step
bcftools index --threads 8 -f $output_bcf  # -f: it forces re-indexing even if an index file already exists

```

2. wrapper: `01_norm_pureSNP_vcf_batch_run.sh`
```bash
#!/bin/bash -l

# List the sample IDs

samples=(
"Sample_XJ-4126-15" 
"Sample_XJ-4126-4-ampure-1" 
"Sample_XJ-4126-3" 
"Sample_XJ-4126-7" 
"Sample_XJ-4126-9" 
"Sample_XJ-4126-8"
 "Sample_XJ-4126-1-ampure"
 "Sample_XJ-4126-2-ampure"
 "Sample_XJ-4126-5"
 "Sample_XJ-4126-6"
 "Sample_XJ-4126-10"
 "Sample_XJ-4126-11"
 "Sample_XJ-4126-12"
 "Sample_XJ-4126-13"
 "Sample_XJ-4126-14"
)


# Loop over each sample ID and submit the SLURM job
for sample in "${samples[@]}"; do
	#for v in snp indel; do
    sbatch 01_norm_pureSNP_vcf.sh $sample snp  # we select only snp
    echo "Submitting job for $sample snp"
	#done
done

```



### Submit the SLURM job
Run the wrapper to submit jobs for all samples at once:
```bash
bash 01_norm_pureSNP_vcf_batch_run.sh
```



## Running time:
~4 mins



## Input 

```bash
analysis/variant_explore/illumina_HC_wgs_XJ4126/024_freebayes_vcf_filter_2

301M	Sample_XJ-4126-10_filtered_indel.vcf
3.0G	Sample_XJ-4126-10_filtered_snp.vcf    # we select only snp
...
```



## Output

did for all samples,   
e.g. in each sample’s folder:

1. after normalising
The multiallelic rows are split into separate rows.

```bash
025_snpEFF_vcf_annota/025-03_normalise_vcf/Sample_XJ-4126-15

565M	Sample_XJ-4126-15_filtered_snp.norm.bcf
588M	Sample_XJ-4126-15_filtered_snp.norm.vcf.gz
```

2. after extract pure SNPs, and index
```bash
459M	Sample_XJ-4126-15_filtered_snp.norm.pureSNP.bcf
180K	Sample_XJ-4126-15_filtered_snp.norm.pureSNP.bcf.csi

482M	Sample_XJ-4126-15_filtered_snp.norm.pureSNP.vcf.gz
```