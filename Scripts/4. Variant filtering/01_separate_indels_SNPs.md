# Separate SNPs and indels 
will analyse SNPs and indels separately in downstream

## Tools
* bcftools


## SLURM bash script

### Create the bash script:

```bash
019_freebayes_vc_merged_separate.sh
```

### Script content

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -n 1          # one task
#SBATCH -c 17         # 1 thread per sample, 17 cores in total
#SBATCH -t 01:00:00
#SBATCH -p shared
#SBATCH -J fb2_separate
#SBATCH --mail-type=All,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se

module load bcftools/1.20

output_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/019_freebayes_vc_merged_separate'

input_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/016_freebayes_vc_merge_index'

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

for vc in "$input_dir"/*.vcf; do
	(
		echo "processing $vc"

		sample=$(basename "$vc")
		sample_name=${sample%%.vcf}

		## Separate SNPs
		bcftools view \
			--types snps \
			--output-type v \
			--output "$output_dir/${sample_name}.snp.vcf" \
			"$vc"

		## Separate indels
		bcftools view \
			--types indels \
			--output-type v \
			--output "$output_dir/${sample_name}.indel.vcf" \
			"$vc"
	) &  # Run each sample in parallel
done

wait  # Wait for all background processes to finish.

# bcftools --types， Types are determined by comparing the REF and ALT alleles in the VCF record not INFO tags

# --output-type， Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v).

```


## Running time:
for all 17 samples: ~4min


## Input 

all 17 samples, e.g. one sample example

```bash
016_freebayes_vc_merge_index> ll

7.1G	Sample_XJ-4126-1-ampure.vcf  # input
1.6G	Sample_XJ-4126-1-ampure.vcf.gz
272K	Sample_XJ-4126-1-ampure.vcf.gz.tbi




## Output

for each sample, output one `.indel.vcf`, one `.snp.vcf`:

```bash
analysis/variant_explore/illumina_HC_wgs_XJ4126/019_freebayes_vc_merged_separate/

443M	Sample_XJ-4126-10.indel.vcf
5.1G	Sample_XJ-4126-10.snp.vcf

375M	Sample_XJ-4126-11.indel.vcf
4.5G	Sample_XJ-4126-11.snp.vcf

422M	Sample_XJ-4126-12.indel.vcf
5.0G	Sample_XJ-4126-12.snp.vcf

473M	Sample_XJ-4126-13.indel.vcf
5.3G	Sample_XJ-4126-13.snp.vcf

476M	Sample_XJ-4126-14.indel.vcf
5.3G	Sample_XJ-4126-14.snp.vcf

373M	Sample_XJ-4126-15.indel.vcf
4.3G	Sample_XJ-4126-15.snp.vcf

16M	Sample_XJ-4126-16.indel.vcf
508M	Sample_XJ-4126-16.snp.vcf

132M	Sample_XJ-4126-17.indel.vcf
2.2G	Sample_XJ-4126-17.snp.vcf

482M	Sample_XJ-4126-1-ampure.indel.vcf
5.4G	Sample_XJ-4126-1-ampure.snp.vcf

472M	Sample_XJ-4126-2-ampure.indel.vcf
5.2G	Sample_XJ-4126-2-ampure.snp.vcf

383M	Sample_XJ-4126-3.indel.vcf
4.3G	Sample_XJ-4126-3.snp.vcf

358M	Sample_XJ-4126-4-ampure-1.indel.vcf
3.9G	Sample_XJ-4126-4-ampure-1.snp.vcf

462M	Sample_XJ-4126-5.indel.vcf
5.1G	Sample_XJ-4126-5.snp.vcf

479M	Sample_XJ-4126-6.indel.vcf
5.3G	Sample_XJ-4126-6.snp.vcf

412M	Sample_XJ-4126-7.indel.vcf
4.6G	Sample_XJ-4126-7.snp.vcf

408M	Sample_XJ-4126-8.indel.vcf
4.6G	Sample_XJ-4126-8.snp.vcf

351M	Sample_XJ-4126-9.indel.vcf
3.8G	Sample_XJ-4126-9.snp.vcf
```


