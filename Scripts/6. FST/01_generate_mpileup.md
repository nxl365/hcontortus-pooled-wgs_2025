# generate mpileup
A multi-sample pileup file was generated from the 15 BAM files (from 15 samples, we already removed samples 16, 17) using SAMtools mpileup


## 1. make a BAM file list
Important! To provide samples in a specific order as input for pileup generation.


```bash
cd '/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-01_generate_mpileup'

# create bamlist_removed_s16_s17.txt
module load nano
nano bamlist_removed_s16_s17.txt

```


```bash
### create the content in bamlist_removed_s16_s17.txt
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-1-ampure.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-2-ampure.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-3.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-4-ampure-1.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-5.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-6.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-7.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-8.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-9.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-10.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-11.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-12.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-13.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-14.sort.markdup.bam
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate/Sample_XJ-4126-15.sort.markdup.bam

```
### Download file locally

```bash
### download this to local

cd variants-investigation-in-pooled-Hc/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-01_generate_mpileup

scp -r ...

```


## 2. generate mpileup


```bash
mkdir -p '/cfs/klemming/projects/snic/naiss2025-23-66/scripts/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-01_generate_mpileup'

# script name
06-01_generate_mpileup.sh

```

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -n 1          # one task
#SBATCH -c 1         # samtools mpileup Only can use one thread
#SBATCH -t 24:00:00
#SBATCH -p shared
#SBATCH -J mplieup
#SBATCH --mail-type=All,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se


refer_genome='/cfs/klemming/projects/snic/naiss2025-23-66/raw_data/Hc_reference_genome/haemonchus_contortus.PRJEB506.WBPS19.genomic.fa'

input_dir='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate'

output_dir='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-01_generate_mpileup'


mkdir -p "$output_dir"

# Load required modules
module load samtools/1.20


# Generate a pileup file:
# 15 samples' BAM files will be used as input one by one to generate a single large mpileup file.

samtools mpileup \
-d 500 \
--min-MQ 20 \
--min-BQ 30 \
-f $refer_genome \
-o "$output_dir/all_samples_removed_s16_s17.mpileup" \
-b "$output_dir/bamlist_removed_s16_s17.txt"  # give a certain order of samples as input



# turns empty fields in the pileup into explicit “!” placeholders
# so downstream tools popoolation2 don’t mis‑count missing columns and error

sed 's/\t\t/\t!\t!/g' "$output_dir/all_samples_removed_s16_s17.mpileup" > "$output_dir/all_samples_removed_s16_s17_fillmiss.mpileup"



## samtools mpileup
# --min-MQ, mapping quality ≥ 20
# --min-BQ, base quality ≥ 30, 
# -d  cap per‑site read depth at 500; any extra reads are ignored to save memory/time

```

### submit the slurm job
```bash
sbatch 06-01_generate_mpileup.sh
```

### running time:
14 hours


### input
the 15 samples' `.sort.markdup.bam` files listed in `bamlist_removed_s16_s17.txt`


### output

Important! The population order starting from the 4th column in the mpileup file should match the sample order specified in `bamlist_removed_s16_s17.txt`.

```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-01_generate_mpileup


275G	all_samples_removed_s16_s17_fillmiss.mpileup

275G	all_samples_removed_s16_s17.mpileup

4.0K	bamlist_removed_s16_s17.txt

```