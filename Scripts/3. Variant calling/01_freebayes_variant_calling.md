# Variants calling using Freebayes

## Tools
* freebayes

- This step uses FreeBayes to call variants from the duplicate-marked BAM files. Variant calling was performed separately for each sample and each chromosome, including six nuclear chromosomes and the mitochondrial genome.

- FreeBayes was used in pooled-continuous mode to call variants from pooled sequencing data.


## SLURM bash script

### Create the bash script:

```bash
015_freebayes_vc.sh
015_freebayes_vc_batch_run.sh # then only need to submit this
```

### Script content

1. Parameterized script: `015_freebayes_vc.sh`

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-***
#SBATCH -n 1          # one task
#SBATCH -c 10         # need 7 threads, but we book 10 in case
#SBATCH -t 10:00:00 
#SBATCH -p shared   
#SBATCH -J fb_vc2
#SBATCH --mail-type=All,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se  


# Use the sample ID passed as the first argument
SAMPLE_ID=$1

refer_genome='.../raw_data/Hc_reference_genome/haemonchus_contortus.PRJEB506.WBPS19.genomic.fa'

input_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate'

output_base='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/015_freebayes_vc'

# Construct file paths assuming BAM files are named like:
# Sample_XJ-4126-1-ampure.sort.markdup.bam, Sample_XJ-4126-2-ampure.sort.markdup.bam, etc.
sample_bam="${input_dir}/${SAMPLE_ID}.sort.markdup.bam"
output_dir="${output_base}/${SAMPLE_ID}"

mkdir -p "$output_dir"

# Load required modules
module load PDC/23.12
module load R/4.4.0
module load freebayes/1.3.7

# List of chromosomes to process
chromosomes_list="hcontortus_chr4_Celeg_TT_arrow_pilon
hcontortus_chr5_Celeg_TT_arrow_pilon
hcontortus_chr2_Celeg_TT_arrow_pilon
hcontortus_chrX_Celeg_TT_arrow_pilon
hcontortus_chr1_Celeg_TT_arrow_pilon
hcontortus_chr3_Celeg_TT_arrow_pilon
mitochondrion"

# Run FreeBayes for each chromosome in parallel
for chr in $chromosomes_list; do
    echo "Processing chromosome $chr for ${SAMPLE_ID}"
    
    freebayes \
        -f "$refer_genome" \
        -r "$chr" \
        -T 0.001 \   ## just same as the default 
        -p 2 \       ## just same as the default
        --pooled-continuous \
        -m 20 \
        -q 30 \
        "$sample_bam" \
        > "$output_dir/${SAMPLE_ID}.${chr}.vcf" \
        2> "$output_dir/${SAMPLE_ID}.${chr}.FB.log" &
done

# Wait for all FreeBayes processes to complete
wait

echo "All FreeBayes jobs for sample ${SAMPLE_ID} completed!"


# freebayes only can use one thread, per sample run 7 chromosomes in parallel in background in loop, 
# so need 7 threads in total

```

2. Wrapper script for submitting all 17 jobs for 17 samples at once:`015_freebayes_vc_batch_run.sh`
```bash
#!/bin/bash -l

# List the sample IDs based on BAM file names (exclude files like metrics or index files)
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
  "Sample_XJ-4126-16"
  "Sample_XJ-4126-17"
)

# Loop over each sample ID and submit the SLURM job
for sample in "${samples[@]}"; do
    echo "Submitting job for sample: $sample"
    sbatch 015_freebayes_vc.sh "$sample"
done

```

### code explanation**
* refer to the source: [Galaxy training - Calling variants in non-diploid systems](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/non-dip/tutorial.html)

```bash
##Population model##
-T 0.001 # default: 0.001, The expected mutation rate or pairwise nucleotide diversity among the population under analysis.
-p 2 # default: 2, Set ploidy of your organism

--pooled-continuous  # generate frequency-based variant calls, Output all alleles which pass input filters, regardless of genotyping outcome or model

## Input filters ## 
-m 20 # default: 1, Exclude alignments from analysis if they have a mapping quality less than Q. FreeBayes only consider reliably aligned reads.
-q 30 # default: 0, Exclude alleles from analysis if their supporting base quality is less than Q. FreeBayes only consider high quality bases.


##  -F and -C , also by default!
   -F --min-alternate-fraction N
                   Require at least this fraction of observations supporting
                   an alternate allele within a single individual in the
                   in order to evaluate the position.  default: 0.05
   -C --min-alternate-count N
                   Require at least this count of observations supporting
                   an alternate allele within a single individual in order
                   to evaluate the position.  default: 2
```


### Submit the SLURM job

1. if you want to run a single sample
```bash
# $1 can be given as a parameter <sample name> when you run it:
sbatch 015_freebayes_vc.sh [Sample1_name]
```

2. Run the wrapper to submit 17 jobs for all 17 samples at once:
```bash
bash 015_freebayes_vc_batch_run.sh
```

## Running time:
for each sample: ~ 4min - 1h40min


## Input 

1. reference genome, index:
```text
haemonchus_contortus.PRJEB506.WBPS19.genomic.fa
haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.fai
```

2. for each sample: bam, bai (index):

```text
6_mark_duplicate/

sample1.sort.markdup.bam      #bam, after marking/removing duplicates

sample1.sort.markdup.bai     # index file

```
* It’s good to have the .bam and .bai named to the same “prefix”!

* size: each sample's bam file after marking duplicates  8.6-14 GB (become a bit bigger )


## Output

each sample has a output folder, 

```bash
.../015_freebayes_vc> du -sh *

6.7G	Sample_XJ-4126-10
6.1G	Sample_XJ-4126-11
6.6G	Sample_XJ-4126-12
7.0G	Sample_XJ-4126-13
7.0G	Sample_XJ-4126-14
5.9G	Sample_XJ-4126-15
689M	Sample_XJ-4126-16
2.9G	Sample_XJ-4126-17
7.1G	Sample_XJ-4126-1-ampure
7.0G	Sample_XJ-4126-2-ampure
5.9G	Sample_XJ-4126-3
5.4G	Sample_XJ-4126-4-ampure-1
6.9G	Sample_XJ-4126-5
7.1G	Sample_XJ-4126-6
6.2G	Sample_XJ-4126-7
6.2G	Sample_XJ-4126-8
5.3G	Sample_XJ-4126-9
```
e.g. in each sample folder, `.vcf` file and `.log` file generated separately for 7 chromosomes 

```bash
Sample_XJ-4126-10.hcontortus_chr1_Celeg_TT_arrow_pilon.FB.log
Sample_XJ-4126-10.hcontortus_chr1_Celeg_TT_arrow_pilon.vcf

Sample_XJ-4126-10.hcontortus_chr2_Celeg_TT_arrow_pilon.FB.log
Sample_XJ-4126-10.hcontortus_chr2_Celeg_TT_arrow_pilon.vcf

Sample_XJ-4126-10.hcontortus_chr3_Celeg_TT_arrow_pilon.FB.log
Sample_XJ-4126-10.hcontortus_chr3_Celeg_TT_arrow_pilon.vcf

Sample_XJ-4126-10.hcontortus_chr4_Celeg_TT_arrow_pilon.FB.log
Sample_XJ-4126-10.hcontortus_chr4_Celeg_TT_arrow_pilon.vcf

Sample_XJ-4126-10.hcontortus_chr5_Celeg_TT_arrow_pilon.FB.log
Sample_XJ-4126-10.hcontortus_chr5_Celeg_TT_arrow_pilon.vcf

Sample_XJ-4126-10.hcontortus_chrX_Celeg_TT_arrow_pilon.FB.log
Sample_XJ-4126-10.hcontortus_chrX_Celeg_TT_arrow_pilon.vcf

Sample_XJ-4126-10.mitochondrion.FB.log
Sample_XJ-4126-10.mitochondrion.vcf
```

* FreeBayes only sends error messages or warnings to stderr. If everything runs smoothly, no messages are printed, leaving your .FB.log file empty.
