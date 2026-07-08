# Statistics of VCF file

## Tools
* RTG

This step calculated summary statistics for each sample-level FreeBayes VCF file.


## SLURM bash script

### Create the bash script:
```bash
018_freebayes_vc_merged_rtg_stats.sh

018_freebayes_vc_merged_rtg_stats_batch_run.sh

```

### Script content

1. parameterized script `018_freebayes_vc_merged_rtg_stats.sh`:

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-***
#SBATCH -n 1          # one task
#SBATCH -c 8          # 16GB RAM per sample, 8 cores
#SBATCH -t 01:00:00
#SBATCH -p shared
#SBATCH -J fb2_stats
#SBATCH --mail-type=All,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se

module load bioinfo-tools
module load rtgtools/3.12.1

sample_ID=$1

output_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/018_freebayes_vc_merged_rtg_stats'

input_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/016_freebayes_vc_merge_index'

mkdir -p "$output_dir"

input_sample="$input_dir/$sample_ID.vcf"

rtg RTG_MEM=16G vcfstats $input_sample > "$output_dir/${sample_ID}_vcf_rtg_stats.out"  # generate stats

```

2. Wrapper script for submitting all 17 jobs for 17 samples at once:`018_freebayes_vc_merged_rtg_stats_batch_run.sh`

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
  "Sample_XJ-4126-16"
  "Sample_XJ-4126-17"
)

# Loop over each sample ID and submit the SLURM job
for sample in "${samples[@]}"; do
    echo "Submitting job for sample: $sample"
    sbatch 018_freebayes_vc_merged_rtg_stats.sh "$sample"
done

```


### Submit the SLURM job

Run the wrapper to submit 17 jobs for all 17 samples at once:
```bash
bash 018_freebayes_vc_merged_rtg_stats_batch_run.sh
```


## Running time:
for each sample/job: ~ 1min - 3min


## Input 
all 17 samples, e.g. one sample example

```bash
016_freebayes_vc_merge_index> ll

Sample_XJ-4126-1-ampure.vcf

```


## Output

17 samples in total; for each sample: 1 output file  

e.g.
```text

Sample_XJ-4126-1-ampure_vcf_rtg_stats.out

```

### Output explanation
one output `_rtg_stats.out` file example:

```bash
Location                     : ...
Failed Filters               : 0
Passed Filters               : 38900199
SNPs                         : 7300226
MNPs                         : 2715197
Insertions                   : 256760
Deletions                    : 311990
Indels                       : 647189
Same as reference            : 27668837
SNP Transitions/Transversions: 2.04 (5180582/2544744)
Total Het/Hom ratio          : 18.57 (10657450/573912)
SNP Het/Hom ratio            : 18.35 (6922897/377329)
MNP Het/Hom ratio            : 22.33 (2598797/116400)
Insertion Het/Hom ratio      : 10.91 (235198/21562)
Deletion Het/Hom ratio       : 9.71 (282869/29121)
Indel Het/Hom ratio          : 20.94 (617689/29500)
Insertion/Deletion ratio     : 0.82 (256760/311990)
Indel/SNP+MNP ratio          : 0.12 (1215939/10015423)
```