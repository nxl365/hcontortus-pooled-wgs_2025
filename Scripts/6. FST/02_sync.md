# generate the synchronised file
generate a synchronised file (containing allele frequency estimates for each sample at every genomic position) using Popoolation2 


## Bash script

### Create the bash script:

```bash
mkdir -p '/cfs/klemming/projects/snic/naiss2025-23-66/scripts/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-02_sync'

01_sync_by_chr.sh

```

### Script content

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -J sync
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -t 10:00:00
#SBATCH -p shared
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se

module load bioinfo-tools popoolation2/1201

mpileup='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-01_generate_mpileup/all_samples_removed_s16_s17_fillmiss.mpileup'

outdir='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-02_sync'

mkdir -p "$outdir" && cd "$outdir"


# Chromosome list as an array
chromosomes_list=(
  hcontortus_chr4_Celeg_TT_arrow_pilon
  hcontortus_chr5_Celeg_TT_arrow_pilon
  hcontortus_chr2_Celeg_TT_arrow_pilon
  hcontortus_chrX_Celeg_TT_arrow_pilon
  hcontortus_chr1_Celeg_TT_arrow_pilon
  hcontortus_chr3_Celeg_TT_arrow_pilon
  mitochondrion
)


## 1. generate a sync file with an entire mpileup file of all samples
echo "[$(date)] Step1: mpileup2sync" ## echo date and time for step1 in log

java -Xmx50g -jar /sw/bioinfo/popoolation2/1201/rackham/mpileup2sync.jar \
  --input "$mpileup" \
  --output all_samples_removed_s16_s17_fillmiss.sync \
  --min-qual 20 \
  --threads $SLURM_CPUS_ON_NODE \
  --fastq-type sanger

# $SLURM_CPUS_ON_NODE is  32, as #SBATCH -c 32



## 2. split the entire sync file into 7 chromosomes .sync files
echo "[$(date)] Step2: split by chromosome"
for chr in "${chromosomes_list[@]}"; do
  grep -P "^${chr}\t" all_samples_removed_s16_s17_fillmiss.sync > "${chr}_rm_s16_s17.sync" # ^ the 1st line, is chr name and \t
done

```
### Submit the SLURM job

```bash
sbatch ...sh
```

## running time
30 mins


## input
one mpileup file
```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-01_generate_mpileup/

all_samples_removed_s16_s17_fillmiss.mpileup
```


## ouput

We convert the entire mpileup file containing all samples (~280 GB) into a single sync file, then split it into 7 chromosome-specific sync files

```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-02_sync> du -sh *

62G	all_samples_removed_s16_s17_fillmiss.sync

9.9G	hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17.sync
11G	hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17.sync
9.6G	hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17.sync
12G	hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17.sync
11G	hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17.sync
9.8G	hcontortus_chrX_Celeg_TT_arrow_pilon_rm_s16_s17.sync
3.2M	mitochondrion_rm_s16_s17.sync
```


