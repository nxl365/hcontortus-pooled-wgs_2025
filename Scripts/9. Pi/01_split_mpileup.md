# Split mpileup files by sample and chromosome

## Bash script

### Create the bash script:
```bash
mkdir -p /cfs/klemming/projects/snic/naiss2025-23-66/scripts/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-01_split_mpileup 


# script
07-01_split_mpileup.sh
```


### Script content

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -J npstats
#SBATCH -n 1
#SBATCH -c 20        # 20 cores
#SBATCH -t 10:00:00
#SBATCH -p shared
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se

mpileup='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-01_generate_mpileup/all_samples_removed_s16_s17_fillmiss.mpileup' 

bam_list="/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-01_generate_mpileup/bamlist_removed_s16_s17.txt"

out_dir="/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-01_split_mpileup"

mkdir -p $out_dir

# Move into output directory
cd "$out_dir"


### 1: Extract position columns and data columns
# Keep columns 1–3: chromosome, position, reference base
cut -f1-3 "$mpileup" > mpileup.positions
cut -f4- "$mpileup" > mpileup.data   # only step which is time consuming, others will be fast


### 2: Define columns per sample
# in our mpileup file, each sample occupies 3 columns (e.g., DP, bases, quality)
sample_cols=3


### 3: get per sample.mpileup
sample_index=1

# Read bamlist.txt line by line
while read bam_path; do
  # Extract sample name from file path (remove folder path and file extension)
  sample_name=$(basename "$bam_path" | sed 's/.sort.markdup.bam//')
  echo "[INFO] Processing sample: $sample_name"

  # Calculate column numbers for this sample's 3 columns
  start_col=$(( (sample_index - 1) * sample_cols + 1 ))
  end_col=$(( start_col + sample_cols - 1 ))

  # Extract this sample’s columns from mpileup.data
  cut -f${start_col}-${end_col} mpileup.data > "${sample_name}.data"

  # Combine position columns + sample data to create full mpileup for this sample
  paste mpileup.positions "${sample_name}.data" > "${sample_name}.mpileup"

  # The index plus one to move to next sample
  ((sample_index++))
done < "$bam_list"  # tells the while loop to read lines from bamlist.txt

### need 1 thead




### 4: Split each sample’s mpileup file by chromosome
chromosomes_list=(
  hcontortus_chr1_Celeg_TT_arrow_pilon
  hcontortus_chr2_Celeg_TT_arrow_pilon
  hcontortus_chr3_Celeg_TT_arrow_pilon
  hcontortus_chr4_Celeg_TT_arrow_pilon
  hcontortus_chr5_Celeg_TT_arrow_pilon
  hcontortus_chrX_Celeg_TT_arrow_pilon
  mitochondrion
)


for mpfile in *.mpileup; do
  (
  # Extract sample name by removing the .mpileup extension
  sample_name=${mpfile%.mpileup}
  echo "[INFO] Splitting sample $sample_name by chromosome..."

  for chr in "${chromosomes_list[@]}"; do
    # Extract lines that start with this chromosome name
    grep -P "^${chr}\t" "$mpfile" > "${sample_name}_${chr}.sub.mpileup"
  done
  ) &  # run all mpfile in parallel
done
wait   # run all mpfile in parallel
### need 15 samples * 1 thread = 15 cores





### 5: Clean up temporary files
rm *.data

echo "[Finished] All samples and chromosomes have be done"
```


### Submit the SLURM job

```bash
sbatch ...sh
```

## running time
5 hours



## input
```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-01_generate_mpileup

275G	all_samples_removed_s16_s17_fillmiss.mpileup

4.0K	bamlist_removed_s16_s17.txt

```


## ouput
there are 15 samples, and each sample has a whole-genome mpileup file and chromosome-specific mpileup files
```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-01_split_mpileup

# one sample's example


# one whole-genome mpileup file

31G	Sample_XJ-4126-9.mpileup

# 7 chromosome-specific mpileup files per sample 

4.9G	Sample_XJ-4126-9_hcontortus_chr1_Celeg_TT_arrow_pilon.sub.mpileup
5.5G	Sample_XJ-4126-9_hcontortus_chr2_Celeg_TT_arrow_pilon.sub.mpileup
5.0G	Sample_XJ-4126-9_hcontortus_chr3_Celeg_TT_arrow_pilon.sub.mpileup
5.8G	Sample_XJ-4126-9_hcontortus_chr4_Celeg_TT_arrow_pilon.sub.mpileup
5.1G	Sample_XJ-4126-9_hcontortus_chr5_Celeg_TT_arrow_pilon.sub.mpileup
4.9G	Sample_XJ-4126-9_hcontortus_chrX_Celeg_TT_arrow_pilon.sub.mpileup
13M	Sample_XJ-4126-9_mitochondrion.sub.mpileup
```