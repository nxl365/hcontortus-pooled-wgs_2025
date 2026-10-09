# Pairwise sliding-window FST

* Calculate FST separately for 7 chromosomes and merge at the end

* Pairwise sliding-window FST analysis was performed using three window sizes (1 kb, 10 kb, 100 kb) with corresponding step sizes (0.5 kb, 5 kb, 50 kb).


## Bash script

1. parameterized script:

```bash

mkdir -p '/cfs/klemming/projects/snic/naiss2025-23-66/scripts/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr'

01_poola_fst_by_chr.sh

```

```bash
########## for: window=1000, step=500;  window=10000, step=5000 

#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -J fst
#SBATCH -n 1
#SBATCH -c 7
#SBATCH -t 3-00:00:00   # 3 days
#SBATCH -p shared
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se


########## for: window=100000, step=50000, need more RAM
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -J fst_100k
#SBATCH -N 1      # number of entire nodes
#SBATCH -t 24:00:00   # main partition has maximum 24h
#SBATCH -p main
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se


module load bioinfo-tools popoolation2/1201

indir='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-02_sync'

outdir_base='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr'


# parameters (from sbatch command line): try different window and step sizes
window=$1
step=$2


## filtering parameters for fst
pool_size=2000
min_count=4
min_coverage=15 # use half of the lowest mean coverage ~30 among all samples (range: 32.8–92.0X): refer to "Analysis_results/2. Mapping and post-mapping processing/05_coverage_statistics/combined_coverage.csv"
max_coverage="2%" # drops each sample’s top 2% deepest sites automatically


## prep and go in the output directory
OUTDIR="${outdir_base}/fst_${window}_window_${step}_step"
mkdir -p $OUTDIR
cd $OUTDIR


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



## 3. Fst calculations per chromosome
echo "[$(date)] Step3: fst-sliding per chromosome"
for chr in "${chromosomes_list[@]}"; do
  (
    perl /sw/bioinfo/popoolation2/1201/rackham/fst-sliding.pl \
      --input "$indir/${chr}_rm_s16_s17.sync" \
      --output "${chr}_rm_s16_s17.fst" \
      --suppress-noninformative \
      --min-count $min_count \
      --min-coverage $min_coverage \
      --max-coverage "$max_coverage" \
      --window-size $window \
      --step-size $step \
      --pool-size $pool_size
    echo "  [$(date)] Done FST $chr"
  ) &
done
wait

## single-threaded Perl script
## do all 7 chromosomes at the same time, so will need 7 threads



## 4. merge all chromosomes.fst to one.fst
echo "[$(date)] Step4: merge results"
cat *.fst > "all_samples_rm_s16_s17_${window}_window_${step}_step.by_chr.fst" # did not give order of chromosomes



echo "[$(date)] finished"

```

2. wrapper script:

```bash
02_poola_fst_by_chr-run_all_windows.sh
```

```bash
#!/bin/bash -l

#window=1000 , step=500 
#window=10000 , step=5000  
#window=100000 , step=50000 

windows=(
  "1000"
  "10000"
  "100000"
)

for window in "${windows[@]}"; do
	step=$((window / 2))
    echo "Calculating Fst : window=$window , step=$step  "
    sbatch 01_poola_fst_by_chr.sh $window $step
done

```

### submit the slurm job 
Run the wrapper to run 3 windows sizes
```bash
bash 02_poola_fst_by_chr-run_all_windows.sh
```

## running time

* window=1000 , step=500: 7h
* window=10000 , step=5000: 6.5h
* window=100000 , step=50000: 5h


## input
7 chromosome-specific sync files

```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-02_sync

9.9G	hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17.sync
11G	hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17.sync
9.6G	hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17.sync
12G	hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17.sync
11G	hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17.sync
9.8G	hcontortus_chrX_Celeg_TT_arrow_pilon_rm_s16_s17.sync
3.2M	mitochondrion_rm_s16_s17.sync
```


## output

* fst files

```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr


# there are 3 folders for 3 window sizes: 

fst_100000_window_50000_step
fst_10000_window_5000_step
fst_1000_window_500_step


# in each window-size folder: all files are filled miss...

|-- all_samples_rm_s16_s17_1000_window_500_step.by_chr.fst # a big file after merging all chromosomes fst files below

|-- hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17.fst
|-- hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params

|-- hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17.fst
|-- hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params

|-- hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17.fst
|-- hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params

|-- hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17.fst
|-- hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params

|-- hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17.fst
|-- hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params

|-- hcontortus_chrX_Celeg_TT_arrow_pilon_rm_s16_s17.fst
|-- hcontortus_chrX_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params

|-- mitochondrion_rm_s16_s17.fst
`-- mitochondrion_rm_s16_s17.fst.params

```


* log

This log shows the exact maximum coverage cutoff 2% `--max-coverage 2%` for each sample (15 samples in total) on each chromosome

```bash
Computing maximum coverage cutoffs from the empirical distributions of coverages
Finished identifying maximum coverages
Result: '--max-coverage 2%' is equivalent to '--max-coverage 480,486,476,484,484,487,488,487,489,485,479,486,488,489,485'
  [Thu Jul  3 16:06:59 CEST 2025] Done FST mitochondrion
  
Computing maximum coverage cutoffs from the empirical distributions of coverages
Finished identifying maximum coverages
Result: '--max-coverage 2%' is equivalent to '--max-coverage 76,85,59,58,72,74,69,55,68,58,39,41,67,78,52'
  [Thu Jul  3 19:54:03 CEST 2025] Done FST hcontortus_chrX_Celeg_TT_arrow_pilon

Computing maximum coverage cutoffs from the empirical distributions of coverages
Finished identifying maximum coverages
Result: '--max-coverage 2%' is equivalent to '--max-coverage 69,75,59,57,68,72,71,58,69,54,35,42,67,80,52'
  [Thu Jul  3 20:25:44 CEST 2025] Done FST hcontortus_chr3_Celeg_TT_arrow_pilon

Computing maximum coverage cutoffs from the empirical distributions of coverages
Finished identifying maximum coverages
Result: '--max-coverage 2%' is equivalent to '--max-coverage 71,77,60,60,69,73,71,58,69,55,36,43,67,80,52'
  [Thu Jul  3 20:34:15 CEST 2025] Done FST hcontortus_chr1_Celeg_TT_arrow_pilon

Computing maximum coverage cutoffs from the empirical distributions of coverages
Finished identifying maximum coverages
Result: '--max-coverage 2%' is equivalent to '--max-coverage 71,77,60,60,70,74,71,58,70,55,36,42,67,82,52'
  [Thu Jul  3 20:58:25 CEST 2025] Done FST hcontortus_chr4_Celeg_TT_arrow_pilon

Computing maximum coverage cutoffs from the empirical distributions of coverages
Finished identifying maximum coverages
Result: '--max-coverage 2%' is equivalent to '--max-coverage 72,80,62,60,71,75,73,59,72,57,37,43,68,83,54'
  [Thu Jul  3 21:03:51 CEST 2025] Done FST hcontortus_chr5_Celeg_TT_arrow_pilon

Computing maximum coverage cutoffs from the empirical distributions of coverages
Finished identifying maximum coverages
Result: '--max-coverage 2%' is equivalent to '--max-coverage 69,74,58,56,68,72,69,56,67,54,35,42,66,80,51'
  [Thu Jul  3 21:03:51 CEST 2025] Done FST hcontortus_chr2_Celeg_TT_arrow_pilon

```