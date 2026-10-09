# Pairwise Per-SNP FST
the process is same to sliding-window FST, but window size is 1 and step size is 1

## Bash script

### Create the bash script:
```bash
mkdir -p '/cfs/klemming/projects/snic/naiss2025-23-66/scripts/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-04_poola_fst_by_chr_everySNP'

01_poola_fst_by_chr_everySNP.sh

```

### Script content

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -J fst_snp
#SBATCH -N 1      # number of entire nodes, need more memory for this job
#SBATCH -t 5-00:00:00  # days
#SBATCH -p long
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se

module load bioinfo-tools popoolation2/1201

indir='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-02_sync'

outdir_base='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-04_poola_fst_by_chr_everySNP'


# parameters: try different window and step sizes
window=1
step=1


## filtering parameters for fst
pool_size=2000
min_count=4
min_coverage=15 # use half of the lowest mean coverage 30
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
rm -f all_samples_removed_s16_s17_fillmiss.by_chr.fst # force remove any old .fst file, in case add new info to it; even no old .fst exist, no error come out
cat *.fst > "all_samples_rm_s16_s17_${window}_window_${step}_step.by_chr.fst" # did not give order of chromosomes



echo "[$(date)] finished"

```

### Submit the SLURM job
```
sbatch ...sh
```

## Running time:
10.5 h



## Input

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




## Output

```bash

/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-04_poola_fst_by_chr_everySNP/fst_1_window_1_step

##
27G	all_samples_rm_s16_s17_1_window_1_step.by_chr.fst # a big file after merging all chromosomes fst files below

4.3G	hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17.fst
4.0K	hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params

5.3G	hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17.fst
4.0K	hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params

4.8G	hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17.fst
4.0K	hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params

4.7G	hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17.fst
4.0K	hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params

5.0G	hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17.fst
4.0K	hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params

2.2G	hcontortus_chrX_Celeg_TT_arrow_pilon_rm_s16_s17.fst
4.0K	hcontortus_chrX_Celeg_TT_arrow_pilon_rm_s16_s17.fst.params
```