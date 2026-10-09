# FET for each SNP
To detect statistically significant differences in allele frequencies between farm populations while accounting for variation in sequencing depth, Fisher’s exact test (FET) was applied to each SNP for all pairwise comparisons


## Bash script

### create the bash script
```bash
mkdir -p '/cfs/klemming/projects/snic/naiss2025-23-66/scripts/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-06_poola_fet_by_chr_everySNP'

01_poola_fet_by_chr_everySNP.sh

```


### Script content
```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -J fet_snp
#SBATCH -N 1      # number of entire nodes, need more memory for this job
#SBATCH -t 24:00:00
#SBATCH -p main
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se

module load bioinfo-tools popoolation2/1201

indir='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-02_sync'
OUTDIR='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-06_poola_fet_by_chr_everySNP'



## filtering parameters
pool_size=2000
min_count=4
min_coverage=15 # use half of the lowest mean coverage 30
max_coverage="2%" # drops each sample’s top 2% deepest sites automatically


## prep and go in the output directory
mkdir -p $OUTDIR
cd $OUTDIR


# Chromosome list as an array
chromosomes_list=(
  hcontortus_chr1_Celeg_TT_arrow_pilon
  hcontortus_chr2_Celeg_TT_arrow_pilon
  hcontortus_chr3_Celeg_TT_arrow_pilon
  hcontortus_chr4_Celeg_TT_arrow_pilon
  hcontortus_chr5_Celeg_TT_arrow_pilon
  hcontortus_chrX_Celeg_TT_arrow_pilon
  mitochondrion
)



## 3. Fet calculations per chromosome
echo "[$(date)] Step3: fet per SNP per chromosome"
for chr in "${chromosomes_list[@]}"; do
  (
    perl /sw/bioinfo/popoolation2/1201/rackham/fisher-test.pl \
      --input "$indir/${chr}_rm_s16_s17.sync" \
      --output "${chr}_rm_s16_s17.fet" \
      --suppress-noninformative \
      --min-count $min_count \
      --min-coverage $min_coverage \
      --max-coverage "$max_coverage"
    echo "  [$(date)] Done FET $chr"
  ) &
done
wait
#  this script provides -log10(P-value) as output.
## will calculate fet per SNP, use windowsize 1bp by default


## 4. merge all chromosomes.fet to one.fet
echo "[$(date)] Step4: merge results"
cat *.fet > "all_samples_rm_s16_s17.by_chr.fet" # did not give order of chromosomes



## 5. extra point, just for IGV: convert the Fisher's exact test results into the .igv file format
perl /sw/bioinfo/popoolation2/1201/rackham/export/pwc2igv.pl \
--input "all_samples_rm_s16_s17.by_chr.fet" \
--output "all_samples_rm_s16_s17.by_chr.fet.igv"



echo "[$(date)] finished"



```

### Submit the SLURM job
```
sbatch ...sh
```

## Running time:
7 hours


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
Note: this script provides -log10(P-value) as output, we need to transfer that to P-value later

```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-06_poola_fet_by_chr_everySNP> du -sh *

27G	all_samples_rm_s16_s17.by_chr.fet
19G	all_samples_rm_s16_s17.by_chr.fet.igv

4.3G	hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17.fet
4.0K	hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17.fet.params

5.3G	hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17.fet
4.0K	hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17.fet.params

4.8G	hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17.fet
4.0K	hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17.fet.params

4.7G	hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17.fet
4.0K	hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17.fet.params

5.0G	hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17.fet
4.0K	hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17.fet.params
```