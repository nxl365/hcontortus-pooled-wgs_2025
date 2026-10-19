# Run the tool NPStat
NPStat will calculate nucleotide diversity (Pi), Watterson’s theta, Tajima’s D per sample per chromosome


## Bash script

### Create the bash script:

```bash
/cfs/klemming/projects/snic/naiss2025-23-66/scripts/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-02_run_npstats 

# script
07-02_run_npstats.sh

```

### Script content

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -J npstats
#SBATCH -n 1
#SBATCH -c 120  
#SBATCH -t 5:00:00
#SBATCH -p shared
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se


input_dir="/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-01_split_mpileup"

out_dir="/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-02_run_npstats"

mkdir -p $out_dir



# Move into output directory
cd "$out_dir"


# debug for : npstat can't read "::env(HOME)
export HOME="$out_dir"

module load bioinfo-tools NPStat/1


for sample_chr in $input_dir/*.sub.mpileup;do
	(
	npstat -n 4000 -l 10000 -maxcov 90 -mincov 15 -minqual 20 -nolowfreq 4 $sample_chr
	) &
done
wait
## 15 samples * 7 chromosome = 105 files in parallel, so need 105 cores


# the ouput files.stats will be in $input_dir automatically

## move .stats to out_dir
mv $input_dir/*.stats $out_dir




```

### NPStat parameters explanation
```bash
npstat
   -n samplesize : haploid sample size
   -l windowlength : window length 
   -mincov minimum_coverage : filter on minimum coverage (default 4; means= per site read depth = 4)
   -maxcov maximum_coverage : filter on maximum coverage (default 100)
   -minqual minimum_base_quality : filter on base quality (default 10)
   -nolowfreq m : filter on minimum allele count mac>m

```
```bash
# NPStat parameters were set to match the FST calculation
# NPStat  ~   FST
-l 10000 ~  10kb windowsize
## there is no step set，we use non-overlapping windows

-mincov 15 ~  min_coverage=15

-nolowfreq 4 ~ min_count=4

-maxcov 90 ~ since max_coverage="2%" was used for FST, based on a `log
shows the exact maximum coverage cutoff 2% `--max-coverage 2%` for each sample`, we choose 90 as the near cutoff

-n: haploid sample size is 4000 ~  since pool_size=2000 for Fst, it is deploid 
## Reference: https://github.com/stephenrdoyle/hcontortus_xqtl/blob/master/03_code/hcontortus_xqtl.workbook.mapping_variants.md; pool_size=200, -n 400

```



### Submit the SLURM job

```bash
sbatch ...sh
```


## running time
30 mins



## input
there are 15 samples in total, e.g. one sample
```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-01_split_mpileup


# 7 chromosome-specific mpileup files per sample 

4.9G	Sample_XJ-4126-9_hcontortus_chr1_Celeg_TT_arrow_pilon.sub.mpileup
5.5G	Sample_XJ-4126-9_hcontortus_chr2_Celeg_TT_arrow_pilon.sub.mpileup
5.0G	Sample_XJ-4126-9_hcontortus_chr3_Celeg_TT_arrow_pilon.sub.mpileup
5.8G	Sample_XJ-4126-9_hcontortus_chr4_Celeg_TT_arrow_pilon.sub.mpileup
5.1G	Sample_XJ-4126-9_hcontortus_chr5_Celeg_TT_arrow_pilon.sub.mpileup
4.9G	Sample_XJ-4126-9_hcontortus_chrX_Celeg_TT_arrow_pilon.sub.mpileup
13M	Sample_XJ-4126-9_mitochondrion.sub.mpileup

```



## ouput
there are 15 samples and each sample has 7 chromosome-specific stats file:

some stats files as examples were uploaded to `Analysis_results`

```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-02_run_npstats

# one sample example
444K	Sample_XJ-4126-9_hcontortus_chr1_Celeg_TT_arrow_pilon.sub.mpileup.stats
460K	Sample_XJ-4126-9_hcontortus_chr2_Celeg_TT_arrow_pilon.sub.mpileup.stats
424K	Sample_XJ-4126-9_hcontortus_chr3_Celeg_TT_arrow_pilon.sub.mpileup.stats
504K	Sample_XJ-4126-9_hcontortus_chr4_Celeg_TT_arrow_pilon.sub.mpileup.stats
468K	Sample_XJ-4126-9_hcontortus_chr5_Celeg_TT_arrow_pilon.sub.mpileup.stats
440K	Sample_XJ-4126-9_hcontortus_chrX_Celeg_TT_arrow_pilon.sub.mpileup.stats
4.0K	Sample_XJ-4126-9_mitochondrion.sub.mpileup.stats


```

### output explanation
```bash
18 columns

window	length	length_outgroup	read_depth	S	Watterson	Pi	Tajima_D	var_S	var_Watterson	unnorm_FayWu_H	FayWu_H	div	nonsyn_pol	syn_pol	nonsyn_div	syn_div	alpha

Each row of the files represents the statistics for a single window. 
For the "_stat" files, the columns represent:

1) window number

2) bases covered by reads (after filtering by the config you set)
(npstat use this to calculate Pi per window)

3) bases covered by reads and with known ancestral allele (after filtering) 
	
4) average read depth
5) number of segregating sites
6) Watterson estimator of theta
7) Tajima (Pi) estimator of theta
8) unnormalized Tajima's D
...
```


