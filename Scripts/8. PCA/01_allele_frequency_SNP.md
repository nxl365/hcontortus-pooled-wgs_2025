# calculate the exact allele frequency differences of SNPs


## Bash script

### create the bash script
```bash
# code file

mkdir -p '/cfs/klemming/projects/snic/naiss2025-23-66/scripts/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-09_poola_AF_perSNP'

01_poola_AF_by_chr_perSNP.sh

```


### Script content

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -J AF_SNP
#SBATCH -N 1      # number of entire nodes, need more memory for this job
#SBATCH -t 24:00:00
#SBATCH -p main
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se

module load bioinfo-tools popoolation2/1201

indir='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-02_sync'
OUTDIR='/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-09_poola_AF_perSNP'



## filtering parameters
pool_size=2000
min_count=4
min_coverage=15 # use half of the lowest mean coverage 30
max_coverage="2%" # drops each sample’s top 2% deepest sites automatically
#window=10000
#step=5000

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



## 3. AF calculations per chromosome, per SNP
echo "[$(date)] Step3: per chromosome"
for chr in "${chromosomes_list[@]}"; do
  (
    perl /sw/bioinfo/popoolation2/1201/rackham/snp-frequency-diff.pl \
      --input "$indir/${chr}_rm_s16_s17.sync" \
      --output-prefix "${chr}_rm_s16_s17-AF" \
      --min-count $min_count \
      --min-coverage $min_coverage \
      --max-coverage $max_coverage
    echo "  [$(date)] Done $chr"
  ) &
done
wait

```

### Submit the SLURM job
```
sbatch ...sh
```


## Running time:
forgot, but probably couple of hours...



## Input
```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-02_sync> du -sh *

## only use the sync file for each chromosome as input
9.9G	hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17.sync
11G	hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17.sync
9.6G	hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17.sync
12G	hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17.sync
11G	hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17.sync
9.8G	hcontortus_chrX_Celeg_TT_arrow_pilon_rm_s16_s17.sync
3.2M	mitochondrion_rm_s16_s17.sync


62G	all_samples_removed_s16_s17_fillmiss.sync # did not use this as input, as it is so big
```



## Output
per chromosome has 3 output files:
will only use rc files in the next steps

```bash
analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-09_poola_AF_perSNP


4.0K	hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17-AF.params
2.3G	hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17-AF_pwc
908M	hcontortus_chr1_Celeg_TT_arrow_pilon_rm_s16_s17-AF_rc

4.0K	hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17-AF.params
2.6G	hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17-AF_pwc
1023M	hcontortus_chr2_Celeg_TT_arrow_pilon_rm_s16_s17-AF_rc

4.0K	hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17-AF.params
2.4G	hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17-AF_pwc
936M	hcontortus_chr3_Celeg_TT_arrow_pilon_rm_s16_s17-AF_rc

4.0K	hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17-AF.params
2.6G	hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17-AF_pwc
1023M	hcontortus_chr4_Celeg_TT_arrow_pilon_rm_s16_s17-AF_rc

4.0K	hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17-AF.params
2.6G	hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17-AF_pwc
1.1G	hcontortus_chr5_Celeg_TT_arrow_pilon_rm_s16_s17-AF_rc

4.0K	hcontortus_chrX_Celeg_TT_arrow_pilon_rm_s16_s17-AF.params
1.1G	hcontortus_chrX_Celeg_TT_arrow_pilon_rm_s16_s17-AF_pwc
397M	hcontortus_chrX_Celeg_TT_arrow_pilon_rm_s16_s17-AF_rc

4.0K	mitochondrion_rm_s16_s17-AF.params
2.4M	mitochondrion_rm_s16_s17-AF_pwc
1004K	mitochondrion_rm_s16_s17-AF_rc
```

### output explanation

1. two important output files having two different extensions:  

* _rc: this file contains the major and minor alleles for every SNP in a concise format

* _pwc: this file contains the differences in allele frequencies for every pairwise comparision of the populations present in the synchronized file
For details see the man pages of the script



2. For details see the man pages of the script:

```bash
## code
man /sw/bioinfo/popoolation2/1201/rackham/snp-frequency-diff.pl
```

```bash

# output rc: exact AF of SNP

===head2 OUTPUT <output-prefix>_rc

 ##chr    pos      rc       allele_count    allele_states   
 #deletion_sum    snp_type         major_alleles(maa)       
 #minor_alleles(mia)       
 #maa_1    maa_2    maa_3    maa_4    maa_5    
 #mia_1    mia_2   mia_3   mia_4   mia_5
 gi|42519920|ref|NC_002978.6      3530     A        2         A/G
0        pop      AAAAA    GGNNG    54/55   76/78   35/35   27/27
17/18   1/55    2/78    0/35    0/27    1/18

 gi|42519920|ref|NC_002978.6      5450     T        2         T/A
0        pop      TTTTT    ANANN    51/54   96/96   9/10    33/33
23/23   3/54    0/96    1/10    0/33    0/23


 col1: reference chromosome (contig)
 col2: referenced position
 col3: reference character
 
 col4: number of alleles found in all populations
 
 col5: allele characters in all populations (sorted by counts in all populations)
 
 col6:  sum  of  deletions in all populations (should be zero, if not the postion may not be reliable)
 
 col7: SNP type: [pop, rc, rc|pop]; 
	 pop.. a SNP within or between the  populations; 
	 rc.. a SNP between the reference sequence character and the 
	      consensus of at least one populaton; 
	 rc|pop..both
 
 col8: most frequent allele in all populations [12345..order]
 
 col9: second most frequent allele in all populations [12345..order]
 
 col10 - col9+n: frequencies of the most frequent allele  
 (major)in the form "allele-count/coverage"
 
 col10+n - col9+2n: frequencies of the second most frequent 
 allele (minor) in the form "allele-count/coverage"
 
 
```

```bash
# output pwc: SNP AF difference 

====head2 OUTPUT <output-prefix>_pwc

 ##chr   pos     rc       allele_count     allele_states    
 #deletion_sum      snp_type          most_variable_allele     
 #diff:1-2     diff:1-3  diff:1-4  diff:1-5         diff:2-3      diff:2-4      diff:2-5   diff:3-4    diff:3-5    diff:4-5
 
 
 
gi|42519920|ref|NC_002978.6      3530     A        2         A/G
0        pop      A        0.007    0.018   0.018   0.037   0.026
0.026   0.030   0.000   0.056   0.056


 gi|42519920|ref|NC_002978.6      5450     T        2         T/A
0        pop      A        0.056    0.044   0.056   0.056   0.100
0.000   0.000   0.100   0.100   0.000

 col1: reference chromosome (contig)
 col2: referenced position
 col3: reference character
 col4: number of alleles found in all populations
 
 col5: 
	 allele characters in all populations (sorted by counts  in
all populations)

 col6:  
	 sum  of  deletions in all populations (should be zero,if not the postion may not be reliable)
	 
 col7: 
 SNP type: [pop, rc, rc|pop]; 
 pop.. a SNP within or between the  populations; 
	rc.. a SNP between the reference sequence character and the consensus of at least one populaton; 
	rc|pop..both

 col8: 
	 allele showing the most variability in terms of frequency
between  all pairwise comparisions; 
	for biallelic SNPs one allele is picked at random
	
 diff:1-2 frequency difference for allele <col8> between  population-1 and population-2
 
 diff:2-4  frequency difference for allele <col8> between population-2 and population-4
 ...


```