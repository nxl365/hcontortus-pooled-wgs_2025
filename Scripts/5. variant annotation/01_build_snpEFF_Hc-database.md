# build a database of H. contortus for snpEFF


## Tools
* snpEff


## SLURM bash script

### Create the bash script:

```bash
01_snpEFF_buildDB.sh
```

### Script content


```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-173
#SBATCH -n 1
#SBATCH -c 30
#SBATCH -t 00:30:00    # 30 mins
#SBATCH -p shared
#SBATCH -J eff_db
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=nanxing.liu@slu.se


module load bioinfo-tools
module load snpEff # version 5.2


## build your database folder, copy and edit config
mkdir -p '/cfs/klemming/projects/snic/naiss2025-23-66/software/snpEFF'

cd /cfs/klemming/projects/snic/naiss2025-23-66/software/snpEFF

cp /sw/bioinfo/snpEff/5.2/rackham/snpEff.config ./

echo "# Haemonchus contortus
HCON_2020.genome : HCON_2020
" >> snpEff.config


## build data folder
mkdir -p './data/HCON_2020'

# link: annotation, reference genome files
## the name must be `sequences.fa`, `genes.gff`
ln -s /cfs/klemming/projects/snic/naiss2025-23-66/raw_data/Hc_annotation/fixed-haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3 ./data/HCON_2020/genes.gff

ln -s /cfs/klemming/projects/snic/naiss2025-23-66/raw_data/Hc_reference_genome/haemonchus_contortus.PRJEB506.WBPS19.genomic.fa ./data/HCON_2020/sequences.fa


java -Xmx32g -jar /sw/bioinfo/snpEff/5.2/rackham/snpEff.jar build \
-c ./snpEff.config \
-gff3 -v HCON_2020 \
-noCheckCds -noCheckProtein # SnpEff to save the new database without checking it (i.e. neither CDS nor protein sequences are checked), otherwise I have errors



# I tried `java -jar`, and use partition 'shared' + c=1 core, it ran out of memory, 
# so now I use `30 nodes` with `java -Xmx32g -jar`


```


## Running time:
15 samples from 1.5 mins


## Input 

```bash
data/

fixed-haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3
haemonchus_contortus.PRJEB506.WBPS19.genomic.fa

snpEff.config
...
```


## Output

```bash
/cfs/klemming/projects/snic/naiss2025-23-66/software/snpEFF/data/HCON_2020

4.0K	genes.gff
8.0K	sequence.bin   # for MT chromosome, will know that if read the log.out 
7.5M	sequence.hcontortus_chr1_Celeg_TT_arrow_pilon.bin
7.3M	sequence.hcontortus_chr2_Celeg_TT_arrow_pilon.bin
6.9M	sequence.hcontortus_chr3_Celeg_TT_arrow_pilon.bin
7.6M	sequence.hcontortus_chr4_Celeg_TT_arrow_pilon.bin
7.2M	sequence.hcontortus_chr5_Celeg_TT_arrow_pilon.bin
7.7M	sequence.hcontortus_chrX_Celeg_TT_arrow_pilon.bin
4.0K	sequences.fa

19M	snpEffectPredictor.bin   # this should be the final database output we will use
...
```