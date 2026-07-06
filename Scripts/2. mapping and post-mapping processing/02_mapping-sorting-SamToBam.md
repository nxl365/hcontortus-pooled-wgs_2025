## mapping, sorting, conversion from sam to bam


### Tools
* bwa
* samtool

This step uses **BWA-MEM** to align trimmed paired-end Illumina reads to the *Haemonchus contortus* reference genome. The SAM output from BWA-MEM is piped to **samtools sort**, which converts it to BAM format and sorts the alignments by genomic coordinates.


### SLURM bash script

Create the bash script:

```bash
module load nano
nano 5_mapping_sort_bam.sh
```

Script content:
1. **test: run one sample with 2 fastq files**
```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-***
#SBATCH -n 1                 # one task, 
#SBATCH --cpus-per-task=58   # each task, use 58 cores
#SBATCH -t 05:00:00
#SBATCH -p shared
#SBATCH -J map_sort_bam
#SBATCH --mail-type=All
#SBATCH --mail-user=nanxing.liu@slu.se   

module load bwa/0.7.18 samtools/1.20

refer_genome='.../raw_data/Hc_reference_genome/haemonchus_contortus.PRJEB506.WBPS19.genomic.fa'

trim_data_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/2_trim_fastp'

output_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/5_mapping_sort_bam'

mkdir -p "$output_dir"


	# r1 is: path/R1.fastq.gz file
for r1 in "$trim_data_dir"/*_fastp_trim_R1.fastq.gz;do
	# derive corresponding R2 file: replace suffix   ${var/oldstring/newstring}
	r2="${r1/_fastp_trim_R1.fastq.gz/_fastp_trim_R2.fastq.gz}"
	

	sample=$(basename "$r1") # extract entire r1 file name including .fastq.gz, remove path
	sample=${sample%%_fastp_trim_R1.fastq.gz}  # only extract sample name, remove the _fastp_trim_R1 
	
	echo "processing $sample"
	
	## bwa mem: mapping to get sam, then pipe to
	## samtools: sort, transfer sam to bam
	bwa mem \
	-t 50 \
	-R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina" \
	"$refer_genome" "$r1" "$r2"\
	2>"$output_dir/${sample}.bwa.log" | \
	
	samtools sort \
	-@ 8 \
	-o "$output_dir/${sample}.sort.bam" -
	
	# If only test one sample, break out of the loop here:
	break 
	
done


# each sample(R1+R2) will use 50+8 threads

echo "completed!"


```

2. **run all 17 samples finally**

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-***
#SBATCH -n 1                 # one task
#SBATCH --cpus-per-task=100  # 100 cores
#SBATCH -t 15:00:00
#SBATCH -p shared
#SBATCH -J map_sort_bam
#SBATCH --mail-type=All
#SBATCH --mail-user=nanxing.liu@slu.se   

module load bwa/0.7.18 samtools/1.20

refer_genome='.../raw_data/Hc_reference_genome/haemonchus_contortus.PRJEB506.WBPS19.genomic.fa'

trim_data_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/2_trim_fastp'

output_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/5_mapping_sort_bam'

mkdir -p "$output_dir"


	# r1 is: path/R1.fastq.gz file
for r1 in "$trim_data_dir"/*_fastp_trim_R1.fastq.gz;do
	# derive corresponding R2 file: replace suffix   ${var/oldstring/newstring}
	r2="${r1/_fastp_trim_R1.fastq.gz/_fastp_trim_R2.fastq.gz}"
	

	sample=$(basename "$r1") # extract entire r1 file name including .fastq.gz, remove path
	sample=${sample%%_fastp_trim_R1.fastq.gz}  # only extract sample name, remove the _fastp_trim_R1 
	
	echo "processing $sample"
	
	## bwa mem: mapping to get sam, then pipe to
	## samtools: sort, transfer sam to bam
	bwa mem \
	-t 80 \
	-R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina" \
	"$refer_genome" "$r1" "$r2"\
	2>"$output_dir/${sample}.bwa.log" | \
	
	samtools sort \
	-@ 20 \
	-o "$output_dir/${sample}.sort.bam" - 
	
	echo "finished $sample"
	
done

echo "All completed!"


# each sample(R1+R2) will use 80+20 threads, 
# in the for loop, each sample will be ran one by one, so each time we only need 100 cores
# in 'shared' partition, each node only have 128 cores

```

**code explanation**
```bash

bwa mem \  # software: BWA; command: mem
-t 2 \ # number of threads
-R '@RG\tID:S1\tSM:S1\tPL:illumina' \ 
# Add ´read group´ information. Some variant callers may require read group information later, such as GATK.
# This information is added to the SAM/BAM header, and same to the VCF header.
# Example: '@RG\tID:sample_ID\tSM:sample_name\tPL:sequencing_platform'
# \t means tab, which allows downstream tools to correctly read the sample information.

../01.ref/genome.fasta \ # reference genome; BWA index files should be in the same directory
../data/S1_1.fq.gz \ # read 1 of sample S1
../data/S1_2.fq.gz \ # read 2 of sample S1

2>S1.bwa.log | \ 
                 # Save log messages separate to a log file from the alignment output， otherwise logs would be written into the output file and lead to mistakes in variant calling process later
                 # | means pipe. It sends the SAM output to the next code.

samtools sort \ 
                # Sort the aligned reads by genomic coordinates for downstream analysis.
                # means: This reorders reads in the SAM/BAM file according to their mapped positions on the reference genome.
                # Without sorting, reads are generally kept in the original input/read order.
-@ 2 \ # number of threads
-m 1G \ # memory per thread
-o S1.sort.bam - # -o specifies the output file name
                 # The final "-" means: samtools sort reads input from the previous | pipe.
                 # you can also skip the "-" here, samtool should recognize input automatically



```

### Submit the SLURM job

```bash
sbatch 5_mapping_sort_bam.sh
```


### Running time:
1. test: run one sample with 2 fastq files: about 25 minutes.
2. run all 17 samples finally: 8 hours




### Input data

1. The decompressed *H. contortus* reference genome in .fa format， and it's index file:
```text
.../raw_data/Hc_reference_genome/

## refer genome
haemonchus_contortus.PRJEB506.WBPS19.genomic.fa


# bwa index output
|-- haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.amb
|-- haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.ann
|-- haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.bwt
|-- haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.pac
`-- haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.sa
``` 

2.  trimmed 17 samples. Each sample had two trimmed paired-end read files, resulting in 34 fastq.gz files in total.  
Each `.fastq.gz` file size: 5–7 GB.

Example input structure:
```text
2_trim_fastp/
│   ├── sample1_fastp_trim_R1.fastq.gz
│   └── sample1_fastp_trim_R2.fastq.gz
│   ├── sample2_fastp_trim_R1.fastq.gz
│   └── sample2_fastp_trim_R2.fastq.gz
...
```


### Output

for each sample: input 2 paired-end trimmed.fastq.gz files (R1 and R2), output one sorted bam file

```text
5_mapping_sort_bam/

Sample1.bwa.log 
Sample1.sort.bam 
```
size: each samples's sorted bam file 8-13 GB


**Code-- check the output bam file**
```bash
module load samtools
samtools view <...bam>|less -S
```