## Quality control after trimming
same process to QC before trimming

### Tools
* FastQC
* MultiQC

Trimmed Illumina paired-end reads were assessed using **FastQC**, and the individual FastQC reports were summarized using **MultiQC**.



### SLURM bash script

Create the bash script:

```bash
module load nano
nano 3_QC_after_trim_fastp.sh
```

Script content:

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-***
#SBATCH -n 1
#SBATCH --cpus-per-task=100
#SBATCH -t 01:00:00
#SBATCH -p shared
#SBATCH -J QCaftertrim_illumina
#SBATCH --mail-type=All
#SBATCH --mail-user=nanxing.liu@slu.se      

 
module load PDC/23.12 multiqc/1.22.3-cpeGNU-23.12 # load multiqc
module load fastqc/0.12.1 #load fastqc


data_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/2_trim_fastp'
output_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/3_QC_after_trim_fastp'

mkdir -p "$output_dir"


# 1. run fastqc for each fastq.gz

fastqc "$data_dir"/*.fastq.gz -t 100 -o $output_dir
# run all fastq files in parallel, can use 100 threads


# 2. run multiqc to summarize all reports of fastqc
multiqc $output_dir -o $output_dir -n XJ4126_after_trim_fastp_multiqc_report.html


echo "FastQC and MultiQC completed!"
```


### Submit the SLURM job

```bash
sbatch 3_QC_after_trim_fastp.sh
```


### Running time:
Approximately 16 minutes.



### Input data

The `$data_dir` contained trimmed 17 samples. Each sample had two trimmed paired-end read files, resulting in 34 fastq.gz files in total.  
Each `.fastq.gz` file size: 5–7 GB.

Example input structure:
```text
data_dir/
│   ├── sample1_fastp_trim_R1.fastq.gz
│   └── sample1_fastp_trim_R2.fastq.gz
│   ├── sample2_fastp_trim_R1.fastq.gz
│   └── sample2_fastp_trim_R2.fastq.gz
...
```



### Output

1. FastQC generated one `.html` report and one `.zip` file for each trimmed FASTQ file.

Example FastQC output for one trimmed sample:

```text
sample1_fastp_trim_R1_fastqc.html
sample1_fastp_trim_R1_fastqc.zip
sample1_fastp_trim_R2_fastqc.html
sample1_fastp_trim_R2_fastqc.zip
```

2. MultiQC generated a merged report for all trimmed samples:

```text
XJ4126_after_trim_fastp_multiqc_report.html   
XJ4126_after_trim_fastp_multiqc_report_data
```

