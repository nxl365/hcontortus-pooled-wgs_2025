## Quality control before trimming

### Tools
* FastQC
* MultiQC

Raw Illumina paired-end reads were assessed using **FastQC**, and the individual FastQC reports were summarized using **MultiQC**.



### SLURM bash script

Create the bash script:

```bash
module load nano
nano 1_QC_before.sh
```

Script content:

```bash
#!/bin/bash -l     # -l: required to get the full environment with modules
                    # #!/bin/bash shebang line: tells the operating system which interpreter(bash shell) to use when executing the script.
 
#SBATCH -A naiss2025-22-***          # NAISS compute project number
#SBATCH -n 1                         # request one task (one main process)
#SBATCH --cpus-per-task=68           # gets 68 CPU cores (same to #SBATCH -c), in shared partition, the number of cores linked to RAM!so , no need to set --mem= 
#SBATCH -t 04:00:00                  # 4 hours will be given to this job
#SBATCH -p shared                    # partition name
#SBATCH -J QC_illumina               # your job's name
#SBATCH --mail-type=All              # Request an email when the job starts and ends.
#SBATCH --mail-user=nanxing.liu@slu.se     

 
module load PDC/23.12 multiqc/1.22.3-cpeGNU-23.12 # load multiqc
module load fastqc/0.12.1 #load fastqc


# Define input and output directories
data_dir='.../raw_data/illumina_HC_wgs_XJ4126'
output_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/1_QC_before'



# 1. run fastqc for each fastq.gz
    # Run each file in the background, Instead of waiting for one file to finish before starting the next, run all files simultaneously.

for x in "$data_dir"/*; do
	echo "Processing $x..."
	fastqc "$x"/*.fastq.gz -t 2 -o $output_dir &  

done   
# each run with 2 threads/cores, so 2*34 (34 fastq files) = 68 threads/cores in total

wait  
# Wait for all background jobs to finish before exiting



# 2. run multiqc to summarize all reports of fastqc
multiqc $output_dir -o $output_dir -n XJ4126_multiqc_report.html # XJ4126 is our project ID

echo "FastQC and MultiQC completed!"
```


### Submit the SLURM job

```bash
sbatch 1_QC_before.sh
```


### Running time:
Approximately 16 minutes.



### Input data

The `$data_dir` contained 17 samples. Each sample had two paired-end read files, resulting in 34 fastq.gz files in total.  
Each `.fastq.gz` file size: 5–7 GB.

Example input structure:
```text
data_dir/
├── sample1/
│   ├── sample1_R1_001.fastq.gz
│   └── sample1_R2_001.fastq.gz
...
```



### Output

1. FastQC generated one `.html` report and one `.zip` file for each FASTQ file.

Example FastQC output for one sample:

```text
sample1_R1_001_fastqc.html
sample1_R1_001_fastqc.zip
sample1_R2_001_fastqc.html
sample1_R2_001_fastqc.zip
```

2. MultiQC generated a merged report for all samples:

```text
XJ4126_multiqc_report.html
XJ4126_multiqc_report_data
```


