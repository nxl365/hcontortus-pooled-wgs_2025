# Trimming

## Tools
* fastp

fastp was used for trimming adapters and removing low quality bases from paired-end illumina reads.

it can identify adaptor sequences automatically (can specify adapter sequences manually), and generates a quality report which can be compared with the **fastqc** reports to see whether the trimming worked.


## SLURM bash script

### Create the bash script:

```bash
module load nano
nano 2_trim_fastp.sh
```

### Script content:

```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-***
#SBATCH -n 1
#SBATCH --cpus-per-task=140
#SBATCH -t 01:00:00
#SBATCH -p shared
#SBATCH -J trim_illumina
#SBATCH --mail-type=All
#SBATCH --mail-user=nanxing.liu@slu.se     

module load fastp/0.24.0

data_dir='.../raw_data/illumina_HC_wgs_XJ4126'
output_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/2_trim_fastp'

mkdir -p "$output_dir"  # Create the directory if it doesn't exist; if exist, no error

# x is: path/sample_1 folder
for x in "$data_dir"/*; do
	sample_name=$(basename "$x")  # Extracts entire sample_1 folder name, remove path
	echo "Processing $sample_name..."
	
	fastp \
	-i "$x"/*_R1_001.fastq.gz \
	-I "$x"/*_R2_001.fastq.gz \
	-o "$output_dir/${sample_name}_fastp_trim_R1.fastq.gz" \
	-O "$output_dir/${sample_name}_fastp_trim_R2.fastq.gz" \
	--thread 8 -g -c -q 20 -l 20 \
	--json "$output_dir/${sample_name}_fastp_trim_report.json" \
	--html "$output_dir/${sample_name}_fastp_trim_report.html" &
done

wait   # & and wait: run all files in background in parallel, then exit
	     # every 2 fastq files(per sample) use 8 threads, so 8 * 17samples= 136 cores in total

echo "finished!"


## explain: 
# -i    input read1
# -I    input read2
# -o    output read1, set path
# -O    output read2, set path
# -g    trim polyG tails
# -c    base correction in overlapping paired data
# -q    Quality threshold per base required
# -l    Reads length threshold
# --html quality reports which we mainly check, set path and name: R1 and R2 will be in one report
# --json another report form, set path and name

			
```



### Submit the SLURM job

```bash
sbatch 2_trim_fastp.sh
```


## Running time:
Approximately 8 minutes.



## Input

The `$data_dir` contained 17 samples. Each sample had two paired-end read files, resulting in 34 fastq.gz files in total.  
Each `.fastq.gz` file size: 5–7 GB.

Example input structure:
```text
raw_data/illumina_HC_wgs_XJ4126/

├── sample1/
│   ├── sample1_R1_001.fastq.gz
│   └── sample1_R2_001.fastq.gz
...
```



## Output

fastp generated 4 output files for each sample.

Example output for one sample:

```text
2_trim_fastp/

sample1_fastp_trim_R1.fastq.gz  #trimmed R1 fastq.gz, size changed slightly after trimming
sample1_fastp_trim_R2.fastq.gz
sample1_fastp_trim_report.html   # QC report after trimming
sample1_fastp_trim_report.json
```

