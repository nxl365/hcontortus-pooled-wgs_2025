# mapping rate statistics


## Tools
* samtools

This step uses `samtools flagstat` to summarize basic mapping statistics for each coordinate-sorted BAM file.


## SLURM bash script

### Create the bash script:

```bash
module load nano
nano 7_mapping_rate_statistics.sh
```


### Script content:
```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-***
#SBATCH -n 1                 # one task
#SBATCH --cpus-per-task=34  
#SBATCH -t 01:00:00
#SBATCH -p shared
#SBATCH -J map_rate
#SBATCH --mail-type=All,TIME_LIMIT_90   # get a warning when I reach 90% of the time limit I put
#SBATCH --mail-user=nanxing.liu@slu.se  


input_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/5_mapping_sort_bam'

output_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/7_mapping_rate_statistics'

mkdir -p $output_dir


module load samtools/1.20

for x in "$input_dir"/*.sort.bam;do
	sample=$(basename "$x")
	sample_name=${sample%%.sort.bam}
	
	echo "processing $sample_name"
		
	samtools flagstat -@ 2 $x > "$output_dir/${sample_name}.sort.bam.flagstat" &
	

	
done

wait


echo "all completed!"

# -@ 1 should be same to not using -@, since 1 thread is the default behavior

# run all 17 samples in parallel in background, 2*17=34 threads/cores in total
```


### Submit the SLURM job

```bash
sbatch 7_mapping_rate_statistics.sh
```

## Running time:
~ 2 mins


## Input
* input is bam file after mapping and sorting, not after marking duplicates.  

17 samples;
for each sample: one sorted bam file (size: 8-13 GB)

```text
5_mapping_sort_bam/

Sample1.sort.bam 
...
```

## Output 

17 samples;
for each sample: one output file

```text

7_mapping_rate_statistics/

Sample1.sort.bam.flagstat

```

### Output explanation

**Mapped %** is the broadest mapping statistic. it indicates the percentage of all reads aligned to the reference genome, including primary, secondary, and supplementary alignments.

**Primary mapped %** is a stricter and commonly used mapping statistic.It indicates the percentage of primary alignments that were mapped to the reference genome, excluding secondary and supplementary alignments.

**Properly Paired %** is a even more stricter statistic.It indicates the percentage of paired-end reads, where both forward and reverse reads were mapped in the expected orientation and distance.



# Extraction of mapping rate statistics


The extraction mainly focuses on: primary mapped, properly paired



## python script

Create the python script:

```python
extract_flagstat_lines.py
```

Script content:

```python
import glob

flagstat_files = glob.glob("Sample_*.sort.bam.flagstat")

output_file = "extracted_mapping_rate_lines.txt"

with open(output_file, "w") as out:
    for file in sorted(flagstat_files):
        matching_lines = []
        with open(file, "r") as f:
            for line in f:
                if "primary mapped" in line or "properly paired" in line:
                    matching_lines.append(line.strip())
        if matching_lines:
            out.write(f"{file}\n")
            for l in matching_lines:
                out.write(f"{l}\n")
            out.write("\n")  # Blank line for separation
                    
print(f"Extracted lines have been written to {output_file}")



```

## Input
all 17 samples' `.flagstat` files

```text
Sample*.sort.bam.flagstat
```

## Output
one output file containing information for all 17 samples: primary mapped %, properly paired %

```text
extracted_mapping_rate_lines.txt 
```