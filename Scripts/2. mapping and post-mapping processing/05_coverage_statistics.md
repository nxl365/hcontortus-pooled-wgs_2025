# Coverage statistics

## Tools
* samtools

`samtools coverage` was used to calculate coverage statistics from each sorted BAM file. The main output statistics include the number of reads, number of covered bases, coverage percentage, mean sequencing depth, mean base quality, and mean mapping quality.


## SLURM bash script

Create the bash script:

```bash
module load nano
nano 8_coverage_after_mapping_statistics.sh
```

Script content:
```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-***
#SBATCH -n 1                 # one task
#SBATCH --cpus-per-task=17  
#SBATCH -t 01:00:00
#SBATCH -p shared
#SBATCH -J coverage_rate
#SBATCH --mail-type=All,TIME_LIMIT_90   # get a warning when I reach 90% of the time limit I put
#SBATCH --mail-user=nanxing.liu@slu.se  


input_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/5_mapping_sort_bam'

output_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/8_coverage_after_mapping_statistics'

mkdir -p $output_dir


module load samtools/1.20

for x in "$input_dir"/*.sort.bam;do
	sample=$(basename "$x")
	sample_name=${sample%%.sort.bam}
	
	echo "processing $sample_name"
		
	samtools coverage $x > "$output_dir/${sample_name}.sort.bam.coverage" &
	

	
done

wait

echo "all completed!"



# -@ 1 should be same to not using -@, since 1 thread is the default behavior
# run al 17 samples in parallel in background, 1*17=17 threads/cores in total
```


Submit the SLURM job

```bash
sbatch 8_coverage_after_mapping_statistics.sh
```

## Running time:
~ 10 mins


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

8_coverage_after_mapping_statistics/

Sample1.sort.bam.coverage

```

### Output explanation

* source: [samtool-coverage manual](https://www.htslib.org/doc/samtools-coverage.html)  

* One output example :
```bash
less Sample.sort.bam.coverage

#rname      startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq
hcontortus_chr4_Celeg_TT_arrow_pilon    1       51826579        28100188        50944075        98.2972 71.418  38.7    49.3
hcontortus_chr5_Celeg_TT_arrow_pilon    1       48868368        26401769        47136636        96.4563 70.4964 38.7    47.1
hcontortus_chr2_Celeg_TT_arrow_pilon    1       47382676        27098637        46477783        98.0902 77.2694 38.7    48.1
hcontortus_chrX_Celeg_TT_arrow_pilon    1       46012677        26480002        43866556        95.3358 79.0603 38.6    43.5
hcontortus_chr1_Celeg_TT_arrow_pilon    1       45774638        24003490        44283137        96.7416 68.2165 38.7    48.7
hcontortus_chr3_Celeg_TT_arrow_pilon    1       43560352        24585081        42549624        97.6797 74.5126 38.7    47.3
mitochondrion   1       14018   84139   14018   100     862.296 39.1    59.1
```

* main statisitcs we should focus on:
    - **coverage**: the percentage of bases that are covered by at least one read (0–100%). it reflects how much of the reference genome is covered by the sequencing data.  

    - **meandepth**: the average number of reads covering each reference base. it represents the average sequencing depth, usually reported as “×” depth.

* Additional:
    - **meanbaseq**: the average Phred base‐quality score across all covered bases, reflecting sequencing accuracy.
    - **meanmapq**: the average Phred mapping‐quality score of the reads used for coverage, reflecting alignment confidence.



# extraction of coverage statistics

## python script

Create the python script:

```python
extract_flagstat_lines.py
```

Script content:

```python
#!/usr/bin/env python3
import os
import glob
import csv

output_csv = os.path.join("./", "./combined_coverage.csv")


# Glob all coverage files in the directory
coverage_files = glob.glob(os.path.join("./", "*.coverage"))

# Define the CSV header
#    The samtools coverage output typically includes columns like:
#      #rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq
#    We'll add "sample" as the first column.
header = [
    "sample",
    "rname",
    "startpos",
    "endpos",
    "numreads",
    "covbases",
    "coverage",
    "meandepth",
    "meanbaseq",
    "meanmapq"
]

with open(output_csv, mode="w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    # Write the header
    writer.writerow(header)

    # Loop through each coverage file
    for cov_file in coverage_files:
        # Extract a "sample name" from the file name, eg,from "Sample_XJ-4126-1.sort.bam.coverage"
        base_name = os.path.basename(cov_file)
        sample_name = base_name.replace(".sort.bam.coverage", "")
        
        with open(cov_file, "r") as infile:
            for line in infile:
                # Skip lines that are comments or empty
                if line.startswith("#") or not line.strip():  # line.strip()--- removes leading and trailing whitespace
                    continue
                
                # Parse the columns in the coverage file
                cols = line.strip().split()
                # Expected order: rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq
                # Insert the sample name as the first column
                row = [sample_name] + cols
                # Write to the combined CSV
                writer.writerow(row)

print(f"Combined coverage written to: {output_csv}")
```

### Input
17 samples;
for each sample: one output file

```text

8_coverage_after_mapping_statistics/

Sample1.sort.bam.coverage

```

### Output
one output file containing coverage statistics information for all 17 samples: 
```text
combined_coverage.csv 
```