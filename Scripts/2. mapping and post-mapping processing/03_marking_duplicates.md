## marking duplicates

### Tools
* picard

Picard `MarkDuplicates` was selected here to identify and mark duplicate reads in the coordinate-sorted BAM files. Duplicate reads were not removed. The output BAM files were indexed during this step.

* compare Samtool and Picard:
`samtools markdup`: faster; more suitable for SNP and indel workflows.
`Picard MarkDuplicates`: slower; widely used in GATK-style workflows; provides detailed duplicate metrics and better compatibility with possible future SV analyses.


### SLURM bash script

Create the bash script:
```bash
module load nano
nano 6_mark_duplicate.sh
```

Script content:
```bash
#!/bin/bash -l
#SBATCH -A naiss2025-22-***
#SBATCH -n 1                 # one task
#SBATCH --cpus-per-task=50   # Each node has 128 cores, which are allocated to me automatically, so not need to specify cores manually here actually, I got 128 anyway
#SBATCH -t 20:00:00
#SBATCH -p main        # `main` partition: I was allocated a entire node with 128 cores at least
#SBATCH -J mark_duplicate
#SBATCH --mail-type=All
#SBATCH --mail-user=nanxing.liu@slu.se   

module load picard/3.3.0

mapped_sorted_data_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/5_mapping_sort_bam'

output_dir='.../analysis/variant_explore/illumina_HC_wgs_XJ4126/6_mark_duplicate'

mkdir -p "$output_dir"


for x in "$mapped_sorted_data_dir"/*.sort.bam; do
	sample_file_name=$(basename "$x")
	sample_name=${sample_file_name%%.sort.bam}
	
	echo "processing $sample_name"

	java -Xmx4g -XX:ParallelGCThreads=2 \
	-jar $PICARD \
	MarkDuplicates I="$x" \
	O="$output_dir/${sample_name}.sort.markdup.bam" \
	CREATE_INDEX=true \
	REMOVE_DUPLICATES=False \
	M="$output_dir/${sample_name}.marked_dup_metrics.txt" &

done

wait

echo "all completed!"

# $x include the path/file
# M= just write dup reads into metrics, we won't use it!
# REMOVE_DUPLICATES=False, marking instead of removing dup

# run each loop/sample in background, so run all samples in parallel, each loop use 2 threads, so 17sample * 2 = 34 threads / CUP cores we need in total

# in PDC Partition ´shared´: 1 core = 2 GB, so using -Xmx4g -XX:ParallelGCThreads=2 should be reasonable
# but I used `main` Partition at the end, as my job was skilled in shared partition,


```

### Submit the SLURM job

```bash
sbatch 6_mark_duplicate.sh
```


### Running time:
~ 54 mins



### Input data

17 samples in total; 
for each sample: one sorted bam file (size: 8-13 GB)

```text
5_mapping_sort_bam/
    
    Sample1.sort.bam 
    Sample2.sort.bam ...

```



### Output
for each sample: 3 output files  

```text

sample1.sort.markdup.bam      #bam, after marking/removing duplicates

sample1.sort.markdup.bai     # index file,  from CREATE_INDEX=true , GATK will require that

sample1.marked_dup_metrics.txt   # just write dup reads into metrics, but we won't use it!

```
* It’s good to have the .bam and .bai named to the same “prefix”!


* size: each sample's bam file after marking duplicates  8.6-14 GB (become a bit bigger )

