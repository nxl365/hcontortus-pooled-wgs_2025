# preprocess the stats files from NPStat


## 1. Preprocess chromosome-specific stats files for each sample

1. add genomic positions
2. remove missing values based on `columns_to_keep`
3. also save summary of missing values per sample as tsv 

### Python script
```python

import os
import glob
import pandas as pd


input_dir = "../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-02_run_npstats"

output_dir = os.path.join(input_dir, "01_preprocessing_npstat_data")
os.makedirs(output_dir, exist_ok=True)

window_size = 10000
columns_to_keep = [
    "window", "length", "read_depth", "S", "Watterson", "Pi", "Tajima_D",
    "window_start", "window_end", "window_mid"
]

#### Summary collection for deleting missing values
summary_stats = []


#### Process Each .stats File 
for input_file in glob.glob(os.path.join(input_dir, "*.stats")):
    basename = os.path.basename(input_file)
    sample_id = basename.replace("_Celeg_TT_arrow_pilon.sub.mpileup.stats", "")
    output_file = os.path.join(output_dir, sample_id + ".with_pos.cleaned.tsv")

    print(f"Processing: {basename}")
    df = pd.read_csv(input_file, sep="\t")

    # Add genomic position columns
    df["window_start"] = (df["window"] - 1) * window_size + 1
    df["window_end"] = df["window_start"] + window_size - 1
    df["window_mid"] = df["window_start"] + (window_size // 2) -1

    # Keep relevant columns
    df = df[columns_to_keep]

    # Check missing values
    missing = df.isna().sum()

    # Drop rows with any missing values
    before = df.shape[0]
    df_clean = df.dropna()
    after = df_clean.shape[0]
    dropped = before - after

    # Save cleaned DataFrame: after dropping missing values
    df_clean.to_csv(output_file, sep="\t", index=False)
    print(f"Cleaned file saved: {basename}")
    


    # Collect summary info about missing values
    summary_row = {
        "sample": sample_id,
        "rows_before": before,
        "rows_after": after,
        "rows_dropped": dropped
    }
    for col in columns_to_keep:
        summary_row[f"NA_{col}"] = missing.get(col, 0)
    
    summary_stats.append(summary_row)


####  Summary of NA removal and column completeness per sample
summary_df = pd.DataFrame(summary_stats)

#####  save summary to file
summary_df.to_csv(os.path.join(output_dir, "missing_values_per_sample_summary.tsv"), sep="\t", index=False)
print(f"\n missing_values_per_sample_summary saved")

```

### input
there are 15 samples and each sample has 6 chromosome-specific stats file (except mitochondrion): 15 * 6 in total
```bash
analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-02_run_npstats

# one sample example
Sample_XJ-4126-9_hcontortus_chr1_Celeg_TT_arrow_pilon.sub.mpileup.stats
Sample_XJ-4126-9_hcontortus_chr2_Celeg_TT_arrow_pilon.sub.mpileup.stats
Sample_XJ-4126-9_hcontortus_chr3_Celeg_TT_arrow_pilon.sub.mpileup.stats
Sample_XJ-4126-9_hcontortus_chr4_Celeg_TT_arrow_pilon.sub.mpileup.stats
Sample_XJ-4126-9_hcontortus_chr5_Celeg_TT_arrow_pilon.sub.mpileup.stats
Sample_XJ-4126-9_hcontortus_chrX_Celeg_TT_arrow_pilon.sub.mpileup.stats
```

### output
one missing value summary; 15 samples in total, each sample has 6 chromosome-specific `.with_pos.cleaned.tsv`:

```bash

missing_values_per_sample_summary.tsv

# one sample example
Sample_XJ-4126-9_hcontortus_chr1.with_pos.cleaned.tsv
Sample_XJ-4126-9_hcontortus_chr2.with_pos.cleaned.tsv
Sample_XJ-4126-9_hcontortus_chr3.with_pos.cleaned.tsv
Sample_XJ-4126-9_hcontortus_chr4.with_pos.cleaned.tsv
Sample_XJ-4126-9_hcontortus_chr5.with_pos.cleaned.tsv
Sample_XJ-4126-9_hcontortus_chrX.with_pos.cleaned.tsv

```



## 2. Merge chromosome-specific tsv files for each sample

1. for each sample, merge all cleaned chromosome-specific tsv files into a single file
2. sort chromosomes in order

### Python script
```python

import os
import glob
import pandas as pd
import re


input_dir = "../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-02_run_npstats/01_preprocessing_npstat_data" 
output_dir = os.path.join(input_dir, "merged_sample_chr")
os.makedirs(output_dir, exist_ok=True)

# 1: group files by sample name
all_files = glob.glob(os.path.join(input_dir, "*.with_pos.cleaned.tsv"))

# a dictionary: { sample_id: [list of files] }
sample_files = {}
for file_path in all_files:
    fname = os.path.basename(file_path)
    match = re.match(r"(Sample_XJ-\d+-\d+).*?_hcontortus_(chr[0-9XYM]+)", fname)
    if match:
        sample_id = match.group(1)
        sample_files.setdefault(sample_id, []).append(file_path)

# 2: for each sample, merge all chromosome files and add chr name col
for sample_id, file_list in sample_files.items():
    dfs = []
    for filepath in file_list:
        df = pd.read_csv(filepath, sep="\t")
        chr_match = re.search(r"chr[0-9XYM]+", os.path.basename(filepath))
        if chr_match:
            df["chr"] = chr_match.group(0)
        dfs.append(df)
    
    merged_df = pd.concat(dfs, ignore_index=True)

    # Sort chromosomes in custom order, then by window
    chr_order = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chrX']
    merged_df['chr'] = pd.Categorical(merged_df['chr'], categories=chr_order, ordered=True)
    merged_df = merged_df.sort_values(by=['chr', 'window'])
    
    # Write to output
    output_path = os.path.join(output_dir, f"{sample_id}_pos_chr_merged.tsv")
    merged_df.to_csv(output_path, sep="\t", index=False)
    print(f"Merged file written: {sample_id}")


```

### input

15 samples in total, each sample has 6 chromosome-specific `.with_pos.cleaned.tsv`:

```bash

# one sample example
Sample_XJ-4126-9_hcontortus_chr1.with_pos.cleaned.tsv
Sample_XJ-4126-9_hcontortus_chr2.with_pos.cleaned.tsv
Sample_XJ-4126-9_hcontortus_chr3.with_pos.cleaned.tsv
Sample_XJ-4126-9_hcontortus_chr4.with_pos.cleaned.tsv
Sample_XJ-4126-9_hcontortus_chr5.with_pos.cleaned.tsv
Sample_XJ-4126-9_hcontortus_chrX.with_pos.cleaned.tsv

```


### output
15 `_merged.tsv` for 15 samples
```bash
analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-02_run_npstats/01_preprocessing_npstat_data/merged_sample_chr 

2.2M	Sample_XJ-4126-10_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-11_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-12_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-13_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-14_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-15_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-1_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-2_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-3_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-4_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-5_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-6_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-7_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-8_pos_chr_merged.tsv
2.2M	Sample_XJ-4126-9_pos_chr_merged.tsv

```