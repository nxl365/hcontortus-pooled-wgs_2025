# plot density curve of each sample ‘s snp and Indels vcf : QUAL and DP
based on the distribution to select filter threshold of QUAL and DP


## R script

### Create the R script:

```bash
020_density_plot_no_mitochondria.R
```

### Script content

```r
## density curve plot for each sample's snp and indels
## remove mitochondria chromosome


library(readr)
library(ggplot2)
library(dplyr)
library(magick)


# define a function to create density plot

make_density_plots <- function(filepath) {
  
  # read the file with all columns as character to avoid parsing issues
  df <- read_tsv(filepath, col_types = cols(.default = col_character()))
  
  ## remove mitochondria chromosome rows: 
  ## if do not wanna remove, just delete this line
  df <- df[!grepl("mitochondrion",df$CHROM, ignore.case = TRUE),]
  
  # change to numeric: QUAL, DP
  df$QUAL <- as.numeric(df$QUAL)
  df$DP   <- as.numeric(df$DP)
  
  # log transform: 
  df$QUAL_log10 <- log10(df$QUAL+1)
  df$DP_log10 <- log10(df$DP+1)
  
  
  # extract sample name: basename will get Sample_....indel_extracted_info.tsv, change _extracted_info.tsv to "" empty, so get Sample_.... = sample_name
  sample_name <- gsub("_extracted_info.tsv$", "",basename(filepath)) 
  
  
  # density curve plot
  
  # for QUAL
  plot_qual <- ggplot(df, aes(x = QUAL_log10)) +
    geom_density(aes(y = ..density..),fill = "steelblue") +
    labs(title = paste("Density Plot of QUAL for", sample_name), x = "QUAL(log10)", y = "Density")
  
  # for DP
  plot_dp <- ggplot(df, aes(x = DP_log10)) +
    geom_density(aes(y = ..density..), fill = "green") +
    labs(title = paste("Density Plot of DP for", sample_name),
         x = "DP (log10)", y = "Density")
  
  
  # Save the plots to the same folder as the input file
  path_name <- dirname(filepath)
  outdir <- file.path(path_name, "density_plot")
  
  ggsave(filename = file.path(outdir, paste0(sample_name, "_QUAL_density_no_mito.png")),
         plot = plot_qual, width = 7, height = 5, dpi = 300)
  ggsave(filename = file.path(outdir, paste0(sample_name, "_DP_density_no_mito.png")),
         plot = plot_dp, width = 7, height = 5, dpi = 300)
  
  cat("Plots saved for", sample_name, "\n")
  
}


### try one sample with the function
# make_density_plots("../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/020_freebayes_extracted_vcf_info/Sample_XJ-4126-1-ampure.indel_extracted_info.tsv")



### run all sample in the folder:
## Process all TSV files in directory that match the pattern.
folder_path <- "../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/020_freebayes_extracted_vcf_info"
tsv_files <- list.files(path = folder_path, pattern = "extracted_info\\.tsv$", full.names = TRUE)

# Apply the function to every file (both indel and snp files)
lapply(tsv_files, make_density_plots)



### combine all images QUAL

plot_path="../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/020_freebayes_extracted_vcf_info/density_plot"

# List all QUAL density PNG files
qual_files <- list.files(path =plot_path, pattern = "QUAL_density_no_mito\\.png$", full.names = TRUE)

# Read in the images
qual_images <- image_read(qual_files)

# Combine the images horizontally (set stack = FALSE) or vertically (stack = TRUE)
combined_qual <- image_append(qual_images, stack = TRUE)


# Save the combined QUAL density image
image_write(combined_qual, path = file.path(plot_path, "combined_QUAL_density_no_mito.png"), format = "png")




### combine all images DP

# List all QUAL density PNG files
qual_files <- list.files(path =plot_path, pattern = "DP_density_no_mito\\.png$", full.names = TRUE)

# Read in the images
qual_images <- image_read(qual_files)

# Combine the images horizontally (set stack = FALSE) or vertically (stack = TRUE)
combined_qual <- image_append(qual_images, stack = TRUE)


# Save the combined QUAL density image
image_write(combined_qual, path = file.path(plot_path, "combined_DP_density_no_mito.png"), format = "png")

```




## Input 

for each sample, one `.indel_extracted_info.tsv`, one `.snp_extracted_info.tsv`:

```bash
analysis/variant_explore/illumina_HC_wgs_XJ4126/020_freebayes_extracted_vcf_info

76M	Sample_XJ-4126-9.indel_extracted_info.tsv
777M	Sample_XJ-4126-9.snp_extracted_info.tsv
...
```
*** Download the files from the HPC cluster to the local machine ***



## Output

* all 17 samples, e.g. for each sample, 4 outputs

```bash
020_freebayes_extracted_vcf_info> ll

Sample_XJ-4126-9.snp_DP_density_no_mito.png
Sample_XJ-4126-9.snp_QUAL_density_no_mito.png
Sample_XJ-4126-9.indel_DP_density_no_mito.png
Sample_XJ-4126-9.indel_QUAL_density_no_mito.png
```

* comebined all samples' plots: 
```bash
# based on these plots to select thresholds
combined_DP_density_no_mito.png
combined_QUAL_density_no_mito.png

```


## Select thresholds
* removed sample 16 and 17 here, due to poor mapping quality and weird distributions of QUAL and DP 
* QUAL: 20   ,x axis of the density plot (log10(20) = 1.3)
* DP:   10   ,x axis of the density plot (log10(10) = 1)