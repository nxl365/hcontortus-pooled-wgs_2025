# Input files and preprocessing

1. original outputs from 03_FST_sliding_window
```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr

# for 3 window sizes
all_samples_rm_s16_s17_1000_window_500_step.by_chr.fst  
all_samples_rm_s16_s17_10000_window_5000_step.by_chr.fst # uploaded to this github `Analysis_results/6. FST`
all_samples_rm_s16_s17_100000_window_50000_step.by_chr.fst

```
2. preprocessed
```bash
preprocessed_all_samples_rm_s16_s17_1000_window_500_step.by_chr.fst
preprocessed_all_samples_rm_s16_s17_10000_window_5000_step.by_chr.fst # uploaded to this github `Analysis_results/6. FST`
preprocessed_all_samples_rm_s16_s17_100000_window_50000_step.by_chr.fst
```

3. For farms A, F, and G, mean FST values of the two biological replicates were calculated
```bash
mean_rep_preprocessed_all_samples_rm_s16_s17_1000_window_500_step.by_chr.fst
mean_rep_preprocessed_all_samples_rm_s16_s17_10000_window_5000_step.by_chr.fst # uploaded to this github `Analysis_results/6. FST`
mean_rep_preprocessed_all_samples_rm_s16_s17_100000_window_50000_step.by_chr.fst
```


# 1. generate heatmaps for genome-wide and chromosome-wide pairwise sliding-window FST within and between biological replicates from three farms
- Lower-triangular heatmaps.
- for three different FST window sizes.
- **Figures generated in the manuscript:** Figure 3 and Figures S1–S2.


## R script

```r
## set working directory same to script directory
library(rstudioapi)
script_path <- getActiveDocumentContext()$path
setwd(dirname(script_path))
cat("Now working in:", getwd(), "\n")

library(stringr)


####### config: run with different inputs and outdir for 3 windowsizes separately

# input <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_1000_window_500_step/preprocessed_all_samples_rm_s16_s17_1000_window_500_step.by_chr.fst'

# outdir <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_1000_window_500_step/02_fst_heatmap_replicates'


input <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_10000_window_5000_step/preprocessed_all_samples_rm_s16_s17_10000_window_5000_step.by_chr.fst'   

outdir <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_10000_window_5000_step/02_fst_heatmap_replicates'


# input <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_100000_window_50000_step/preprocessed_all_samples_rm_s16_s17_100000_window_50000_step.by_chr.fst'

# outdir <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_100000_window_50000_step/02_fst_heatmap_replicates'


####### start common process: 
if (!dir.exists(outdir)){dir.create(outdir, recursive = TRUE)} # create outdir if does not exist


fn <- basename(input) # e.g. "preprocessed_all_samples_rm_s16_s17_1000_window_500_step.by_chr.fst"
desc     <- str_match(fn, ".*_([0-9]+_window_[0-9]+_step)\\.by_chr\\.fst$")[,2]
title_chr <- str_glue(" Heatmap of Chromosome-Wise Pairwise FST among Replicate Samples ({desc})")
title_gen <- str_glue(" Heatmap of Genome-Wide Pairwise FST among Replicate Samples ({desc})")
output_name_chr  <- paste0("fst_heatmap_replicates_chromosomewise_", desc, ".png")
output_name_gen  <- paste0("fst_heatmap_replicates_genomewide_", desc, ".png")

output_name_merge  <- paste0("fst_heatmap_replicates-merged_chr_geno-", desc, ".pdf")


## data: read once
data <- read.table(input, header = T)


########## 1. select only replicates samples：##########
keep_samples<- c(1,6,2,10,7,9) # replicates

keep_samples_string <- paste(keep_samples, collapse="|")  #get single string "1|2|6|7|9|10", '|' will serve as “OR”
keep_samples_string

# build the regax pattern:  paste0--paste with no separator, wanna pattern: s..._VS_s...
pat <- paste0(
  "^s(",
  keep_samples_string,
  ")_VS_s(",
  keep_samples_string,
  ")$"
  )

# check how many names 
matched_names <- grep(pat, names(data),value = T)
matched_names

# get selected df
df <- data[,c("CHR","POS","SNPs","Fraction_window","Aver_min_coverage",matched_names)]



chromosomes <- unique(df$CHR)
chromosomes


############ 2. chromosome-wide FST############

############ get weighted mean FST for each chromosome ##########


w_fst_chr_df <- data.frame(matrix(nrow = length(chromosomes), ncol = length(matched_names)))
colnames(w_fst_chr_df)<- matched_names
rownames(w_fst_chr_df)<- chromosomes

for (chr in chromosomes) {
  chr_data <- df[df$CHR == chr, ]
  
  for (col in matched_names) {
    w_mean <- weighted.mean(chr_data[[col]], w=chr_data$SNPs, na.rm=T)
    w_fst_chr_df[chr,col] <- w_mean
  }
  
}
w_fst_chr_df


### mapping between old IDs and new names
map <- c(
  s1  = "A1",
  s2  = "F1",
  s3  = "C",
  s4  = "D",
  s5  = "E",
  s6  = "A2",
  s7  = "G1",
  s8  = "B",
  s9  = "G2",
  s10 = "F2",
  s11 = "H",
  s12 = "N",
  s13 = "I",
  s14 = "L",
  s15 = "K",
  s16 = "M",
  s17 = "J"
)



colnames(w_fst_chr_df) <- sapply(colnames(w_fst_chr_df), function(x) {
  parts<- strsplit(x, '_VS_')[[1]]
  paste0(map[parts[1]],'_VS_',map[parts[2]])
})


  


########## preprocess fst data for each chromosome for ggplot ##########
library(dplyr)
library(tidyverse)

w_fst_chr_df$CHR <- rownames(w_fst_chr_df)  # add a CHR col




#convert to long format
fst_chr_long <- pivot_longer(w_fst_chr_df, cols = -CHR, names_to = "Pair", values_to = "FST") %>%
  separate(Pair, into = c("Sample1", "Sample2"), sep = "_VS_")


# add rows about reverse pairwise fst
fst_chr_long_rev <- fst_chr_long %>%
  relocate(Sample2, .after = CHR) %>%
  relocate(Sample1, .after = Sample2)%>%
  rename(Sample1 = Sample2, Sample2 = Sample1) # switch sample1 and sample2 cols and rename

# change chr to Chr
fst_chr_long_new <- bind_rows(fst_chr_long,fst_chr_long_rev) %>% # combine
  arrange(CHR, Sample1, Sample2)



## set factor
#sample_levels<- c("s1", "s6", "s2", "s10", "s7", "s9") # replicates
sample_levels<- c("A1", "A2", "F1", "F2", "G1", "G2") # replicates

fst_chr_long_new <- fst_chr_long_new %>%
  mutate(Sample1 = factor(Sample1, levels = sample_levels),
         Sample2 = factor(Sample2, levels = sample_levels))%>%  # turn Sample1 and 2 cols into ordered factors
  filter(as.integer(Sample1) > as.integer(Sample2)) # compare the order based on sample_levels you set


fst_chr_long_new$CHR <- factor(
  fst_chr_long_new$CHR,
  levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chrX"),
  labels = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrX")
)

# set group by farms
sample_groups <- c(
  "F1" = "Skåne F",
  "F2" = "Skåne F",
  "G1" = "Västra Götaland G",
  "G2" = "Västra Götaland G",
  "A1" = "Uppland A",
  "A2" = "Uppland A"
)

farm_colors <- c(
  "Skåne F" = "#66C2A5",
  "Västra Götaland G" = "#8DA0CB",
  "Uppland A" = "#FF6699"
)
axis_colors <- farm_colors[sample_groups[sample_levels]] # colors in the same order of samples_level



########### plot with ggplot2 using facet per chromosome ##########
library(ggplot2)
library(ggtext)
library(grid)

## set a fixed theme for genome and chromosome plots both
theme_bmc_heat <- theme_minimal(base_family = "Arial", base_size = 10) +
  theme(
    axis.title.x = element_text(size = 11, colour = "black",margin = margin(t = 12)),
    axis.title.y = element_text(size = 11, colour = "black",margin = margin(r = 12)),
    axis.text.x = element_text(size = 9, face = "bold",colour = axis_colors,margin = margin(t = 4)),
    axis.text.y = element_text(size = 9, face = "bold",colour = axis_colors,margin = margin(r = 4)),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 11, face = "bold", colour = "black", margin = margin(b = 4)),
    
    legend.title = element_text(size = 10, colour = "black"),
    legend.text  = element_text(size = 9, colour = "black"),
    legend.key.height = unit(8, "mm"),
    legend.key.width  = unit(3, "mm"),
    
    panel.spacing.x = unit(5, "mm"),
    panel.spacing.y = unit(8, "mm"),
    
    plot.title = element_blank(),
    plot.caption = element_markdown(size = 10, hjust = 0),
    plot.margin = margin(8, 8, 8, 8)
  )





# Plot per chromosome
p_heat_c <- ggplot(fst_chr_long_new, aes(x = Sample1, y = Sample2, fill = FST)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = round(FST, 3)), size = 2.2, na.rm = TRUE, 
            family = "Arial", colour = "black") + # the text size 
  scale_fill_gradient(low = "white", high = "red", name = expression(F[ST]), na.value = "white") +
  facet_wrap(~CHR) +
  scale_x_discrete(limits = sample_levels, drop = FALSE) +
  scale_y_discrete(limits = sample_levels, drop = FALSE) + # don’t drop any sample levels, even if there’s no data for them
  labs(x = "Samples", y = "Samples") +
  theme_minimal(base_size=12) +
  theme_bmc_heat +
  theme(plot.margin = margin(t = 14, r = 8, b = 8, l = 8))

p_heat_c





############ 3. genome-wide FST ############

########## get weighted mean FST per genome ##########

w_fst_gen_df <- data.frame(matrix(nrow = 1, ncol = length(matched_names)))
colnames(w_fst_gen_df)<- matched_names
rownames(w_fst_gen_df)<- 'genome'

for (col in matched_names) {
  w_mean <- weighted.mean(df[[col]], w=df$SNPs, na.rm=T)
  w_fst_gen_df[col] <- w_mean
}

w_fst_gen_df


colnames(w_fst_gen_df)<- sapply(colnames(w_fst_gen_df),function(x){
  parts <- strsplit(x, '_VS_')[[1]]
  paste0(map[parts[1]],'_VS_',map[parts[2]])
})


########## preprocess fst data for genome for ggplot ##########

w_fst_gen_df$genome <- rownames(w_fst_gen_df)  # add a CHR col

#convert to long format
fst_gen_long <- pivot_longer(w_fst_gen_df, cols = -genome, names_to = "Pair", values_to = "FST") %>%
  separate(Pair, into = c("Sample1", "Sample2"), sep = "_VS_")

# add rows about reverse pairwise fst
fst_gen_long_rev <- fst_gen_long %>%
  relocate(Sample2, .after = genome) %>%
  relocate(Sample1, .after = Sample2)%>%
  rename(Sample1 = Sample2, Sample2 = Sample1) # switch sample1 and sample2 cols and rename

fst_gen_long_new <- bind_rows(fst_gen_long,fst_gen_long_rev) %>% # combine
  arrange(genome, Sample1, Sample2)

## set factor
fst_gen_long_new <- fst_gen_long_new %>%
  mutate(Sample1 = factor(Sample1, levels = sample_levels),
         Sample2 = factor(Sample2, levels = sample_levels))%>%  # turn Sample1 and 2 cols into ordered factors
  filter(as.integer(Sample1) > as.integer(Sample2)) # compare the order based on sample_levels you set




# Plot for genome
p_heat_g <- ggplot(fst_gen_long_new, aes(x = Sample1, y = Sample2, fill = FST)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(FST, 3)), size = 3.5, na.rm = TRUE,
            family = "Arial", colour = "black") +
  scale_fill_gradient(low = "white", high = "red", name = expression(F[ST]), na.value = "white") +
  #facet_wrap(~CHR) +
  scale_x_discrete(limits = sample_levels, drop = FALSE) +
  scale_y_discrete(limits = sample_levels, drop = FALSE) + # don’t drop any sample levels, even if there’s no data for them
  labs(x = "Samples", y = "Samples") +
  theme_bmc_heat +
  theme(plot.margin = margin(t = 8, r = 8, b = 14, l = 8))
  
p_heat_g 



### merge plot by patchwork
library(patchwork)

p_merged <- p_heat_g / p_heat_c +
  plot_layout(heights = c(1, 1.9)) + # height assignment
  plot_annotation(tag_levels = "A") &   # add A and B section  title
  theme(
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 14,
      colour = "black"
    ),
    plot.tag.position = c(0.01, 0.99)
  )


  # save it to final figure dir
final_figures_outdir= '../../../../../final_figures'

## save as a full page plot
ggsave(
  filename = output_name_merge,
  plot = p_merged,
  path = final_figures_outdir,
  width = 6.69,  wide
  height = 8.5,  legend
  units = "in", # inch
  device = cairo_pdf  # as pdf vector figure
)

```




# 2. generate heatmaps for genome-wide and chromosome-wide pairwise sliding-window FST among 12 farms
- Lower-triangular heatmaps.
- arrange the cells in order of descending fst.
- for three different FST window sizes.
- **Figures generated in the manuscript:** Figure 4 and Figures S3–S4.


## R script
```r
## set working directory same to script directory
library(rstudioapi)
script_path <- getActiveDocumentContext()$path
setwd(dirname(script_path))
cat("Now working in:", getwd(), "\n")


library(pheatmap)
library(RColorBrewer)
library(ggplotify)   # for as.ggplot()
library(patchwork)   # for wrap_plots() + plot_annotation()
library(ggplot2)    # for theme(), element_text()

library(stringr)


####### config: run with different inputs and outdir for 3 windowsizes separately

#input = "../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_1000_window_500_step/03_fst_heatmap_12farms/mean_rep_epreprocessed_all_samples_rm_s16_s17_1000_window_500_step.by_chr.fst"

input ="../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_10000_window_5000_step/03_fst_heatmap_12farms/mean_rep_preprocessed_all_samples_rm_s16_s17_10000_window_5000_step.by_chr.fst"

#input="../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_100000_window_50000_step/03_fst_heatmap_12farms/mean_rep_preprocessed_all_samples_rm_s16_s17_100000_window_50000_step.by_chr.fst"


#outdir <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_1000_window_500_step/03_fst_heatmap_12farms'

outdir <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_10000_window_5000_step/03_fst_heatmap_12farms'

#outdir <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_100000_window_50000_step/03_fst_heatmap_12farms'


####### start common process: 
if (!dir.exists(outdir)){dir.create(outdir, recursive = TRUE)} # create outdir if does not exist


fn <- basename(input) # e.g. "mean_rep_preprocessed_all_samples_rm_s16_s17_1000_window_500_step.by_chr.fst"
desc     <- str_match(fn, ".*_([0-9]+_window_[0-9]+_step)\\.by_chr\\.fst$")[,2] # get window and step
title_chr <- str_glue(" Heatmap of Chromosome-Wise Pairwise FST among 12 farms ({desc})")
title_gen <- str_glue(" Heatmap of Genome-Wide Pairwise FST among 12 farms ({desc})")
output_name_chr  <- paste0("fst_heatmap_12farms_chromosomewise_", desc, ".png")
output_name_gen  <- paste0("fst_heatmap_12farms_genomewide_", desc, ".png")

output_name_merge  <- paste0("fst_heatmap_12farms-merged_chr_geno-", desc, ".pdf")


## data: read once
df <- read.table(input, header = T)


## for changing old sample ID to new
map <- c(
  A  = "A",
  F  = "F",
  s3  = "C",
  s4  = "D",
  s5  = "E",
  G  = "G",
  s8  = "B",
  s11 = "H",
  s12 = "N",
  s13 = "I",
  s14 = "L",
  s15 = "K",
  s16 = "M",
  s17 = "J"
)

idx_after6 <- 6:ncol(df)

colnames(df)[idx_after6] <- sapply(colnames(df)[idx_after6], function(x) {
  parts <- strsplit(x, "_VS_")[[1]]
  paste0(map[parts[1]], "_VS_", map[parts[2]])
})


## all unique samples names

fst_cols <- names(df)[6:ncol(df)]


chromosomes <- unique(df$CHR)
chromosomes


samps <- unique(unlist(strsplit(fst_cols,"_VS_")))
samps
length(samps)






########## 1. chromosome-wide FST##########

############ get weighted mean FST for each chromosome ##########


w_fst_chr_df <- data.frame(matrix(nrow = length(chromosomes), ncol = length(fst_cols)))
colnames(w_fst_chr_df)<- fst_cols
rownames(w_fst_chr_df)<- chromosomes

for (chr in chromosomes) {
  chr_data <- df[df$CHR == chr, ]
  
  for (col in fst_cols) {
    w_mean <- weighted.mean(chr_data[[col]], w=chr_data$SNPs, na.rm=T)
    w_fst_chr_df[chr,col] <- w_mean
  }
  
}
w_fst_chr_df


########## make a function for heatmap of each chr


make_chr_heatmap <- function(chr){
  
  ## title: change chr -> Chr
  chr_title<- sub("^chr", "Chr", chr)
  
  ### build an empty matrix and fill it symmetrically
  mat <- matrix(NA, nrow = length(samps), ncol = length(samps), dimnames = list(samps,samps))
  
  # for chr1 as an example: 
  for (col in fst_cols){
    parts <- strsplit(col, "_VS_")[[1]]  # get "s1", "s2"
    sample1<- parts[1]
    sample2<- parts[2]
    mat[sample1,sample2] <- w_fst_chr_df[chr, col]
    mat[sample2,sample1] <- w_fst_chr_df[chr, col]}
  
  
  #### order rows and cols by decreasing mean‐FST
  ord <- order(colMeans(mat, na.rm=TRUE), decreasing=TRUE)
  #mat <- mat[ord, ord]  # arrange mat in the same decreasing order
  mat <- mat[rev(ord), rev(ord)]
  mat[upper.tri(mat)] <- NA  # remove the upper triangle
  #mat


    
  ### drawing heatmap: lower triangle, arrange the cells in order of descending fst

  reds <- colorRampPalette(brewer.pal(9, "Reds"))(50) #build red palette
  
  text.cols <- ifelse(mat > 0.115, "white", "black") # create a color matrix: if mat > 0.14, text is white(True), otherwise, black(False); and NA will be NA, so nothing here
  
  ph<-pheatmap(mat,
               cluster_rows    = F,
               cluster_cols    = F,
               display_numbers = TRUE,
               number_format   = "%.3f",      # round to 3 decimals
               number_color    = text.cols,  # a matrix/vector of colors matching mat
               na_col          = "white",
               #color           = reds,
               color           = colorRampPalette(c("white", "mistyrose", "red", "darkred"))(100),
               breaks          = seq(0, 0.2, length.out = 101), # set fixed legend scale for all chr plots， based on the value range, 0-0.204
               angle_col       = 0,         # make x-axis labels vertical
               border_color    = "grey90",   # lighten the grid lines
               #main            = chr_title,   # title
               fontsize = 7, # general text size
               fontsize_number = 4.3,   # cell text size
               fontsize_row    = 7,   # y axis size
               fontsize_col    = 8,    # x axis size
               fontfamily      = "Arial",
               #legend = (chr == "chr4")   # False or True: only chr4 shows legend; 
  )

  

  
  ### convert the pheatmap gtable into ggplot object for patchwork
  ph_chr_gg<- as.ggplot(ph$gtable) +
    ggtitle(chr_title) +     ## add plot title by ggplot2
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(
        family = "Arial",
        face = "bold",
        size = 10,
        hjust = 0.5    # in the middle
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
  return(ph_chr_gg)
}


## use function to generate all heatmap plots for chrs

plots_chr <- lapply(chromosomes,make_chr_heatmap)


## patchwork: only for ggplots
merge_plots_chr <- wrap_plots(plots_chr, ncol = 2, nrow = 3) +
  plot_annotation(
    title   = title_chr,
    theme   = theme(
      plot.title   = element_text(face="bold"),
      plot.caption = element_text(size=10, hjust=0)
    )
  )

merge_plots_chr






######### 2. genome-wide FST##########

########## get weighted mean FST per genome ##########

w_fst_gen_df <- data.frame(matrix(nrow = 1, ncol = length(fst_cols)))
colnames(w_fst_gen_df)<- fst_cols
rownames(w_fst_gen_df)<- 'genome'

for (col in fst_cols) {
  w_mean <- weighted.mean(df[[col]], w=df$SNPs, na.rm=T)
  w_fst_gen_df[col] <- w_mean
}

w_fst_gen_df



### build an empty matrix and fill it symmetrically
mat <- matrix(NA, nrow = length(samps), ncol = length(samps), dimnames = list(samps,samps))
mat
for (col in fst_cols){
  parts <- strsplit(col, "_VS_")[[1]]  # get "s1", "s2"
  sample1<- parts[1]
  sample2<- parts[2]
  mat[sample1,sample2] <- w_fst_gen_df[1, col]
  mat[sample2,sample1] <- w_fst_gen_df[1, col]}


#### order rows and cols by decreasing mean‐FST
ord <- order(colMeans(mat, na.rm=TRUE), decreasing=TRUE)
#mat <- mat[ord, ord]  # arrange mat in the same decreasing order
mat <- mat[rev(ord), rev(ord)]
mat[upper.tri(mat)] <- NA  # remove the upper triangle
#mat


### drawing heatmap: lower triangle, arrange the cells in order of descending fst

reds <- colorRampPalette(brewer.pal(9, "Reds"))(50) #build red palette

text.cols <- ifelse(mat > 0.12, "white", "black") # create a color matrix: if mat > 0.14, text is white(True), otherwise, black(False); and NA will be NA, so nothing here

ph_gen<-pheatmap(mat,
             cluster_rows    = F,
             cluster_cols    = F,
             display_numbers = TRUE,
             number_format   = "%.3f",      # round to 3 decimals
             number_color    = text.cols,  # a matrix/vector of colors matching mat
             na_col          = "white",
             #color           = reds,
             color           = colorRampPalette(c("white", "mistyrose", "red", "darkred"))(100),
             breaks          = seq(0, 0.16, length.out = 101), # set fixed legend scale， based on the value range, 0-0.16
             angle_col       = 0,         # make x-axis labels vertical
             border_color    = "grey90",   # lighten the grid lines
             #main            =  title_gen,  # title
             fontsize = 7, # general text size
             fontsize_number = 6,   # cell text size
             fontsize_row    = 7,   # y axis size
             fontsize_col    = 8,    # x axis size
             fontfamily      = "Arial"
)

### convert the pheatmap gtable into ggplot object for patchwork
ph_gen_gg<- as.ggplot(ph_gen$gtable) +
  theme(
    plot.margin = margin(8, 8, 8, 8)
  )


plot_gen_final <- ph_gen_gg +
  plot_annotation(
    #caption = "Farms: A = Uppland A; F = Skåne F; G = Västra Götaland G; s14 = Skåne L; s8 = Västra Götaland B; s4 = Västra Götaland D; \n  s11 = Västra Götaland H; s13 = Västra Götaland I; s15 = Västra Götaland K; s3 = Kronoberg C; s5 = Östergötland E; s12 = Småland N",
    theme   = theme(
      plot.caption = element_text(size=10, hjust=0)
    )
  )
## \n : change line

plot_gen_final

  



############ merge plot by patchwork

library(patchwork)

p_merged <- plot_gen_final / plot_spacer() / merge_plots_chr +
  plot_layout(heights = c(1, 0.1, 3)) + # height assignment
  plot_annotation(tag_levels = "A") &   # add A and B section  title
  theme(
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 14,
      colour = "black"
    ),
    plot.tag.position = c(0.01, 0.99)
  )



# save it to final figure dir
final_figures_outdir= '../../../../../final_figures'

## save as a full page plot
ggsave(
  filename = output_name_merge,
  plot = p_merged,
  path = final_figures_outdir,
  width = 6.69,
  height = 8.5,
  units = "in", # inch
  device = cairo_pdf  # as pdf vector figure
)


```