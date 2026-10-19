# Nucleotide Diversity (Pi) Analysis
Summary of the R script:

1. **Filter windows of Pi with insufficient coverage**
   - Exclude windows with fewer than 500 bp length covered by reads (5% of the 10 kb window size).
   - Reason: Some windows with extreme high Pi values had only a few bases covered by reads.

2. **Calculate genome-wide and chromosome-wide Pi for each sample/population**
   - Calculate the weighted mean Pi, using `length` (the number of bases covered by reads) as the weight (as Npstat use that to calculate Pi in each window)
   - Generate bar plots of genome-wide and chromosome-wide Pi.
   - For farms with replicates: Use the arithmetic mean of the two replicates to represent the farm.

3. **Generate line plots of Pi along genomic positions**
   - Chromosome-wide: Performed for each chromosome.
   - Genome-wide: Not performed, as the plot would be too long and cluttered.


4. **Examine the distribution of window-based Pi**
   - **Box plots:** Compare the median Pi values among samples.
     - Genome-wide: Performed.
     - Chromosome-wide: Performed.
   - **Wilcoxon rank sum tests:** Compare the distributions of window-based Pi values between samples.
     - Genome-wide: Performed, with a heatmap showing pairwise p-values.
     - Chromosome-wide: Not performed due to the large number of comparisons.
     - For farms with replicates: The geometric mean of the two p-values was used to summarize the results.



```r
## set working directory same to script directory
library(rstudioapi)
script_path <- getActiveDocumentContext()$path
setwd(dirname(script_path))
cat("Now working in:", getwd(), "\n")

library(patchwork)  
library(ggplot2)    

library(stringr)

library(dplyr)
library(viridis)

library(tidyr)
library(purrr)

####### 
input_dir <- "../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-02_run_npstats/01_preprocessing_npstat_data/merged_sample_chr/"

outdir='../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/07-Pi_Tajima_removed_s16_s17/07-02_run_npstats/02_pi_summarize_plot'


####### start common process: 
if (!dir.exists(outdir)){dir.create(outdir, recursive = TRUE)} # create outdir if does not exist

```

```r
####################### A.calculate genome wide and chromosome wide Pi and make bar plot #####################################

######### 1. run all files and calculate genome wide and chromosome wide Pi with weighted mean (weight = `length`)
files <- list.files(input_dir, pattern = "_pos_chr_merged\\.tsv$", full.names = TRUE)

# read first file to get chromosome names (assume all files have same chr set)
df <- read.table(files[1], header=TRUE, sep="\t", stringsAsFactors=FALSE)
chrs <- unique(df$chr)


# result table
pi_summary <- data.frame(sample=character(),
                         Pi_genome=numeric(),
                         matrix(nrow=0, ncol=length(chrs)))
colnames(pi_summary)[3:(2+length(chrs))] <- paste0("Pi_", chrs)

# for log of filtering
filter_log <- data.frame()


for (f in files) {
  df1 <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  sample_id <- sub("_pos_chr_merged\\.tsv$", "", basename(f))
  sample_id <- sub("^Sample_XJ-4126-", "", sample_id)
  
  ####### filter windows: based on the checking below in B section
  ##### remove outliers : only keep windows with length > 500 bp (5%), as our window size is 10kbp
  
  #check how many rows df1$length < 500 and percentage
  removed_n <- sum(df1$length < 500)
  total_n <- nrow(df1)
  removed_pct = removed_n / total_n *100
  
  
  # store filtering stats for this sample
  filter_log <- rbind(filter_log, data.frame(
    sample = sample_id,
    total_windows = total_n,
    removed_windows = removed_n,
    removed_percentage = removed_pct,
    kept_windows = total_n - removed_n
  ))
  
  
  # now filter out windows
  df1 <- df1[df1$length >= 500, ]
  
  
  # get weighted mean Pi of each window in whole genome
  pi_genome <- weighted.mean(df1$Pi, w=df1$length)
  
  # get weighted mean Pi of each window in each chromosome
  pi_chr <- numeric(length(chrs))
  names(pi_chr) <- chrs
  for (cc in chrs) {
    d <- df1[df1$chr == cc, ]
    pi_chr[cc] <- weighted.mean(d$Pi, w=d$length)
  }
  
  # add one row
  new_row <- data.frame(sample=sample_id, Pi_genome=pi_genome, t(pi_chr),
                        check.names=FALSE)
  colnames(new_row)[3:(2+length(chrs))] <- paste0("Pi_", chrs)
  
  pi_summary <- rbind(pi_summary, new_row)
}
pi_summary

filter_log # check if you wanna know how much you filtered


############ 2. change sample id:
#### mapping between old IDs and new names
map <- c(
  "1"  = "A1",
  "6"  = "A2",
  "8"  = "B",
  "3"  = "C",
  "4"  = "D",
  "5"  = "E",
  "2"  = "F1",
  "10" = "F2",
  "7"  = "G1",
  "9"  = "G2",
  "11" = "H",
  "13" = "I",
  #"17" = "J",
  "15" = "K",
  "14" = "L",
  #"16" = "M",
  "12" = "N"
)

# change sample values
pi_summary$sample <- map[as.character(pi_summary$sample)]
pi_summary

# order rows following the map order
pi_summary$sample <- factor(pi_summary$sample, levels = unname(map)) ## Make 'sample' a factor so we can control the row order
pi_summary <- pi_summary[order(pi_summary$sample), ]
pi_summary$sample <- as.character(pi_summary$sample) # Change factor back to normal character column

pi_summary


################ 3. averaged replicates if you want
# groups to merge
groups <- list(
  A = c("A1", "A2"),
  F = c("F1", "F2"),
  G = c("G1", "G2")
)

num_cols <- setdiff(colnames(pi_summary), "sample")

new_tab <- data.frame()

# 1) add merged rows (A, F, G)
for (gname in names(groups)) {
  members <- groups[[gname]]
  sub <- pi_summary[pi_summary$sample %in% members, ]
  
  avg_vals <- colMeans(sub[, num_cols], na.rm = TRUE)
  
  new_row <- data.frame(sample = gname, t(avg_vals), check.names = FALSE)
  colnames(new_row) <- c("sample", num_cols)
  
  new_tab <- rbind(new_tab, new_row)
}

# 2) add the samples that are NOT in any merge group
merged_members <- unlist(groups)
keep_rows <- pi_summary[!pi_summary$sample %in% merged_members, ]

new_tab <- rbind(new_tab, keep_rows)
new_tab

# 3) order rows nicely
order_levels <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N")
new_tab$sample <- factor(new_tab$sample, levels = order_levels)
new_tab <- new_tab[order(new_tab$sample), ]
new_tab$sample <- as.character(new_tab$sample)

new_tab


## round to 6 decimal places
pi_summary[] <- lapply(pi_summary, function(x) 
  if(is.numeric(x)) round(x, 5) else x)

new_tab[] <- lapply(new_tab, function(x) 
  if(is.numeric(x)) round(x, 5) else x)

# save
# write.table(pi_summary,
#             file = file.path(outdir, "15samples-Pi_genome_chr_weightedmean_summary.tsv"),
#             sep = "\t", quote = FALSE, row.names = FALSE)
# 
# write.table(new_tab,
#             file = file.path(outdir, "averaged-replicates-Pi_genome_chr_weightedmean_summary.tsv"),
#             sep = "\t", quote = FALSE, row.names = FALSE)
# 


############## 4. make a bar plot of all farms' genome wide Pi ##############
####using new tab, and averaged replicates already

# try one single plot-genome wide:
# reorder data based on Pi-genome  from small to big
new_tab$sample <- reorder(new_tab$sample, new_tab$Pi_genome)

bar_plot_gen <- ggplot(new_tab, aes(x = sample, y = Pi_genome)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(
    aes(label = sprintf("%.4f", Pi_genome)),
    vjust = -0.5,
    size = 4.3
  ) +
  coord_cartesian(
    ylim = c(0.02, max(new_tab$Pi_genome) + 0.001) # y scale from 0.02
  ) +
  labs(
    x = "Farms",
    y = "Nucleotide diversity (Pi)",
    title = "genome-wide Pi (calculated by weighted mean of windows)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text  = element_text(size = 13),
  )


bar_plot_gen

# ggsave(file.path(outdir, "genomewide_Pi_barplot.png"),
#   plot = bar_plot_gen, width = 12, height = 7, dpi = 300)



################# finally make all bar plots: genome wide + chromosome wide
cols_to_plot <- c("Pi_genome","Pi_chr1", "Pi_chr2",
                  "Pi_chr3", "Pi_chr4", "Pi_chr5", "Pi_chrX")

bar_plot_chr_list <- list()



for (col in cols_to_plot) {
  
  # reorder x based on this column
  new_tab$sample <- reorder(new_tab$sample, new_tab[[col]])
  
  # y axis maximum
  y_max <- max(as.matrix(new_tab[, cols_to_plot]), na.rm = TRUE)
  
  # remove title for genome plot, keep chromosome titles, and change chr -> Chr
  this_title <- sub("^Pi_", "", col)
  if (col == "Pi_genome") {
    this_title <- NULL
  } else {
    this_title <- gsub("^chr", "Chr", this_title)
  }
  
  ## set how to show values: chrx -- on the bar, other chrs --- in the bar
  if (col == "Pi_chrX") {
    label_y     <- new_tab[[col]] + 0.001
    label_color <- "black"
    label_hjust <- -0
  } else {
    label_y     <- new_tab[[col]] - 0.003
    label_color <- "white"
    label_hjust <- -0
  }
  
  ## make the plot
  p <- ggplot(new_tab, aes(x = sample, y = .data[[col]])) +
    geom_col(fill = "steelblue", width = 0.7) +
    
    geom_text(aes(y = label_y,
                  label = sprintf("%.4f", .data[[col]])),
              size = 2,
              family = "Arial",
              #fontface = "bold",
              angle = 90,
              colour = label_color,
              vjust = 0.5,
              hjust = label_hjust
              ) +
    
    
    coord_cartesian(
      ylim = c(0, y_max + 0.003),
      clip = 'off' ## will not cut off the things higher than y axis
      
      ) +
    
    scale_y_continuous(expand = c(0, 0)) + # remove the space above and below the y axis
  
    labs(
      x = NULL,
      y = "Pi",
      title = this_title
    ) +
    
    theme_bw(base_family = "Arial", base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.grid.major = element_line(colour = "gray90", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "gray95", linewidth = 0.2),
      
      plot.margin = margin(6, 8, 6, 8)
    )
  
  bar_plot_chr_list[[col]] <- p
}
#bar_plot_chr_list

########## merge plot by patchwork

p_chr_block <- (bar_plot_chr_list[[2]] + bar_plot_chr_list[[3]]) /
  (bar_plot_chr_list[[4]] + bar_plot_chr_list[[5]]) /
  (bar_plot_chr_list[[6]] + bar_plot_chr_list[[7]])


p_gen_bar <- bar_plot_chr_list[[1]]


## final merged plot
p_all_bar <- p_gen_bar /
  plot_spacer() /
  p_chr_block +
  plot_layout(heights = c(1, 0.15, 3)) +  # modify the space p_gen_bar, space, p_chr_block
  plot_annotation(tag_levels = list("A")) &
  theme(
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 14,
      colour = "black"
    ),
    plot.tag.position = c(0.01, 0.99)
  )

p_all_bar


# save it to final figure dir
final_figures_outdir= '../../../../../final_figures'

# save
ggsave(
  file.path(final_figures_outdir, "Pi_barplot_genome_chr_merged.pdf"),
  plot = p_all_bar,
  width = 6.69,   
  height = 6, 
  units = "in",
  device = cairo_pdf
)



```

```r
###################### B.line plots,show genome location vs Pi ##############################################################

############ make Line plots of Pi for each chromosome
##### since it would be too small to show Pi across the whole genome

# read all files
all_dat <- data.frame()

######## 1. combines all files into one big data frame by stacking rows from every file
for (f in files) {
  df_p <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  sample_id <- sub("_pos_chr_merged\\.tsv$", "", basename(f))
  sample_id <- sub("^Sample_XJ-4126-", "", sample_id)
  
  df_p$sample <- map[as.character(sample_id)]  # write sample name
  
  df_p$mid <- df_p$window_mid / 1000  # window mid position within chr (kbp)
  
  all_dat <- rbind(all_dat, df_p)
}


######## 2. set sample order (controls legend order)
sample_order <- unname(map)
all_dat$sample <- factor(all_dat$sample, levels = sample_order)

# chromosomes list
chrs <- unique(all_dat$chr)



######## 3. filter Pi to get cleaned_all_data, which will be used for all plots and analysis below !! #############
######### check outliers in the plots: check the highest Pi window (peaks)
### Focus on 'length'= Number of base pairs covered by sequences: The number of base pairs in a window that passed coverage and quality filters and were used to calculate π.
### If 'length' is very small, π is estimated from only a few sites, making the estimate unstable and potentially producing artificially high peaks, especially when the number of segregating sites is small.


for (cc in chrs) {
  sub <- all_dat[all_dat$chr == cc, ]
  
  # top 10 highest Pi windows
  top <- sub[order(-sub$Pi), ][1:10, c("sample","chr","window_start","window_end","length","read_depth","mid","Pi","S")]
  cat("\n===== ", cc, " =====\n")  # print chromosome title
  print(top)                     # print the top table
}



##### remove outliers in line plots: only keep windows with length > 500 bp (5%), as our window size is 10kbp
cleaned_all_data <- all_dat[all_dat$length >= 500, ]



######## 4. make one line plot per chromosome ###############

##### 4.1 set same color for the same sample in plots
unique_samples <- sort(unique(cleaned_all_data$sample))
unique_samples

# generate a named vector of rainbow colors (in order) using viridis (same palette style as before)
sample_cols <- setNames(
  viridisLite::viridis(length(unique_samples), option = "turbo", begin = 0.05, end = 0.95),
  unique_samples
)
sample_cols



##### 4.2 scale y axis: global y range across ALL chromosomes (and all samples)
y_lim <- range(cleaned_all_data$Pi, na.rm = TRUE);y_lim   # so the maximum should be like 0.1 


##### 4.3 make line plots
plot_line_list <- list()

for (cc in chrs) {
  sub <- cleaned_all_data[cleaned_all_data$chr == cc, ]
  
  # title: change chr -> Chr
  this_title <- gsub("^chr", "Chr", cc)

  p_line <- ggplot(sub, aes(x = mid, y = Pi, group = sample, color = sample)) +
    geom_line(linewidth = 0.3, alpha = 0.65) +
    guides(color = guide_legend(override.aes = list(linewidth = 1))) + #control the legend lines' width
    
    # 1) clear distinct colors for many groups
    #scale_color_viridis_d(option = "turbo", begin = 0.05, end = 0.95, drop = FALSE) +  # use colors for samples randomly
    scale_color_manual(values = sample_cols, drop = FALSE) + # fixed each sample's color manually
    
    # 2) more detailed axes
    scale_x_continuous(breaks = scales::breaks_pretty(n = 12)) +
    scale_y_continuous(limits = y_lim, breaks = scales::breaks_pretty(n = 8)) +
    
    labs(title = this_title, x = "Position (kb)", y = "Pi", color = "Sample") +
    
    theme_bw(base_family = "Arial", base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5), # hjust = 0.5, title is in the middle
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      axis.text.x  = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
      axis.text.y  = element_text(size = 7.7),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      
      panel.grid.major = element_line(colour = "gray90", linewidth = 0.3), 
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(6, 8, 6, 8)
      
    )
  
  # store each plot in list (key = chromosome name)
  plot_line_list[[cc]] <- p_line
  
  ### save each plot if you want
  # ggsave(
  #   filename = file.path(outdir, paste0("Pi_lineplot_", cc, "_kbp.png")),
  #   plot = p, width = 12, height = 4.5, dpi = 300)
}


## merge plot by patchwork
#combine 6 plots with patchwork
p_line_chr_all <- (plot_line_list[[1]] + plot_line_list[[2]]) /
  (plot_line_list[[3]] + plot_line_list[[4]]) /
  (plot_line_list[[5]] + plot_line_list[[6]]) +
  plot_layout(guides = "collect") & theme(legend.position = "right") # merged legends to one
#p_line_chr_all


# save it to final figure dir
final_figures_outdir= '../../../../../final_figures'

# save
ggsave(
  file.path(final_figures_outdir, "Pi_lineplot_ALL_chr_kbp.pdf"),
  plot = p_line_chr_all,
  width = 12, 
  height = 6.69, 
  units = "in",
  device = cairo_pdf
)


```


```r

###################### c.box plots for distribution of Pi ##############################################################
## use  cleaned_all_data

### 1. compute median per pop within entire genome
med_df <- cleaned_all_data %>%
  group_by(sample) %>%
  summarise(medPi = median(Pi, na.rm = TRUE), .groups = "drop")
med_df

### 2. if you wanna order samples by median on the plot; delete that section if you do not wanna order...
# 1) order levels by median (small -> large)
sample_order <- med_df %>%
  arrange(medPi) %>%
  pull(sample)       ## extract the sample column as a vector (this will be the factor level order)

# 2) convert data sample into an ordered factor using sample_order
# here: make data copies (do NOT change cleaned_all_data for next analysis)
cleaned_all_data_ordered <- cleaned_all_data %>%
  mutate(sample = factor(sample, levels = sample_order)) 

med_df <- med_df %>%
  mutate(sample = factor(sample, levels = sample_order)) # # ensure med_df uses the same sample order,



######## 3. geonome:  boxplot + all window points in the genome + median labels

#if you wanna order samples by median on the plot, use data: cleaned_all_data_ordered
#if you do not wanna order samples , use data: cleaned_all_data


### 3.1 generate a named vector of rainbow colors (in order) using viridis ; used ordered samples by median within entire genome
sample_cols_ordered <- setNames(
  viridisLite::viridis(length(sample_order), option = "turbo", begin = 0.05, end = 0.95),
  sample_order
)
sample_cols_ordered


### 3.2 make plot
p_box_genome <- ggplot(cleaned_all_data_ordered, aes(x = sample, y = Pi, color = sample)) +
  
  # scatter points for each window (distribution)
  geom_point(
    position = position_jitter(width = 0.3, height = 0), # spread points in x
    alpha = 0.4,
    size  = 1
  ) +
  
  # boxplot on top (shows median as the thick line)
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    linewidth = 0.25,
    color = "black"    
  ) +
  
  # median value label
  geom_text(
    data = med_df,
    aes(x = sample, y = medPi + 0.005, label = sprintf("%.4f", medPi)),    # y =, the y coordinate of the value: median + ...
    inherit.aes = FALSE, 
    #angle = 90,   # flip 90 degrees
    #nudge_x = -0.42,   # move label to the left of the box
    hjust = 0.5,
    vjust = 0.5,
    color = "black",
    size = 3,
    family = "Arial",
    fontface = "bold"
  ) +
  
  
  # colors
  #scale_color_viridis_d(option = "turbo", begin = 0.05, end = 0.95, drop = FALSE) +  # use colors for samples randomly
  scale_color_manual(values = sample_cols_ordered, drop = FALSE) + # fixed each sample's color manually
  coord_cartesian(ylim = y_lim) +   # set y scale same for all plots
  
  labs(
    x = NULL,
    y = "Pi"
    #color = "Samples"  # legend title
  ) +
  
  theme_bw(base_family = "Arial") +
  theme(
    legend.position = "none",                      
    #legend.title = element_text(size = 13),
    #legend.text  = element_text(size = 13),
    
    plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 9),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.3), 
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(8, 8, 8, 8)
    
    )
#p_box_genome


######## 4. per chromosome boxplot + all window points per chromosome (6 chromosomes) + median labels
plot_list_box <- list()   

for (cc in chrs) {
  
  sub <- cleaned_all_data[cleaned_all_data$chr == cc, ]
  
  # title: change chr -> Chr
  this_title <- gsub("^chr", "Chr", cc)
  
  

  
  # 4.1 compute median per pop within this chromosome
  med_df <- sub %>%
    group_by(sample) %>%
    summarise(medPi = median(Pi, na.rm = TRUE), .groups = "drop")
  
  ### 4.2. if you wanna order samples by median on the plot; delete that section if you do not wanna order
  # 1) order levels by median (small -> large)
  sample_order <- med_df %>%
    arrange(medPi) %>%
    pull(sample)       ## extract the sample column as a vector (this will be the factor level order)
  
  # 2) convert data sample into an ordered factor using sample_order
  sub <- sub %>%
    mutate(sample = factor(sample, levels = sample_order)) 
  
  med_df <- med_df %>%
    mutate(sample = factor(sample, levels = sample_order)) # # ensure med_df uses the same sample order,
  
  
  # 4.3 boxplot + jitter points + median labels
  p_box_chr <- ggplot(sub, aes(x = sample, y = Pi, color = sample)) +
    
    # all windows as jittered points
    geom_point(
      position = position_jitter(width = 0.3, height = 0),
      alpha = 0.4, size = 1
    ) +
    
    # boxplot (median is the thick line)
    geom_boxplot(
      width = 0.8,
      outlier.shape = NA,
      alpha = 0.15,
      linewidth = 0.25,
      color = "black"
    ) +
    
    # add median value as text label
    geom_text(
      data = med_df,
      aes(x = sample, y = medPi + 0.0032, label = sprintf("%.4f", medPi)),    # y =, the y coordinate of the value: median + ...
      inherit.aes = FALSE, 
      #angle = 90,   # flip 90 degrees
      #nudge_x = -0.42,   # move label to the left of the box
      hjust = 0.5,
      vjust = 0.5,
      color = "black",
      size = 2,
      family = "Arial",
      fontface = "bold"
    ) +
    
    
    #scale_color_viridis_d(option = "turbo", begin = 0.05, end = 0.95, drop = FALSE) +
    scale_color_manual(values = sample_cols_ordered, drop = FALSE) + # fixed each sample's color manually
    coord_cartesian(ylim = y_lim) +   # set y scale same for all plots
    
    labs(title = this_title, 
         x = NULL, 
         y = "Pi"
         ) +
    
    theme_bw(base_family = "Arial") +
    theme(
      legend.position = "none",                      
      #legend.title = element_text(size = 13),
      #legend.text  = element_text(size = 13),
      
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text  = element_text(size = 9),
      panel.grid.major = element_line(colour = "gray90", linewidth = 0.3), 
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(8, 8, 8, 8)
    
  )
  
  plot_list_box[[cc]] <- p_box_chr
  
  # optional: save each chromosome plot
  # ggsave(
  #   filename = file.path(outdir, paste0("Pi_boxplot_", cc, ".png")),
  #   plot = p_box_chr, width = 10, height = 5, dpi = 300
  # )
}

# Example: show one chromosome plot
#plot_list_box[[1]]



### 5. merge genome and chromosome box plots by patchwork

p_box_chr_block <- (plot_list_box[[1]] + plot_list_box[[2]]) /
    (plot_list_box[[3]] + plot_list_box[[4]]) /
    (plot_list_box[[5]] + plot_list_box[[6]])

p_all_boxes <- p_box_genome / 
  plot_spacer() /
  p_box_chr_block +
  plot_layout(heights = c(1, 0.15, 3)) +
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
  file.path(final_figures_outdir, "Pi_boxplot_genome_chr_merged.pdf"),
  plot = p_all_boxes,
  width = 12,  
  height = 8.8,  
  units = "in",
  device = cairo_pdf
)





```


```r

###################### D.statistics test ##############################################################
## use cleaned_all_data

########## 1. check normality: confirm the absence of normality for each sample's Pi
### Shapiro-Wilk normality test: 
## results:  all populations do not follow normal distribution  (all P values < 0.05)

set.seed(1)  # reproducible random sampling

shapiro_all_pops <- cleaned_all_data %>%
  group_by(sample) %>%
  summarise(
    n_total = n(),                       # total windows in this pop
    n_test  = min(n_total, 5000),        # sample size used in Shapiro：extract 5000 Pi
    p_value = {
      v <- Pi
      v_sub <- sample(v, size = min(length(v), 5000), replace = FALSE)
      shapiro.test(v_sub)$p.value
    },
    .groups = "drop"
  ) %>%
  arrange(p_value)

shapiro_all_pops




########## 2. check overall difference across 15 pops (genome-wide, chromosome-wide):This is an important step before performing pairwise comparisons to see any pop is sig different 
#### Kruskal–Wallis test: not require normal distribution
## The test assess overall differences in Pi among populations. 
## determines whether at least one population differs from the others in distribution or median Pi, but does not specify which populations are different.


## genome wide test：
# result:    p-value < 2.2e-16, yes, at least one population differs significantly
kw_genome <- kruskal.test(Pi ~ sample, data = cleaned_all_data)
kw_genome


## per chromosome test
# result:    all chromosome p-value < 0.05, yes, at least one population differs significantly per chromosome
kw_by_chr <- cleaned_all_data %>%
  group_by(chr) %>%
  summarise(
    p_value = kruskal.test(Pi ~ sample)$p.value,
    .groups = "drop"
  ) %>%
  arrange(p_value)

kw_by_chr



########## 3. genome wide pairwise.wilcox.test for all samples

###The Wilcoxon rank sum test (Mann–Whitney U) is a non-parametric test used to compare two independent groups. 
# It tests whether one group tends to have larger or smaller values than the other (i.e., a shift in the distribution/median), 
# without assuming the data are normally distributed.

### 1st: for replicates: we use geometric mean to combine P value :  average on the log scale
#### 2nd: then use Benjamini–Hochberg correction


###### a. do pairwise test for all samples includes replicates ###### 
## genome wide: do pairwise.wilcox.test for all samples pairs
pw <- pairwise.wilcox.test(
  x = cleaned_all_data$Pi,
  g = cleaned_all_data$sample, # compare pairwise of g
  p.adjust.method = "none",   # No multiple testing correction for now (later will use Benjamini–Hochberg
  exact = FALSE             # to avoid warning
)

## extract P value of each unique pair to dataframe
res_allsample_genome <- as.data.frame(as.table(pw$p.value)) %>%
  filter(!is.na(Freq)) %>%
  transmute(
    group1 = Var1,
    group2 = Var2,
    p_value = Freq
  )
res_allsample_genome

# check p value
range(res_allsample_genome$p_value)
any(res_allsample_genome$p_value == 0)
any(is.na(res_allsample_genome$p_value))




###### b.  combine p value of 2 replicates of 3 farms each###### 
## a function: get p for a pair (order-insensitive, no matter A vs B or B vs A)
get_p <- function(res, a, b) {
  out <- res %>%
    dplyr::filter((group1 == a & group2 == b) | (group1 == b & group2 == a)) %>% ## extract a vs b OR b vs a 's P-value
    dplyr::pull(p_value)
  if (length(out) != 1) stop("Missing or duplicated pair: ", a, " vs ", b)
  out
}


## function: geometric mean to combine replicate p-values (average on the log scale)
geo_mean_p <- function(pvec) {
  exp(mean(log(pvec)))
}


## set replicates list ( A/F/G farms have two replicates each)
rep_sets <- list(    # A, F G are set name,
  A = c("A1", "A2"),
  F = c("F1", "F2"),
  G = c("G1", "G2")
)



## a function； Combine a set of 2 replicate p-values (e.g. A1/A2) against all other samples

combine_one_set <- function(set_name) { # set_name= A, F or G
  
  reps <- rep_sets[[set_name]]              # reps = the two replicates (e.g. A1, A2)
  targets <- setdiff(unique_samples, reps)  # all other samples
  
  purrr::map_dfr(targets, function(x) { # map_dfr--a loop：for each 'other sample name' as x, do function(x), and row bind all output dfs to one dataframe
    
    p1 <- get_p(res_allsample_genome, reps[1], x) # get p value of A1 and x
    p2 <- get_p(res_allsample_genome, reps[2], x) # get p value of A2 and x
    
    tibble::tibble( # create a output df
      replicate_set = set_name,
      target = x,
      p1 = p1,
      p2 = p2,
      p_geo = geo_mean_p(c(p1, p2)) )
  })
}

# Apply the function to each replicate set (A, F, G) and row-bind all results
combined_raw <- purrr::map_dfr(names(rep_sets), combine_one_set) # names(rep_sets) is A, F or G
combined_raw


## Now, continue using geometric mean to combine

p_A_G <- geo_mean_p(
  combined_raw$p_geo[combined_raw$replicate_set == "A" & combined_raw$target %in% c("G1", "G2")]
)
p_A_G

p_A_F <- geo_mean_p(
  combined_raw$p_geo[combined_raw$replicate_set == "A" & combined_raw$target %in% c("F1", "F2")]
)
p_A_F

p_F_G <- geo_mean_p(
  combined_raw$p_geo[combined_raw$replicate_set == "F" & combined_raw$target %in% c("G1", "G2")]
)
p_F_G

# add new rows
combined_raw <- rbind(
  combined_raw,
  data.frame(replicate_set="A", target="G", p1=NA, p2=NA, p_geo=p_A_G),
  data.frame(replicate_set="A", target="F", p1=NA, p2=NA, p_geo=p_A_F),
  data.frame(replicate_set="F", target="G", p1=NA, p2=NA, p_geo=p_F_G)
)


# delete old rows
all_reps <- unname(unlist(rep_sets))# all replicates A1/A2/F1/F2/G1/G2
all_reps

combined_raw <- combined_raw[!combined_raw$target %in% all_reps, ]
combined_raw

# drop p1/p2, and rename columns
combined_raw <- combined_raw[, c("replicate_set", "target", "p_geo")]
names(combined_raw) <- c("group1", "group2", "p_value")

combined_raw


###### c. add pairwise p-values for non-replicate farms ######

# identify farms that are NOT replicates
all_reps
unique_samples
non_rep_farms <- setdiff(unique_samples, all_reps)
non_rep_farms

# extract pairwise p-values among non-replicate farms
res_nonrep <- res_allsample_genome %>%
  dplyr::filter(group1 %in% non_rep_farms & group2 %in% non_rep_farms)
res_nonrep

# append to combined_raw
combined_raw_all <- rbind(combined_raw, res_nonrep)

combined_raw_all ## 66 rows, 12 unique farms in total


### d. multiple testing correction: Benjamini–Hochberg
combined_final_padj <- combined_raw_all %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(p_adj_BH) # sort

# scientific notation with 3 decimal
combined_final_padj$p_adj_BH <- formatC(combined_final_padj$p_adj_BH, format = "e", digits = 3)
combined_final_padj$p_value  <- formatC(combined_final_padj$p_value,  format = "e", digits = 3)

combined_final_padj



####### d. make a heatmap plot for genome wide pairwise.wilcox.test for all samples
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(patchwork)

df_heatmap <- combined_final_padj %>%
  mutate(
    p_adj_BH_num = as.numeric(p_adj_BH),      # convert from chr -> numeric
    p_show = formatC(p_adj_BH_num, format="e", digits=3),  # numbers to display
    log_val = -log10(p_adj_BH_num)           # use -log10 values for heatmap colors： e.g -log10(7.009e-84) =~ 84
  )

# all sample names
samps <- sort(unique(c(df_heatmap$group1, df_heatmap$group2)))


### build an empty matrix and fill it symmetrically (for -log10 values for color)
mat <- matrix(NA, nrow=length(samps), ncol=length(samps), dimnames=list(samps, samps))

### another matrix for p values (strings)
mat_num <- matrix("", nrow=length(samps), ncol=length(samps), dimnames=list(samps, samps))

for (i in 1:nrow(df_heatmap)) {
  s1 <- df_heatmap$group1[i]
  s2 <- df_heatmap$group2[i]
  
  mat[s1, s2] <- df_heatmap$log_val[i]
  mat[s2, s1] <- df_heatmap$log_val[i]
  
  mat_num[s1, s2] <- df_heatmap$p_show[i]
  mat_num[s2, s1] <- df_heatmap$p_show[i]
}



#### order rows and cols by decreasing mean (same logic as our FST heatmap example)
ord <- order(colMeans(mat, na.rm=TRUE), decreasing=TRUE)
mat <- mat[rev(ord), rev(ord)]
mat_num <- mat_num[rev(ord), rev(ord)]


# keep only lower triangle
mat[upper.tri(mat)] <- NA
mat_num[upper.tri(mat_num)] <- ""

### purple palette
purples <- colorRampPalette(brewer.pal(9, "Purples"))(50)

### create a color matrix: if mat >..., text is white(True), otherwise, black(False); and NA will be NA, so nothing here
text.cols <- ifelse(mat > quantile(mat, 0.75, na.rm=TRUE), "white", "black")

title_gen <- "Heatmap of genome-wide Pi pairwise Wilcoxon test (BH-adjusted p)"


## make heatmap: 
## the legend means -log10(P value), no title of legend, so write it in plot's caption
## cells show P values
ph_gen <- pheatmap(
  mat,                         # use -log10 value to make the heatmap, and show colors
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  display_numbers = mat_num,     # show formatted p-values on the heatmap
  number_color    = text.cols,   # a matrix/vector of colors matching mat
  na_col          = "white",
  color           = purples,
  angle_col       = 0,         # make x-axis labels vertical
  border_color    = "grey90", # lighten the grid lines
  #main            = title_gen, # title
  fontsize = 10, # general text size, can control like legend
  fontsize_number = 9.5, # cell text size
  fontsize_row    = 12, # y axis size
  fontsize_col    = 12,  # x axis size
  fontfamily ="Arial"
)


### convert pheatmap to ggplot for patchwork (when you need)
ph_gen_gg <- as.ggplot(ph_gen$gtable)


# save it to final figure dir
final_figures_outdir= '../../../../../final_figures'

## save as a full page plot
ggsave(
  filename = "Pi_genomewide_wilcoxtest_allsample.pdf",
  plot = ph_gen_gg,
  path = final_figures_outdir,
  width = 10,
  height = 6, 
  units = "in", # inch
  device = cairo_pdf  # as pdf vector figure
)




```