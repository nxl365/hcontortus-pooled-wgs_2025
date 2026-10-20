# 1. generate genome-wide Manhattan plots for selected farm pairs of interest
- Manhattan plots of sliding-window FST using three different window sizes.
- can specify different FST cutoffs.
- **Figures generated in the manuscript:** Figure 5 （main manhattan plots）and Figures S5–S6.

## R script
```r
########## set working directory same to script directory ########## 
library(rstudioapi)
script_path <- getActiveDocumentContext()$path
setwd(dirname(script_path))
cat("Now working in:", getwd(), "\n")




# install.packages("CMplot")
library(CMplot)

# help： ?CMplot
# manual： https://github.com/YinLiLin/CMplot/blob/master/User%20Manual%20for%20CMplot.pdf


library(stringr)


#################  config: run with different inputs and outdir for 3 windowsizes separately

#input = "../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_1000_window_500_step/03_fst_heatmap_12farms/mean_rep_epreprocessed_all_samples_rm_s16_s17_1000_window_500_step.by_chr.fst"

input ="../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_10000_window_5000_step/03_fst_heatmap_12farms/mean_rep_preprocessed_all_samples_rm_s16_s17_10000_window_5000_step.by_chr.fst"

#input="../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_100000_window_50000_step/03_fst_heatmap_12farms/mean_rep_preprocessed_all_samples_rm_s16_s17_100000_window_50000_step.by_chr.fst"


#outdir <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_1000_window_500_step/04_fst_manhattan_plot/final_figure'

outdir <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_10000_window_5000_step/04_fst_manhattan_plot/final_figure'

#outdir <- '../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_100000_window_50000_step/04_fst_manhattan_plot/final_figure'


####### select the cutoff_name for output pdf #########
#cutoff_name <- "-Cut-3sd"
#cutoff_name <- "-Cut-5sd"
#cutoff_name <- "-Cut-top5"
cutoff_name <- "-Cut-top1"  # use this at the end
#cutoff_name <- "-Cut-top0.1"
#cutoff_name <- "-Cut-top0.5"



####### 1.start data preprocess  ########## 
if (!dir.exists(outdir)){dir.create(outdir, recursive = TRUE)} # create outdir if does not exist

fn <- basename(input) # e.g. "mean_rep_preprocessed_all_samples_rm_s16_s17_1000_window_500_step.by_chr.fst"
desc     <- str_match(fn, ".*_([0-9]+_window_[0-9]+_step)\\.by_chr\\.fst$")[,2] # get window and step



data <- read.table(input, header = T)



## CHR column change all chr.. to Chr

data$CHR <- gsub("^chr", "Chr", data$CHR)


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

idx_after6 <- 6:ncol(data)

colnames(data)[idx_after6] <- sapply(colnames(data)[idx_after6], function(x) {
  parts <- strsplit(x, "_VS_")[[1]]
  paste0(map[parts[1]], "_VS_", map[parts[2]])
})


########## 2. make manhattan plots ##########

##########  make a function for manhattan plot:  with CMplot ########## 
make_fst_plot <- function(colname){
  
  
  # build the 4‐col df for per pairwise fst which CMplot needs
  Fst <-data.frame(ID = paste(data[,1], data[,2], sep=":"),
                   chr = data[,1],
                   pos = (data[,2]),       # midpoint of the window  →  used for plotting !!
                   val = data[[colname]]    # from the 6th col: are the fst values
  )
  
  
  ##########  choose one genome wide significance cutoff: ########## 
  
  ## 1: use mean+3σ (sd: standard deviations) fst as threshold
  #cutoff <- (mean(Fst$val) + 3*sd(Fst$val))
  #cutoff_name <- "-Cut-3sd"

  
  ## 2: use mean+5σ fst as threshold
  #cutoff <- (mean(Fst$val) + 5*sd(Fst$val))
  #cutoff_name <- "-Cut-5sd"

  
  ## 3: use top % fst value as threshold 
  #cutoff <- quantile(Fst$val,0.95) # top 5%, generally people use
  #cutoff_name <- "-Cut-top5"
  
  cutoff <- quantile(Fst$val,0.99) # top 1%  ##### use this at the end
  cutoff_name <- "-Cut-top1"

  # cutoff <- quantile(Fst$val,0.999) # top 0.1% = 0.001
  # cutoff_name <- "-Cut-top0.1"

  
  #cutoff <- quantile(Fst$val,0.995) # top 0.5% = 0.005
  #cutoff_name <- "-Cut-top0.5"

  print(cutoff_name)
  print(cutoff)
  

  
  ## plot name:
  plot_name=colname
  
  cols6 <- c(
    "dodgerblue3",
    "darkorange2",
    "forestgreen",
    "firebrick3",
    "mediumpurple3",
    "goldenrod2"
  )
  
  
  
  
  ## change working directory for output file's path
  workwd <- getwd()
  setwd(outdir)
  
  # manhattan plot of fst
  CMplot(Fst,
         type = "p",          # "p" (point), "l" (cross line), "h" (vertical lines)
         plot.type = "m",     # Manhattan plot
         LOG10 = F,           # whether to change the p-value into log10(p-value) scale.
         col = cols6,      # the color for chromosomes，  if the length of the "col" is shorter than the number of chromosomes, then colors will be applied circularly
         cex = 0.3,        # the size of points
         band = 1,       # the size of space between chromosomes
         
         #ylim = c(0, max(Fst$val)), # vector (c(min, max)) or list, CMplot will only plot the points among this interval
         ylim = c(0,1),
         
         threshold = cutoff ,      ## significant threshold: more than one significant line can be added
         threshold.col = 'black',   ## the color for the line of threshold levels
         threshold.lty = 2,       ## the type for the line of threshold levels
         threshold.lwd = 3,       ##  the width for the line of threshold levels, 
         
         amplify = T,             ## amplify the significant points，default: amplify=TRUE
         signal.cex= 1,	  ## a number, if amplify=TRUE, set the size of significant points.
         signal.pch= 1,    # a number, if amplify=TRUE, set the shape of significant points.
         #signal.col=,	  # a character, if amplify=TRUE, set the colour of significant points
         points.alpha= 100, #  transparency of point: 0- 255
         
         ylab=expression(F[ST]),
         ylab.pos = 4,   # the distance between ylab and yaxis.
         axis.cex= 2 , # both axis size
         axis.lwd= 2 , # both axis thickness
         lab.cex= 2 ,  # label size 
         lab.font= 1.3, # label thickness

        
         # main=plot_name,  ## title
         # main.cex= 2.5,  ## size of title
         # main.font=1.2,  ## 	font of title
         
         file.output =T,      #output the plot results, T or F, print to current path
    
         # file = "png", # "jpg", "pdf", "tiff", and "png" can be selected. The "png" has a transparent background.
         
         file = "pdf",
         file.name = plot_name,
         height = 7,   # should be inch
         width = 13
         
         ) 
  
  
  ##change back to working directory
  setwd(workwd)
  getwd() # check
  
  message("Export: ", plot_name)
}


############# choose some interesting pairwise samples to  make manhattan plots  ########## 

## we choose 3 pairs
cols_check <- c("D_VS_K","K_VS_D",
                "D_VS_G",  "G_VS_D",
                "D_VS_C", "C_VS_D")

check_vec <- setNames(cols_check %in% colnames(data), cols_check ) # assigns names/T or F to a vector
check_vec

exists_col <- names(check_vec)[check_vec] # give the names which is True
exists_col # these are the 3 pairs names in the data： "K_VS_D" "G_VS_D" "C_VS_D"


################### save option1: save each plot separately in outdir: ########## 
# if wanna print to outdir: must set CMplot (file.output = T) ; if file.output = F, only print plot in R studio

####### In the final paper: saved the final plots separately, then merged by Affinity

for (col_name in exists_col) {
  make_fst_plot(col_name)
}



####################  save option2: save all plots in one pdf in outdir ########## 
# must set: CMplot (file.output = F)

# pdf_name <- paste0("merged_manhattan_plots_chrs", desc, cutoff_name, ".pdf")
# pdf(file.path(outdir, pdf_name), width=15, height=8)
# for (col_name in exists_col) {
#   p <- make_fst_plot(col_name)
#   print(p)
# }
# dev.off()
# while (length(dev.list()) > 0) dev.off() # make sure the device is null and call it off


```





# 2. generate additional Manhattan plots showing the genomic coordinates of FST peak regions
- focuses on sliding-window FST peak regions on chromosomes 4, 5, and X.
- only uses the 10 kb window size.
- top 1% FST cutoff.
- includes only the 3 most differentiated farm pairs (D–K, D–G, and D–C).
- **Figures generated in the manuscript:** Figure 5 （enlarged panels）


## R script
```r

########## set working directory same to script directory ##########
library(rstudioapi)
script_path <- getActiveDocumentContext()$path
setwd(dirname(script_path))
cat("Now working in:", getwd(), "\n")

library(ggplot2)
library(stringr)
library(grid)

################# config #################
## only use the windowsize 10kb and 5kb step size

input <- "../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_10000_window_5000_step/03_fst_heatmap_12farms/mean_rep_preprocessed_all_samples_rm_s16_s17_10000_window_5000_step.by_chr.fst"

outdir <- "../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-03_poola_fst_by_chr/fst_10000_window_5000_step/04_fst_manhattan_plot/final_figure/peak_manhattan"


if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


####### start data preprocess  ########## 

data <- read.table(input, header = TRUE, check.names = FALSE)

## CHR column: chr -> Chr
data$CHR <- gsub("^chr", "Chr", data$CHR)

## rename old sample IDs to new IDs
map <- c(
  A   = "A",
  F   = "F",
  s3  = "C",
  s4  = "D",
  s5  = "E",
  G   = "G",
  s8  = "B",
  s11 = "H",
  s12 = "N",
  s13 = "I",
  s14 = "L",
  s15 = "K",
  s16 = "M",
  s17 = "J"
)

idx_after6 <- 6:ncol(data)

colnames(data)[idx_after6] <- sapply(colnames(data)[idx_after6], function(x) {
  parts <- strsplit(x, "_VS_")[[1]]
  paste0(map[parts[1]], "_VS_", map[parts[2]])
})


########## choose comparisons and chr 4,5,x ##########

#  3 pairs : 

cols_check <- c(
  "D_VS_K", "K_VS_D",
  "D_VS_G", "G_VS_D",
  "D_VS_C", "C_VS_D"
)

check_vec <- setNames(cols_check %in% colnames(data), cols_check)
check_vec

exists_col <- names(check_vec)[check_vec]
exists_col  # these are the 3 pairs names in the data： "K_VS_D" "G_VS_D" "C_VS_D"


# only keep the 3 pairs in data, 
data <- data[, c(colnames(data)[1:5],
                 "K_VS_D", "G_VS_D", "C_VS_D")]

# only keep chr4,5,X
unique(data$CHR) #"Chr1" "Chr2" "Chr3" "Chr4" "Chr5" "ChrX"

data <- data[data$CHR %in% c("Chr4", "Chr5", "ChrX"), ]



########## setting ##########


## same colours as previous manhattan code
chr_cols <- c(
  "Chr4" = "firebrick3",
  "Chr5" = "mediumpurple3",
  "ChrX" = "goldenrod2"
)



########## plotting function ##########
####### only plot the peaks of chr4,5,x: xmin_mb, xmax_mb is the region coordinator you should select

make_single_chr_plot <- function(colname, chr, xmin_mb, xmax_mb, step_mb = 2) {
  
  ## top 1% cutoff for this population pair
  cutoff <- as.numeric(quantile(data[[colname]], 0.99, na.rm = TRUE))
  
  ## select chromosome + region
  Fst <- data.frame(
    CHR = data$CHR,
    POS = data$POS,
    FST = data[[colname]]
  )
  
  Fst <- Fst[
    Fst$CHR == chr &
      Fst$POS >= xmin_mb * 1e6 &
      Fst$POS <= xmax_mb * 1e6,
  ]
  
  Fst$POS_Mb <- Fst$POS / 1e6
  Fst$is_high <- Fst$FST > cutoff
  
  #y_lim <- ceiling(max(Fst$FST, na.rm = TRUE) * 100) / 100
  y_lim <- 1 # set the same y lim， the maximum is 1 among all 3 pairs: sapply(data[, c("K_VS_D", "G_VS_D", "C_VS_D")], max, na.rm = TRUE)
  
  
  p <- ggplot(Fst, aes(x = POS_Mb, y = FST)) +
    
    ## all points
    geom_point(
      colour = chr_cols[[chr]],
      size = 2,
      alpha = 0.6
    ) +
    
    ## high FST points
    geom_point(
      data = subset(Fst, is_high),
      colour = "black",
      #colour = chr_cols[[chr]],
      size = 2.8,
      alpha = 0.7
    ) +
    
    ## top 1% line
    geom_hline(
      yintercept = cutoff,
      linetype = "dashed",
      linewidth = 0.8,
      colour = "black"
    ) +
    
    scale_x_continuous(
      limits = c(xmin_mb, xmax_mb),
      breaks = seq(xmin_mb, xmax_mb, by = step_mb),
      expand = expansion(add = c(0.4, 0.4)) ## make the start and end of x-axis a bit expanded
    ) +
    
    # scale_y_continuous(
    #   limits = c(0, y_lim),
    #   expand = expansion(mult = c(0, 0.2))
    # ) +
    
    scale_y_continuous(  # y axis setting
      limits = c(0, y_lim),   
      breaks = seq(0, y_lim, by = 0.25),   
      expand = expansion(mult = c(0, 0.1))
    ) +
    
    labs(
      title = chr,
      x = "Position in chromosome (Mbp)",
      y = expression(F[ST])
    ) +
    
    theme_bw(
      base_size = 10,
      base_family = "Arial"
    ) +
    
    theme(
      panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      
      plot.title = element_text(
        size = 25,
        face = "bold",
        hjust = 0.5
      ),
      
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      axis.text.x  = element_text(size = 23.2),
      axis.text.y  = element_text(size = 24),
      axis.ticks   = element_line(linewidth = 1),
      
      plot.margin = margin(8, 12, 8, 12)
    )
  
  return(p)
}

########## 9 separate plots ##########
## 3 pairs * 3 chromosomes
# selected peak region... of each chr

## K vs D
p_K_chr4 <- make_single_chr_plot("K_VS_D", "Chr4", 38, 52)
p_K_chr5 <- make_single_chr_plot("K_VS_D", "Chr5", 35, 47)
p_K_chrX <- make_single_chr_plot("K_VS_D", "ChrX",  8, 37)


## G vs D
p_G_chr4 <- make_single_chr_plot("G_VS_D", "Chr4", 41, 52)
p_G_chr5 <- make_single_chr_plot("G_VS_D", "Chr5", 36, 40)
p_G_chrX <- make_single_chr_plot("G_VS_D", "ChrX",  8, 37)


## C vs D
p_C_chr4 <- make_single_chr_plot("C_VS_D", "Chr4", 43, 52)
p_C_chr5 <- make_single_chr_plot("C_VS_D", "Chr5", 34, 47)
p_C_chrX <- make_single_chr_plot("C_VS_D", "ChrX",  8, 37)

# print out plots
p_K_chr4
p_K_chr5
p_K_chrX

p_G_chr4
p_G_chr5
p_G_chrX

p_C_chr4
p_C_chr5
p_C_chrX


## ggasave 9 plots

plots <- list(
  K_VS_D_Chr4 = p_K_chr4,
  K_VS_D_Chr5 = p_K_chr5,
  K_VS_D_ChrX = p_K_chrX,
  
  G_VS_D_Chr4 = p_G_chr4,
  G_VS_D_Chr5 = p_G_chr5,
  G_VS_D_ChrX = p_G_chrX,
  
  C_VS_D_Chr4 = p_C_chr4,
  C_VS_D_Chr5 = p_C_chr5,
  C_VS_D_ChrX = p_C_chrX
)

for (nm in names(plots)) {
  ggsave(
    filename = file.path(outdir, paste0(nm, "_peak.pdf")),
    plot = plots[[nm]],
    width = 7.5,
    height = 3.5,
    units = "in",
    device = cairo_pdf
  )
}

####### In the final paper: saved the final plots separately, then merged by Affinity
```