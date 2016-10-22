#------------------------------------------------------------------------
# Title: Consensus of Modules for Individual Networks
# Date: June 1 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# CONTENTS:
  # Choose the soft-thresholding power
  # Run network (functions)
  # Format cluster data into data.frame
  # Save data

# NOTES:
#------------------------------------------------------------------------


# Input Command Line Args
#------------------------------------------------------------------------

library(optparse)

# set arguments
option_list <- list(
  make_option(c("--tissue"), help = "What tissue to analyze")
)

# process arguments
opt <- parse_args(OptionParser(option_list = option_list))

# set input parameters
tissue <- opt$tissue


# Load Libraries
#------------------------------------------------------------------------
library(reshape2)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(WGCNA)


# Load Data
#------------------------------------------------------------------------

# set wd
setwd("E:/MICROARRAY DATA/INDIVIDUAL_NETWORKS/")

# obtain file names
files <- str_subset( list.files(tissue, full.names = TRUE), "networkConstruction")

# open files and assign names
for(f in files){
  load(f)
  assign(str_extract(f, "(chow|highfat)_(female|male)"), network_data)
} 


# Load Data
#------------------------------------------------------------------------

map_overlap <- function(set1, set2, title){
  
  # obtain module color for each gene
  module_colors1 <- set1$moduleColors
  module_colors2 <- set2$moduleColors
  
  # isolate module labels 
  module_labels1 <- unique(module_colors1)
  module_labels2 <- unique(module_colors2)
  
  # number of modules
  nmod1 <- length(module_labels1)
  nmod2 <- length(module_labels2)
  
  # initialize table of counts and p-values
  p_table <- matrix(0, nrow = nmod1, ncol = nmod2)
  count_table <- matrix(0, nrow = nmod1, ncol = nmod2)
  
  # pairwise comparisons of modules between two sets
  for(mod1 in 1:nmod1){
    for(mod2 in 1:nmod2){
      
      # count number that matched to module
      mod1_members <- module_colors1 == module_labels1[mod1]
      mod2_members <- module_colors2 == module_labels2[mod2]
      
      # test whether there are significant differences between two counts
      p_table[mod1, mod2] <- -log10( fisher.test(mod1_members, mod2_members, alternative = "greater")$p.value )
      
      # table of counts belonging to cell
      count_table[mod1, mod2] <- sum(mod1_members & mod2_members)
      
    }
  }
  
  # truncate p values smaller than 10^{-50} to 10^{-50}
  p_table[is.infinite(p_table)] = 1.3*max(p_table[is.finite(p_table)]);
  p_table[p_table>50 ] = 50 ;
  
  # Marginal counts (really module sizes)
  set1_totals <- rowSums(count_table)
  set2_totals <- colSums(count_table)
  
  
  
  png( paste0(tissue, "/preliminary consensus ", title, ".png"), height = 1000, width = 1200)
  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  print(
    labeledHeatmap(Matrix = p_table,
                 xLabels = paste(" ", module_labels2),
                 yLabels = paste(" ", module_labels1),
                 colorLabels = TRUE,
                 # xSymbols = paste("Set2 ", module_labels2, ": ", nmod2, sep=""),
                 # ySymbols = paste("Set1 ", module_labels1, ": ", nmod1, sep=""),
                 # textMatrix = count_table,
                 colors = greenWhiteRed(100)[50:100],
                 main = title,
                 cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE
                 )
  )
  dev.off()
  
  
  return()
  
}


map_overlap(chow_male, highfat_male, "chow male vs highfat male")
map_overlap(highfat_male, highfat_female, "highfat male vs highfat female")


