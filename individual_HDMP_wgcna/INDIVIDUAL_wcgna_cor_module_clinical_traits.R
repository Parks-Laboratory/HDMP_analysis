#------------------------------------------------------------------------
# Title: Correlate Eigengenes with Clinical Traits
# Date: June 7 2016
# Author: Jenny Nguyen
# Email: jnnguyen2@wisc.edu
#------------------------------------------------------------------------

# CONTENTS:


# NOTES:
  # files are big - need to multiple pages of plot output 
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
library(RODBC)
library(WGCNA)


# Load Data
#------------------------------------------------------------------------

source("E:/MICROARRAY DATA/wgcna_pipeline/wgcna_network_pipeline.R")
setwd( paste0("E:/MICROARRAY DATA/INDIVIDUAL_NETWORKS/", tissue) )


# Load Data: Clinical Traits
#------------------------------------------------------------------------

db <- odbcDriverConnect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')

# clinical trait data from DB
data_sources <- sqlTables(db) %>% 
  subset(TABLE_SCHEM == "dbo") %>% 
  subset( str_detect(TABLE_NAME, "avg_clinical_traits_by_strain") )
queries <- paste0("select * from dbo.", data_sources$TABLE_NAME)
clinical_traits <- lapply(queries, function(q) sqlQuery(db, q))

# filter clinical trait columns using columns in common
col_to_keep <- sqlQuery(db, "select * from dbo.clinical_traits_in_common where highfat != '' and notes != 'ignore'") %>% 
  mutate_each(funs(as.character))

# extract columns in common
clinical_traits <- lapply(clinical_traits, function(d){
  
  # find the appropriate columns to keep (ie is it chow or highfat?)
  keep_var1 <- Filter(function(x) x %in% colnames(d), col_to_keep$chow)
  keep_var2 <- Filter(function(x) x %in% colnames(d), col_to_keep$highfat)
  
  # keep the ones with the largest number of matches
  if(length(keep_var1) > length(keep_var2)) keep_var <- keep_var1 else keep_var <- keep_var2
  rd <- dplyr::select(d, strain, one_of(keep_var))
  
  # rename columns to common name
  sub_col_to_keep <- subset(col_to_keep, chow %in% keep_var | highfat %in% keep_var)
  try_default( default = NULL, setnames(rd, sub_col_to_keep$chow, sub_col_to_keep$common_name), quiet = TRUE )
  try_default( default = NULL, setnames(rd, sub_col_to_keep$highfat, sub_col_to_keep$common_name), quiet = TRUE )
  
  # return data
  return(rd)
})

# assign names to clinical trait data
names(clinical_traits) <- str_extract(queries, "(chow|highfat)_(male|female)")

odbcClose(db)


# Load Data: Strain Names & Mod Eigengenes
#------------------------------------------------------------------------

# load the module eigengenes assign strain name
options <- c("chow_male", "highfat_male", "highfat_female")
network <- lapply(options, function(x) fread(paste0(x, "_module_eigengenes.txt")))
names(list) <- options


# Process Data
#------------------------------------------------------------------------

cor_data <- lapply(options, function(opt){
  # subset clinical traits
  CTs <- subset(clinical_traits[[opt]], strain %in% network[[opt]]$strain)
  MEs <- subset(network[[opt]], strain %in% CTs$strain)
  
  # return as list
  return(list(clinical_trait = CTs, mod_eigengene = MEs))
})
names(cor_data) <- options


# Run Correlation Heatmaps
#------------------------------------------------------------------------

# correlate clinical traits
l_ply(options, function(opt){
  
  # print out plot - split up plots by 25 modules per plot
  pdf( paste0("clinical_trait_cor_", tissue, "_", opt, ".pdf"), height = 8, width = 14)
  
  # print in iterations of 25 (last column is strain)
  loop_vars <- c( seq(1, ncol(cor_data[[opt]]$mod_eigengene), by = 14), ncol(cor_data[[opt]]$mod_eigengene) )
  for(i in 1:(length(loop_vars)-1)){
    keep_col <- c(loop_vars[i]:(loop_vars[i+1]-1))
    print( with(cor_data[[opt]], clinical_trait_cor(mod_eigengene[, c(keep_col, ncol(mod_eigengene))], clinical_trait)) )
  }
  
  dev.off()
})


# 
# # colors corresponding to cholesterol
# # chow male: black
# # highfat male: lightyellow
# # highfat female: greenyellow
# # reactive: black
# 
# # chow-male
# ct <- cor_data$chow_male$clinical_trait
# me <- cor_data$chow_male$mod_eigengene %>% dplyr::select(strain, MEblack) %>% mutate(me = MEblack)
# plot_data1 <- merge(ct, me, "strain") %>% mutate(kind = "chow_male")
# 
# # highfat male
# ct <- cor_data$highfat_male$clinical_trait
# me <- cor_data$highfat_male$mod_eigengene %>% dplyr::select(strain, MElightyellow) %>% mutate(me = MElightyellow)
# plot_data2 <- merge(ct, me, "strain") %>% mutate(kind = "highfat_m")
# 
# # highfat female
# ct <- cor_data$highfat_female$clinical_trait
# me <- cor_data$highfat_female$mod_eigengene %>% dplyr::select(strain, MEgreenyellow) %>% mutate(me = MEgreenyellow)
# plot_data3 <- merge(ct, me, "strain") %>% mutate(kind = "highfat_f")
# 
# # combine
# plot_data <- as.data.frame(rbindlist(list(plot_data1, plot_data2, plot_data3), fill = TRUE))
# 
# 
# # for(i in colnames(ct)[-1]){
# pdf("test.pdf", height = 7, width = 7)
# l <- lapply(colnames(ct)[-1], function(i){
#   plot_data$plot_col <- plot_data[, i]
#   print( ggplot(data = plot_data, aes(me, plot_col, color = kind)) + geom_smooth(method = "lm") + ylab(i) )
# })
# dev.off()
# 
